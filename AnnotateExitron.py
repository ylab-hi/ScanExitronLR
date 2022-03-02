#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 29 12:01:00 2022

@author: Josh Fry @YangLab UMN
"""
__version__ = 'v0.1'
import os
import sys
import pybedtools
import gffutils
import pysam
import argparse
import subprocess
import traceback
import pandas as pd
from Bio.Seq import Seq
from shutil import rmtree


#===============================================================================
# Helper Methods
#===============================================================================

def parse_args():
    parser = argparse.ArgumentParser(
        description="%(prog)s",
        epilog="ScanExitronLR v0.1: detecting exitron splicing events using RNA-Seq data",
    )
    parser.add_argument(
        "-i",
        "--input",
        action="store",
        dest="input",
        help="Input exitron file.",
        required=True,
    )
    parser.add_argument(
        "-g",
        "--genome",
        action="store",
        dest="genome_ref",
        help="Path to reference genome.",
        default=None,
        required=True
    )
    parser.add_argument(
        "-r",
        "--reference-annotations",
        action="store",
        dest="annotation_ref",
        help="Path to reference annotations.",
        default=None,
        required=True
    )
    parser.add_argument(
        "-o",
        "--output",
        action="store",
        dest="out",
        help="Output filename. (default: $input_filename_annotated.exitron)",
        default=None,
    )
    parser.add_argument(
        "-b",
        "--bam-file",
        action="store",
        dest="bam_file",
        help="Input bamfile. (optional)",
        default=None,
        required=True
    )
    parser.add_argument(
        "-pd",
        "--prot-domains",
        action="store",
        dest="prot_domains",
        help="Input tab separated file of pfam protein domains. (optional)",
        default=os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data/prot_domains_pfam.tsv'),
    )
    parser.add_argument(
        "-v", "--version", action="version", version="%(prog)s {}".format(__version__)
    )
    args = parser.parse_args()
    return args



def read_exitron_file(filename):
    """


    Parameters
    ----------
    filename : str
        path to exitron file.

    Returns
    -------
    List of dictionaries, where each dict is an exitron record.

    """
    exitrons = []
    with open(filename, 'r') as f:
        header = f.readline().rstrip().split('\t')
        for line in f:
            exitron = {}
            entries = line.rstrip().split('\t')
            for col, ent in zip(header, entries):
                exitron[col] = ent
            exitrons.append(exitron)
    return exitrons


#=============================================================================
# Modules
#=============================================================================



def get_nmd_status(exitron, frameshift_pos, bamfile):
    e_reads = exitron['reads'].split(',')
    nmd = 0
    tot = 0

    match_or_deletion = {0, 2, 7, 8} # only M/=/X (0/7/8) and D (2) are related to genome position
    BAM_CREF_SKIP = 3

    for read in bamfile.fetch(exitron['chrom'], int(exitron['start']), int(exitron['end'])):
        ej_pos = []
        # bamfiles are indexed only by position, not name, so we just have to
        # find our spliced reads by position
        if read.query_name in e_reads:
            base_position = read.pos
            read_position = 0
            # iterate through cigar string looking for N
            for i, (tag, nt) in enumerate(read.cigartuples):
                # if (0, X), keep track of base_position.
                # if (3, X), which corresponds to N,
                # look at match before and after
                if tag in match_or_deletion:
                    base_position += nt
                    read_position += nt
                elif read.cigartuples[i][0] == BAM_CREF_SKIP:
                    junc_start = base_position
                    base_position += nt
                    ej_pos.append((junc_start, base_position))
            if exitron['strand'] == '+':
                nmd += any(x > frameshift_pos + 50 for x in map(lambda x: x[0], ej_pos))
                tot += 1

            else:
                nmd += any(x < frameshift_pos - 50 for x in map(lambda x: x[0], ej_pos))
                tot += 1


    return (nmd, tot)


def get_gene_exitron_seq(exitron, db, genome_fn):
    seq = ''
    e_start = int(exitron['start'])
    e_end = int(exitron['end'])
    strand = exitron['strand']
    transcript = exitron['transcript_id']
    try:
        start_codon = next(db.children(db[transcript], featuretype = ['start_codon']))
        stop_codon = next(db.children(db[transcript], featuretype = ['stop_codon']))
    except StopIteration: # this occurs when the annotation does not provide a start/stop codon
        # see tags 'mRNA_start_NF', 'cds_start_NF' etc
        return None, None, None, None, None, None, None
    ej_after_start_codon = []
    ej_before_start_codon = []
    after_start_codon = False
    before_start_codon = True
    seq_pos = []

    for exon in db.children(db[transcript], featuretype = 'exon', order_by = 'start'):
        if exon.start <= start_codon.start <= exon.end:
            start_codon_pos = len(seq) + start_codon.start - exon.start
            after_start_codon = True
            before_start_codon = False
        if exon.start <= stop_codon.start <= exon.end:
            stop_codon_pos = len(seq) + stop_codon.start - exon.start
        if exon.start < e_start < exon.end and exon.start < e_end < exon.end:
            cds_seq = exon.sequence(genome_fn).upper() if strand == '+' else str(Seq(exon.sequence(genome_fn).upper()).reverse_complement())
            try:
                exitron_pos = len(seq) + e_start - exon.start - start_codon_pos if strand == '+' else len(seq) + e_start - exon.start - stop_codon_pos
            except UnboundLocalError:
                # occurs when neither stop nor start codon was encountered before exitron
                # this can happen if the inferred transcript is one where the
                # exitron is located in a UTR region instead of CDS
                return None, None, None, None, None, None, None
            seq += cds_seq[:e_start - exon.start + 1] + cds_seq[e_end - exon.start:]
            seq_pos.extend([*range(exon.start, e_start + 1)])
            seq_pos.extend([*range(e_end, exon.end + 1)])

            # if we encountered the start or stop codon in the same exon as the
            # exitron, the start and stop codons will be positioned relative to
            # unspliced sequence.  thus we need to subtract length of exitron
            if exon.start <= start_codon.start <= exon.end and strand == '-':
                start_codon_pos -= int(exitron['length'])
            if exon.start <= stop_codon.start <= exon.end and strand == '+':
                stop_codon_pos -= int(exitron['length'])

        else:
            # if strand is -, exon.sequence returns reverse complement
            # however, it's easier to build forward strand exon by exon and
            # only afterwards reverse complement to get gene/prot sequence.
            seq += exon.sequence(genome_fn).upper() if strand == '+' else str(Seq(exon.sequence(genome_fn).upper()).reverse_complement())
            seq_pos.extend([*range(exon.start, exon.end + 1)])
            if after_start_codon: ej_after_start_codon.append(len(seq) - 1)
            if before_start_codon: ej_before_start_codon.append(len(seq) - 1)

    return seq, start_codon_pos, stop_codon_pos, exitron_pos, seq_pos, ej_after_start_codon, ej_before_start_codon

def get_gene_seq(exitron, db, genome_fn):
    seq = ''
    transcript = exitron['transcript_id']

    for cds in db.children(db[transcript], featuretype = 'CDS', order_by = 'start'):
        seq += cds.sequence(genome_fn).upper() if exitron['strand'] == '+' else str(Seq(cds.sequence(genome_fn).upper()).reverse_complement())
    return seq

def categorize_exitron(exitron, transcript, bamfile, db, genome_fn):
    """


    Parameters
    ----------
    exitron : exitron

    Returns
    -------
    categorizes exitron as frameshift, truncation, truncation + missense.

    if frameshift, returns framshifted protein and whether it triggers NMD
    if truncated + missense, returns the missense mutation plus protein sequence
    if just truncated, returns protein sequence

    (description, detail, dna_seq, prot_seq)
    """
    exitron['transcript_id'] = transcript
    seq, start_codon_pos, stop_codon_pos, exitron_pos, seq_pos, ej_after_start_codon, ej_before_start_codon = get_gene_exitron_seq(exitron, db, genome_fn)
    if not seq:
        # without start/stop codons, all we can infer is truncation vs frameshift
        exitron['exitron_prot_position'] = 'NA'
        exitron['type'] = 'truncated' if int(exitron['length']) % 3 == 0 else 'frameshift'
        exitron['substitution'] = '.'
        exitron['nmd_status_predicted'] = 'NA'
        exitron['nmd_status_read_percentage'] = 'NA'
        exitron['downstream_inframe_AUG'] = 'NA'
        exitron['start_proximal_PTC'] = 'NA'
        return exitron, None, None

    seq = Seq(seq) if exitron['strand'] == '+' else Seq(seq).reverse_complement()
    if exitron['strand'] == '-':
        # recalibrate position
        # we keep track of original start_codon_pos in backward seq in order to
        # test exon-exon borders
        orig_start_codon_pos = start_codon_pos
        start_codon_pos = len(seq) - start_codon_pos - 3

        stop_codon_pos = len(seq) - stop_codon_pos - 3
        exitron_pos = len(seq[start_codon_pos:stop_codon_pos+3]) - exitron_pos
        seq_pos = seq_pos[::-1]

    # exitron is frameshift iff length % 3 != 0
    if int(exitron['length']) % 3 != 0:

        frameshift_prot = seq[start_codon_pos:].translate(to_stop = True)
        frameshift_pos = seq_pos[start_codon_pos + len(frameshift_prot)*3 - 1]
        nmd, tot = get_nmd_status(exitron, frameshift_pos,  bamfile)

        dna_seq = seq[start_codon_pos:start_codon_pos + len(frameshift_prot)*3 + 3]
        # use ej_after_start_codon[:-1] here because the last exon boundry does not count as an EJC
        # TODO: make sure this works with '-' strand...


        if exitron['strand'] == '+':
            nmd_pred = 'NMD' if any(x > len(frameshift_prot)*3 + 50 for x in ej_after_start_codon[:-1]) else 'no_NMD'
        else:
            nmd_pred = 'NMD' if any(x < orig_start_codon_pos - len(frameshift_prot)*3 - 50 for x in ej_before_start_codon) else 'no_NMD'

        exitron['exitron_prot_position'] = exitron_pos//3
        exitron['type'] = 'frameshift'
        exitron['substitution'] = '.'
        exitron['nmd_status_predicted'] = nmd_pred
        exitron['nmd_status_read_percentage'] = nmd/tot
        exitron['downstream_inframe_AUG'] = 'M' in seq[start_codon_pos + len(frameshift_prot)*3 + 3:].translate()
        exitron['start_proximal_PTC'] = exitron_pos <= 200

        return exitron, str(dna_seq), str(frameshift_prot) + '*'

    # exitron is truncated + missense iff (exitron_pos - 1) % 3 != 0
    elif (exitron_pos - 1) % 3 != 0:
        seq_full = get_gene_seq(exitron, db, genome_fn) if exitron['strand'] == '+' else str(Seq(get_gene_seq(exitron, db, genome_fn)).reverse_complement())
        # new aa is at aa where exitron begins in spliced protein sequence
        new_aa = str(seq[start_codon_pos:stop_codon_pos + 3].translate())[exitron_pos//3]
        # old aa is at aa where exitron begins in non-spliced protein sequence
        old_aa = str(Seq(seq_full).translate())[exitron_pos//3]

        exitron['exitron_prot_position'] = exitron_pos//3
        exitron['type'] = 'truncated+missense'
        exitron['substitution'] = f'{old_aa}->{new_aa}'
        exitron['nmd_status_predicted'] = '.'
        exitron['nmd_status_read_percentage'] = '.'
        exitron['downstream_inframe_AUG'] = '.'
        exitron['start_proximal_PTC'] = '.'

        return (exitron,
                str(seq[start_codon_pos:stop_codon_pos+3]),
                str(seq[start_codon_pos:stop_codon_pos+3].translate()))
    # otherwise, exitron is truncation
    else:

        exitron['exitron_prot_position'] = exitron_pos//3
        exitron['type'] = 'truncated'
        exitron['substitution'] = '.'
        exitron['nmd_status_predicted'] = '.'
        exitron['nmd_status_read_percentage'] = '.'
        exitron['downstream_inframe_AUG'] = '.'
        exitron['start_proximal_PTC'] = '.'

        return (exitron,
                str(seq[start_codon_pos:stop_codon_pos+3]),
                str(seq[start_codon_pos:stop_codon_pos+3].translate()))



# seq, start_codon_pos, stop_codon_pos, exitron_pos, after_start_codon = get_gene_exitron_seq(exitron)
# seq_full = get_gene_seq(exitron)
# bamfile = pysam.AlignmentFile('data/test_mt111421_pb.bam')

# for exitron in read_exitron_file('data/test.exitron'):
#     print(categorize_exitron(exitron, bamfile))

def get_pfam_domains(exitron, prot_df):
    df_gene = prot_df[prot_df['Transcript stable ID version'] == exitron['transcript_id']]
    if exitron['type'][:9] == 'truncated':
        e_start = int(exitron['exitron_prot_position'])
        e_end = e_start + int(exitron['length'])//3
        pf_ids = df_gene[(((df_gene['Pfam start'] >= e_start) & (df_gene['Pfam start'] <= e_end)) |
                          ((df_gene['Pfam end'] >= e_start) & (df_gene['Pfam end'] <= e_end)) |
                          ((df_gene['Pfam start'] <= e_start) & (df_gene['Pfam end'] >= e_end)))]
        pf_ids = set(pf_ids['Pfam ID'])
        if pf_ids:
            exitron['prot_domains'] = ','.join(pf_ids)
        else:
            exitron['prot_domains'] = '.'
    elif exitron['type'] == 'frameshift':
        e_start = int(exitron['exitron_prot_position'])
        pf_ids = df_gene[(df_gene['Pfam start'] > e_start)]
        pf_ids = set(pf_ids['Pfam ID'])
        if pf_ids:
            exitron['prot_domains'] = ','.join(pf_ids)
        else:
            exitron['prot_domains'] = '.'


#=============================================================================
# Main
#=============================================================================


def main(tmp_path):
    # Get arguments
    args = parse_args()

    # Check to see if bamfile can be opened and there is an index.
    if args.bam_file:
        try:
            bamfile = pysam.AlignmentFile(args.bam_file, 'rb', require_index = True)
        except FileNotFoundError:
            try:
                print('Building bam index file')
                pysam.index(args.input)
                bamfile = pysam.AlignmentFile(args.input, 'rb', require_index = True)
            except FileNotFoundError:
                print(f'There is a problem opening bam file at: {args.input}')
    try:
        db = gffutils.create_db(args.annotation_ref,
                                dbfn = args.annotation_ref + '.db',
                                disable_infer_genes = True,
                                disable_infer_transcripts = True)
        db = gffutils.FeatureDB(args.annotation_ref + '.db')
    except:
        print(f'Using annotation databse {args.annotation_ref + ".db"}')
        db = gffutils.FeatureDB(args.annotation_ref + '.db')

    if not args.out:
        args.out = args.input + '.annotated'

    with open(args.out, 'w') as out:
        header = ['chrom',
                  'start',
                  'end',
                  'name',
                  'region',
                  'ao',
                  'strand',
                  'gene_name',
                  'gene_id',
                  'length',
                  'splice_site',
                  'transcript_id',
                  'pso',
                  'dp',
                  'consensus_prop',
                  'exitron_prot_position',
                  'type',
                  'substitution',
                  'nmd_status_predicted',
                  'nmd_status_read_percentage',
                  'downstream_inframe_AUG',
                  'start_proximal_PTC',
                  'prot_domains',
                  'reads']
        for column in header:
            out.write(column + '\t')
        out.write('\n')
        prot_df = pd.read_csv(args.prot_domains, delimiter='\t')
        exitrons = read_exitron_file(args.input)
        for exitron in exitrons:
            transcripts = exitron['transcript_id'].split(';')
            for transcript in transcripts:
                t_id = transcript.split(',')[0]
                abundance = transcript.split(',')[1]
                # dna / prot seq output not yet implemented
                res, _, _= categorize_exitron(exitron.copy(), t_id, bamfile, db, args.genome_ref)

                # idetnify pfam domains
                if res['exitron_prot_position'] != 'NA':
                    get_pfam_domains(res, prot_df)
                else:
                    res['prot_domains'] = 'NA'

                # update abundance
                res['transcript_id'] += f',{abundance}'
                for column in header:
                    out.write(str(res[column]) + '\t')
                out.write('\n')


if __name__ == '__main__':
    # Set tmp directory
    this_dir = os.path.dirname(os.path.realpath(__file__))
    tmp_path = os.path.join(this_dir, f'scanexitron_tmp{os.getpid()}')
    try: os.mkdir(tmp_path)
    except FileExistsError: pass
    pybedtools.helpers.set_tempdir(tmp_path)
    try:
        main(tmp_path)
    except Exception as e:
        if e.__class__.__name__ == 'InterruptedError':
            sys.stderr.write("User interrupt!")
        else:
            traceback.print_tb(e.__traceback__)
        pybedtools.helpers.cleanup(remove_all=True)
        rmtree(tmp_path)
        sys.exit(1)






