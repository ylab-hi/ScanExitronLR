#!/usr/bin/env python3
#
# @author: Josh Fry @YangLab, Hormel Institute, UMN
#
# ===============================================================================

__version__ = 'v1.1.6'
import os
import sys
import gffutils
import pysam
import argparse
import traceback
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


# ===============================================================================
# Helper Methods
# ===============================================================================

def parse_args():
    usage = """
    \t\t-i STR\t\tREQUIRED: Input exitron file, generated from selr extract
    \t\t-g STR\t\tREQUIRED: Input genome reference (e.g. hg38.fa)
    \t\t-r STR\t\tREQUIRED: Input *sorted* and *bgzip'd* annotation reference (e.g. gencode_v38_sorted.gtf.gz).
    \t\t\t\tIf your annotation file is not sorted, use the following command: awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k4,4n -k5,5n"}' in.gtf > out_sorted.gtf
    \t\t\t\tIf your annotation file is not gziped, use the following command: bgzip -c in.gtf > out.gtf.gz
    \t\t-o STR\t\tOutput filename (e.g. bam_filename.exitron.annotation <- this is default)
    \t\t-b/--bam-file STR\t\tIf specified, annotation includes read supported NMD status directly from alignments.
    \t\t-arabidopsis\tUse this flag if using alignments from Arabidopsis. Interrupted Arabidopsis protein domains will be reported.
    """
    parser = argparse.ArgumentParser(
        description="annotate",
        usage=usage
    )
    parser.add_argument(
        "-i",
        "--input",
        action="store",
        dest="input",
        required=True,
    )
    parser.add_argument(
        "-g",
        "--genome",
        action="store",
        dest="genome_ref",
        default=None,
        required=True
    )
    parser.add_argument(
        "-r",
        "--reference-annotations",
        action="store",
        dest="annotation_ref",
        default=None,
        required=True
    )
    parser.add_argument(
        "-o",
        "--output",
        action="store",
        dest="out",
        default=None,
    )
    parser.add_argument(
        "-b",
        "--bam-file",
        action="store",
        dest="bam_file",
        default=None,
    )
    parser.add_argument(
        "-arabidopsis",
        "--arabidopsis",
        action="store_true",
        dest="arabidopsis",
    )
    parser.add_argument(
        "-fasta",
        "--fasta",
        action="store_true",
        dest="fasta",
    )
    parser.add_argument(
        "-v", "--version", action="version", version=f"%(prog)s {__version__}"
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
    with open(filename) as f:
        header = f.readline().rstrip().split('\t')
        for line in f:
            exitron = {}
            entries = line.rstrip().split('\t')
            for col, ent in zip(header, entries):
                exitron[col] = ent
            exitrons.append(exitron)
    return exitrons


# =============================================================================
# Modules
# =============================================================================


def get_nmd_status(exitron, frameshift_pos, bamfile):
    """


    Parameters
    ----------
    exitron : dict
        single exitron dict.
    frameshift_pos : int
        DESCRIPTION.
    bamfile : pysam.AlignmentFile

    Returns
    -------
    nmd : int
        number of reads supporting NMD
    tot : int
        number of reads total with exitrons

    """
    e_reads = exitron['reads'].split(',')
    nmd = 0
    tot = 0

    # only M/=/X (0/7/8) and D (2) are related to genome position
    match_or_deletion = {0, 2, 7, 8}
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
                nmd += any(x > frameshift_pos +
                           50 for x in map(lambda x: x[0], ej_pos))
                tot += 1

            else:
                nmd += any(x < frameshift_pos -
                           50 for x in map(lambda x: x[0], ej_pos))
                tot += 1

    return nmd, tot


def get_gene_exitron_seq(exitron, db, genome_fn, arabidopsis):
    """


    Parameters
    ----------
    exitron : dict
    db : gffutils database
    genome_fn : str
    arabidopsis : bool

    Returns
    -------
    seq : str
        complete gene sequence (including UTRs)
    start_codon_pos: int
        start_codon_pos *within* seq
    stop_codon_pos: int
        stop_codon_pos *within* seq
    exitron_pos: int
        exitron_pos *within* seq
    seq_pos: list
        this is a list of genome positions for each nt in seq
    ej_after_start_codon: list of int
        exon-exon junction locations after start_codon, used for NMD identification
    ej_before_start_codon: list of int
        exon-exon junction locations before start_codon, used for NMD identification

    """
    seq = ''
    e_start = int(exitron['start'])
    e_end = int(exitron['end'])
    strand = exitron['strand']
    transcript = exitron['transcript_id']
    exitron_pos = None
    try:
        if not arabidopsis:
            try:
                start_codon = next(db.children(
                    db[transcript], featuretype=['start_codon']))
                stop_codon = next(db.children(
                    db[transcript], featuretype=['stop_codon']))
                start_codon_start = start_codon.start
                stop_codon_start = stop_codon.start
            except StopIteration:  # this occurs when the annotation does not provide a start/stop codon
                # see tags 'mRNA_start_NF', 'cds_start_NF' etc
                # this is why we do not infer start and stop codons in human annotation
                # some transcripts explicitly leave them out to indicate that the transcript
                # is partial. If the transcript is partial, just report frameshift/trunc
                return None, None, None, None, None, None, None
        else:
            if strand == '+':
                start_codon_start = next(db.children(
                    db[transcript], featuretype=['CDS'])).start
                stop_codon_start = next(db.children(db[transcript], featuretype=[
                                        'CDS'], order_by='start', reverse=True)).start
            else:
                start_codon_start = next(db.children(db[transcript], featuretype=[
                                         'CDS'], order_by='start', reverse=True)).end - 2
                stop_codon_start = next(db.children(
                    db[transcript], featuretype=['CDS'])).start
    except:
        return None, None, None, None, None, None, None
    ej_after_start_codon = []
    ej_before_start_codon = []
    after_start_codon = False
    before_start_codon = True
    seq_pos = []

    for exon in db.children(db[transcript], featuretype='exon', order_by='start'):
        if exon.start <= start_codon_start <= exon.end:
            start_codon_pos = len(seq) + start_codon_start - exon.start
            after_start_codon = True
            before_start_codon = False
        if exon.start <= stop_codon_start <= exon.end:
            stop_codon_pos = len(seq) + stop_codon_start - exon.start
        if exon.start < e_start < exon.end and exon.start < e_end < exon.end:
            cds_seq = exon.sequence(genome_fn).upper(
            ) if strand == '+' else str(Seq(exon.sequence(genome_fn).upper()).reverse_complement())
            exitron_pos = len(seq) + e_start - exon.start - start_codon_pos if strand == '+' else len(
                seq) + e_start - exon.start - stop_codon_pos
            seq += cds_seq[:e_start - exon.start + 1] + \
                cds_seq[e_end - exon.start:]
            seq_pos.extend([*range(exon.start, e_start + 1)])
            seq_pos.extend([*range(e_end, exon.end + 1)])

            # if we encountered the start or stop codon in the same exon as the
            # exitron, the start and stop codons will be positioned relative to
            # unspliced sequence.  thus we need to subtract length of exitron
            if exon.start <= start_codon_start <= exon.end and strand == '-':
                start_codon_pos -= int(exitron['length'])
            if exon.start <= stop_codon_start <= exon.end and strand == '+':
                stop_codon_pos -= int(exitron['length'])

        else:
            # if strand is -, exon.sequence returns reverse complement
            # however, it's easier to build forward strand exon by exon and
            # only afterwards reverse complement to get gene/prot sequence.
            seq += exon.sequence(genome_fn).upper() if strand == '+' else str(
                Seq(exon.sequence(genome_fn).upper()).reverse_complement())
            seq_pos.extend([*range(exon.start, exon.end + 1)])
            if after_start_codon:
                ej_after_start_codon.append(len(seq) - 1)
            if before_start_codon:
                ej_before_start_codon.append(len(seq) - 1)

    if not exitron_pos:
        # This can happen if LIQA detects exitron in a transcript for which the
        # start or end point of the exitron is not present.  This happens if LIQA
        # just gets the transcript wrong, or if the exitron is similar to an intron
        # in another transcript--LIQA may assign some abundance to that transcript
        # If this happens, just return None and move on.
        return None, None, None, None, None, None, None
    return seq, start_codon_pos, stop_codon_pos, exitron_pos, seq_pos, ej_after_start_codon, ej_before_start_codon


def get_gene_seq(exitron, db, genome_fn):
    """


    Parameters
    ----------
    exitron : dict
    db : gffutils database
    genome_fn : str

    Returns
    -------
    seq : str
        full *coding sequence* of exitron spliced gene

    """
    seq = ''
    transcript = exitron['transcript_id']

    for cds in db.children(db[transcript], featuretype='CDS', order_by='start'):
        seq += cds.sequence(genome_fn).upper() if exitron['strand'] == '+' else str(
            Seq(cds.sequence(genome_fn).upper()).reverse_complement())
    return seq


def categorize_exitron(exitron, transcript, bamfile, db, genome_fn, arabidopsis):
    """


    Parameters
    ----------
    exitron : exitron

    Returns
    -------
    categorizes exitron as frameshift, truncation, truncation + substitution.

    if frameshift, returns framshifted protein and whether it triggers NMD
    if truncated + substitution, returns the missense mutation plus protein sequence
    if just truncated, returns protein sequence

    (description, detail, dna_seq, prot_seq)
    """
    exitron['transcript_id'] = transcript
    seq, start_codon_pos, stop_codon_pos, exitron_pos, seq_pos, ej_after_start_codon, ej_before_start_codon = get_gene_exitron_seq(
        exitron, db, genome_fn, arabidopsis)
    if not seq:
        # without start/stop codons, all we can infer is truncation vs frameshift
        exitron['exitron_prot_position'] = '.'
        exitron['type'] = 'truncated' if int(
            exitron['length']) % 3 == 0 else 'frameshift'
        exitron['substitution'] = '.'
        exitron['nmd_status_predicted'] = '.'
        exitron['nmd_status_read_percentage'] = '.'
        exitron['downstream_inframe_AUG'] = '.'
        exitron['start_proximal_PTC'] = '.'
        return exitron, None, None

    seq = Seq(
        seq) if exitron['strand'] == '+' else Seq(seq).reverse_complement()
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
        # add trailing Ns so that sequence is divisible by 3 (otherwise biopython complains)
        n_adjust = 'N'*(3 - len(seq[start_codon_pos:]) % 3)
        frameshift_seq = seq[start_codon_pos:] + n_adjust
        frameshift_prot = frameshift_seq.translate(to_stop=True)
        frameshift_pos = seq_pos[start_codon_pos + len(frameshift_prot)*3 - 1]
        nmd, tot = get_nmd_status(
            exitron, frameshift_pos,  bamfile) if bamfile else (-1, -1)

        dna_seq = seq[start_codon_pos:start_codon_pos +
                      len(frameshift_prot)*3 + 3]

        # use ej_after_start_codon[:-1] here because the last exon boundry does not count as an EJC
        if exitron['strand'] == '+':
            nmd_pred = 'NMD' if any(x > len(
                frameshift_prot)*3 + 50 for x in ej_after_start_codon[:-1]) else 'no_NMD'
        else:
            nmd_pred = 'NMD' if any(x < orig_start_codon_pos - len(
                frameshift_prot)*3 - 50 for x in ej_before_start_codon) else 'no_NMD'

        exitron['exitron_prot_position'] = exitron_pos//3
        exitron['type'] = 'frameshift'
        exitron['substitution'] = '.'
        exitron['nmd_status_predicted'] = nmd_pred
        exitron['nmd_status_read_percentage'] = nmd/tot if tot > 0 else '.'
        exitron['downstream_inframe_AUG'] = 'M' in (
            seq[start_codon_pos + len(frameshift_prot)*3 + 3:] + n_adjust).translate()
        exitron['start_proximal_PTC'] = exitron_pos <= 200

        return exitron, dna_seq, frameshift_prot + '*'

    # exitron is truncated + substitution iff (exitron_pos - 1) % 3 != 0
    elif (exitron_pos - 1) % 3 != 0:
        seq_full = get_gene_seq(exitron, db, genome_fn) if exitron['strand'] == '+' else str(
            Seq(get_gene_seq(exitron, db, genome_fn)).reverse_complement())
        # new aa is at aa where exitron begins in spliced protein sequence
        try:
            new_aa = str(seq[start_codon_pos:stop_codon_pos +
                         3].translate())[exitron_pos//3]
        except:
            # This can happen if LIQA predicts that exitron occurs in transcript
            # at 3' or 5' UTR.  This can happen if exitron was originally called in
            # CDS region.  If this happens, just continue
            exitron['exitron_prot_position'] = '.'
            exitron['type'] = '.'
            exitron['substitution'] = '.'
            exitron['nmd_status_predicted'] = '.'
            exitron['nmd_status_read_percentage'] = '.'
            exitron['downstream_inframe_AUG'] = '.'
            exitron['start_proximal_PTC'] = '.'
            return exitron, None, None
        # old aa is at aa where exitron begins in non-spliced protein sequence
        old_aa = str(Seq(seq_full).translate())[exitron_pos//3]

        exitron['exitron_prot_position'] = exitron_pos//3
        exitron['type'] = 'truncated+substitution'
        exitron['substitution'] = f'{old_aa}->{new_aa}'
        exitron['nmd_status_predicted'] = '.'
        exitron['nmd_status_read_percentage'] = '.'
        exitron['downstream_inframe_AUG'] = '.'
        exitron['start_proximal_PTC'] = '.'

        return (exitron,
                seq[start_codon_pos:stop_codon_pos+3],
                seq[start_codon_pos:stop_codon_pos+3].translate())
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
                seq[start_codon_pos:stop_codon_pos+3],
                seq[start_codon_pos:stop_codon_pos+3].translate())


def get_pfam_domains(exitron, prot_df):
    """


    Parameters
    ----------
    exitron : dict
    prot_df : pandas dataframe
        columns: transcript, PFAM_domain, start, end

    Returns
    -------
    None. Mutates exitron dict.

    """
    df_gene = prot_df[prot_df['transcript'] == exitron['transcript_id']]
    if exitron['type'][:9] == 'truncated':
        e_start = int(exitron['exitron_prot_position'])
        e_end = e_start + int(exitron['length'])//3
        pf_ids = df_gene[(((df_gene['start'] >= e_start) & (df_gene['start'] <= e_end)) |
                          ((df_gene['end'] >= e_start) & (df_gene['end'] <= e_end)) |
                          ((df_gene['start'] <= e_start) & (df_gene['end'] >= e_end)))]
        pf_ids = set(pf_ids['PFAM_domain'])
        if pf_ids:
            exitron['prot_domains'] = ','.join(pf_ids)
        else:
            exitron['prot_domains'] = '.'
    elif exitron['type'] == 'frameshift':
        e_start = int(exitron['exitron_prot_position'])
        pf_ids = df_gene[(df_gene['start'] > e_start)]
        pf_ids = set(pf_ids['PFAM_domain'])
        if pf_ids:
            exitron['prot_domains'] = ','.join(pf_ids)
        else:
            exitron['prot_domains'] = '.'


# =============================================================================
# Main
# =============================================================================


def main(args):
    # Check to see if bamfile can be opened and there is an index.
    if args.bam_file:
        try:
            bamfile = pysam.AlignmentFile(
                args.bam_file, 'rb', require_index=True)
        except FileNotFoundError:
            try:
                print('Building bam index file')
                pysam.index(args.bam_file)
                bamfile = pysam.AlignmentFile(
                    args.input, 'rb', require_index=True)
            except FileNotFoundError:
                print(f'There is a problem opening bam file at: {args.input}')
    else:
        bamfile = None
    try:
        db = gffutils.create_db(args.annotation_ref,
                                dbfn=args.annotation_ref + '.db',
                                disable_infer_genes=True,
                                disable_infer_transcripts=True)
        db = gffutils.FeatureDB(args.annotation_ref + '.db')
    except:
        print(f'Using annotation databse {args.annotation_ref + ".db"}')
        db = gffutils.FeatureDB(args.annotation_ref + '.db')

    if not args.out:
        args.out = args.input + '.annotated'

    if args.fasta:
        dna_seqs = []
        prot_seqs = []

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
                  'cluster_purity',
                  'exitron_prot_position',
                  'type',
                  'substitution',
                  'nmd_status_predicted',
                  'nmd_status_read_percentage',
                  'downstream_inframe_AUG',
                  'start_proximal_PTC',
                  'prot_domains',
                  'reads']
        out.write('\t'.join(header) + '\n')
        prot_df = pd.read_csv(f'{os.path.dirname(os.path.realpath(__file__))}/human_pfam.tsv', delimiter='\t') if not args.arabidopsis else pd.read_csv(
            f'{os.path.dirname(os.path.realpath(__file__))}/arabidopsis_pfam.tsv', delimiter='\t')
        exitrons = read_exitron_file(args.input)
        for exitron in exitrons:
            transcripts = exitron['transcript_id'].split(';')
            for transcript in transcripts:
                t_id = transcript.split(',')[0]
                abundance = transcript.split(',')[1]
                # dna / prot seq output not yet implemented
                res, dna, prot = categorize_exitron(
                    exitron.copy(), t_id, bamfile, db, args.genome_ref, args.arabidopsis)
                if args.fasta and dna and prot:
                    dna_seqs.append(SeqRecord(
                        dna, id=f"{exitron['gene_name']}_{t_id}", description=exitron['name']))
                    prot_seqs.append(SeqRecord(
                        prot, id=f"{exitron['gene_name']}_{t_id}", description=exitron['name']))
                # idetnify pfam domains
                if res['exitron_prot_position'] != '.':
                    get_pfam_domains(res, prot_df)
                else:
                    res['prot_domains'] = '.'

                # update abundance
                res['transcript_id'] += f',{abundance}'
                out.write('\t'.join([str(res[column])
                          for column in header]))
                out.write('\n')

    if args.fasta:
        input_fn = args.input.split('.')[:-1]
        input_fn = '.'.join(input_fn)
        SeqIO.write(dna_seqs, f'{input_fn}_dna_exitrons.fa', 'fasta')
        SeqIO.write(prot_seqs, f'{input_fn}_prot_exitrons.fa', 'fasta')


if __name__ == '__main__':
    # Get arguments
    args = parse_args()
    try:
        main(args)
    except Exception as e:
        if e.__class__.__name__ == 'InterruptedError':
            sys.stderr.write("User interrupt!")
        else:
            traceback.print_exc()
        sys.exit(1)
