#!/usr/bin/env python3
#-*- coding: utf-8 -*-
#
# ScanExitron v1 written by Tingyou Wang@Yang Lab, Hormel In
#
# ScanExitronLR written by Josh Fry@Yang Lab, University of Minnesota
#
#===============================================================================
__version__ = 'v0.2'
import sys
import os
import argparse
import pysam
import pybedtools
import subprocess
import gffutils
import re
import traceback
import multiprocessing as mp
import pandas as pd
# mp.set_start_method("spawn")
from shutil import rmtree
# from statistics import median
from Bio import pairwise2
from collections import Counter, defaultdict

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
        help="Input BAM/CRAM file (if no index is found, an index file is created).",
        required=True,
    )
    parser.add_argument(
        "-arabidopsis",
        "--arabidopsis",
        action="store_true",
        dest="arabidopsis",
        help="If set, runs exitron calling assuming genome and annotation is from arabidopsis thaliana",
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
        help="Output filename. (default: $bam_file_name.exitron)",
        default=None,
    )
    parser.add_argument(
        '-c',
        '--cores',
        action='store',
        dest='cores',
        type=int,
        help="Number of cores for parallel processing. If 0 then no parallel processing is used (default: %(default)s)",
        default=1)
    parser.add_argument(
        '-m',
        '--mapq',
        action='store',
        dest='mapq',
        type=int,
        help="Consider only reads with MAPQ >= cutoff (default: %(default)s)",
        default=50)
    parser.add_argument(
        '-j',
        '--jitter',
        action='store',
        dest='jitter',
        type=int,
        help="Exitrons with splice sites +/- jitter away will be counted as a single splicing event.  (default: %(default)s)",
        default=10)
    parser.add_argument(
        "-a",
        "--ao",
        action="store",
        dest="ao_min",
        type=int,
        help="AO cutoff (default: %(default)s)",
        default=1,
    )
    parser.add_argument(
        "-p",
        "--pso",
        action="store",
        dest="pso_min",
        type=float,
        help="PSO cutoff (default: %(default)s)",
        default=0.01,
    )
    parser.add_argument(
        "-sr",
        "--skip-realign",
        action="store_true",
        dest="skip_realign",
        help="If specified, candidate misaligned exitrons will be not be realigned. (default: %(default)s"
    )
    parser.add_argument(
        "-sa",
        "--save-abundance",
        action="store_true",
        dest="save_abundance",
        help="If specified, saves two files OUTPUT_FILENAME.isoform.normals and OUTPUT_FILENAME.isoform..exitrons to use in downstream LIQA differential expression analysis.",
    )
    parser.add_argument(
        "-md",
        "--meta_data",
        action="store",
        dest="meta_data",
        help="If specified, metadata will be written to the file (default: %(default)s",
        default=None,
    )
    parser.add_argument(
        "-s",
        "--stranded",
        action="store",
        dest="stranded",
        help="Determines how read strand is inferred. Options are 'no', 'fr-firststrand', 'fr-secondstrand'.  If 'no' for unstranded reads, the XS tag will be used. Otherwise, strand will be inferred based on library prep. See https://chipster.csc.fi/manual/library-type-summary.html for details. (default: %(default)s",
        type=str,
        choices=['no', 'fr-firststrand','fr-secondstrand'],
        default='no',
    )
    parser.add_argument(
        "-v", "--version", action="version", version="%(prog)s {}".format(__version__)
    )
    args = parser.parse_args()
    return args

    #TODO: use for future implementations of NCBI reference
    # chrms_dict = {'1':'chr1', '2':'chr2', '3':'chr3', '4':'chr4', '5':'chr5',
    #         '6':'chr6', '7':'chr7', '8':'chr8', '9':'chr9', '10':'chr10',
    #         '11':'chr11','12':'chr12', '13':'chr13', '14':'chr14', '15':'chr15',
    #         '16':'chr16','17':'chr17', '18':'chr18', '19':'chr19', '20':'chr20',
    #         '21':'chr21', '22':'chr22', 'X':'chrX', 'Y':'chrY', 'MT':'chrM'}

    # reverse_chrms_dict = dict((chrms_dict[i], i) for i in chrms_dict)

tab = str.maketrans("ACTG", "TGAC")
def rc(seq):
    return seq.upper().translate(tab)[::-1]
# No longer using config.ini

# def config_getter(config_file='config.ini'):
#     """


#     Parameters
#     ----------
#     config_file : location of config file

#     Returns
#     -------
#     dict : genome reference and annotation reference locations.

#     """
#     this_dir = os.path.dirname(os.path.realpath(__file__))
#     config_default = os.path.join(this_dir, config_file)
#     config = ConfigParser(os.environ)
#     config.read(config_default)
#     genome_ref = config.get("fasta", "genome")
#     annotation_ref = config.get("sorted GENCODE annotation", "annotation")
#     return genome_ref, annotation_ref


#===============================================================================
# Modules
#===============================================================================

def find_introns(read_iterator, stranded):
    """

    Parameters
    ----------
    read_iterator : iterator of reads from a pysam.AlignmentFile.
        Expected that the iterator will be an entire chromosome. See exitron_caller


    Returns
    -------
    introns -- counter of (intron_start, intron_end, strand)
    reads -- dictionary of reads that support the junctions in introns.  This
        will be used in later filtering steps.
    meta_data -- dictionary consisting of metadata collected along the way.
        This is used in later steps.

    """
    BAM_CREF_SKIP = 3

    introns = Counter()
    meta_data = defaultdict(list)
    reads = defaultdict(list)

    match = {0, 7, 9} # only M/=/X (0/7/8) and D (2) are related to genome position
    for r in read_iterator:
        base_position = r.pos
        read_position = 0
        # if cigarstring is * (r.cigartuples == None), unmatched, continue
        if r.cigartuples == None:
            continue
        # iterate through cigar string looking for N
        for i, (tag, nt) in enumerate(r.cigartuples):
            # if (0, X), keep track of base_position.
            # if (3, X), which corresponds to N,
            # look at match before and after
            if tag in match or (tag == 2 and nt < 30):
                base_position += nt
                read_position += nt
            elif tag == BAM_CREF_SKIP or (tag == 2 and nt >= 30):
                junc_start = base_position
                base_position += nt
                junc_end = base_position
                if stranded == 'no':
                    try:
                        if r.get_tag('ts') == '+':
                            strand = '-' if r.is_reverse else '+'
                        elif r.get_tag('ts') == '-':
                            strand = '+' if r.is_reverse else '-'
                    except KeyError:
                        # ts tag is not present, thus we don't know strand, continue
                        continue
                    if r.cigartuples[i - 1][0] == 2: junc_start -= r.cigartuples[i - 1][1]
                    if r.cigartuples[i + 1][0] == 2: junc_end += r.cigartuples[i + 1][1]
                    introns[(junc_start, junc_end, strand)] += 1
                    reads[(junc_start, junc_end, strand)].append(r.query_name)
                else:
                    if stranded == 'fr-firststrand':
                        strand = '+' if (r.is_read2 and not r.is_reverse) or \
                                        (r.is_read1 and r.is_reverse) else '-'
                    elif stranded == 'fr-secondstrand':
                        strand = '+' if (r.is_read1 and not r.is_reverse) or \
                                        (r.is_read2 and r.is_reverse) else '-'
                    introns[(junc_start, junc_end, strand)] += 1
                    reads[(junc_start, junc_end, strand)].append(r.query_name)

    return introns, reads, meta_data


def exitron_caller(bamfile, referencename, chrm, db, arabidopsis = False, stranded = 'no', mapq = 50, jitter = 10):
    """


    Parameters
    ----------
    bamfile : pysam.AlignmentFile
    referencename : str
    chrms : list
    mapq : int, optional
        Only considers reads from bamfile with quality >= mapq. The default is 50.

    Returns
    -------
    List of unfiltered exitrons.  Each exitron is a dictionary of features.

    """
    # add jitter.  This will later be a parameter
    # jitter = 10

    # intron_bed = []

    introns, reads, meta_data = find_introns(
        (read for read in bamfile.fetch(chrm) if read.mapping_quality >= mapq),
        stranded
        )

    # for start, stop, strand in introns:
    #     #-1 and +2 so that we can capture the ends and beginning of adjacent transcripts
    #     #this allows us to determine whether there is a known donor or acceptor site
    #     #in very rare cases, inron may be at pos 0, in which case jitter may cause overflow.
    #     intron_bed.append((chrm, max(start - 1 - jitter, 0), max(stop + 1 + jitter + 1,0), introns[(start, stop, strand)], 0, strand))

    # if not bool(introns):
    #     # No introns were found.
    #     return ([], reads, meta_data)

    # intron_bed_file = pybedtools.BedTool(intron_bed)
    # gtf_anno_sorted = pybedtools.BedTool(referencename)
    # # To deal with memory issues and multi-processing, we require that annotation is sorted.
    # # One line bash command: awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k4,4n -k5,5n"}' in.gtf > out_sorted.gtf
    # intersection = gtf_anno_sorted.intersect(intron_bed_file.sort(), s = True, wo = True, sorted = True)

    # del gtf_anno_sorted
    # del intron_bed_file
    gtf = pysam.TabixFile(referencename, parser = pysam.asGTF())
    exitrons = []
    exitrons_added = []
    known_splices = set()

    # load blacklist of known false-positive hot spots
    this_dir = os.path.dirname(os.path.realpath(__file__))
    blacklist_path = os.path.join(this_dir, 'blacklist.tsv')
    try:
        with open(blacklist_path, 'r') as b:
            b.readline()
            blacklist = [l.split('\t')[1].rstrip() for l in b]
    except:
            blacklist = []

    for intron in introns:
        intron_start = intron[0] - 1 - jitter
        intron_end = intron[1] + 1 + jitter + 1
        intron_witnesses = introns[intron]
        intersection = gtf.fetch(chrm, intron_start, intron_end)
        for feature in intersection:
            # Check for intron within coding exon.
            if feature.strand != intron[2]: continue # gtf.fetch is strand agnostic
            region_type = feature.feature
            region_start = feature.start
            region_end = feature.end
            if arabidopsis and (region_type not in ['CDS', 'exon', 'three_prime_UTR', 'five_prime_UTR']): continue
            gene_name = feature.gene_name if not arabidopsis else feature.attrs['Parent'].split('.')[0]
            gene_id = feature.gene_id if not arabidopsis else feature.attrs['Parent'].split('.')[0]

            # Use the ends to check for known donors or acceptors
            if region_type == 'exon':
                transcript_id = feature.transcript_id if not arabidopsis else feature.attrs['Parent'].split(',')[0]
                if intron_start in range(region_end - jitter*2, region_end + 1) and \
                    (chrm, intron_start + 1 + jitter, intron_end - 2  - jitter + 1, 'D') not in known_splices and \
                    feature.end != db[transcript_id].end:
                    # intron matches a known donor
                    known_splices.add((chrm, intron_start + 1 + jitter, intron_end - 2  - jitter + 1, 'D'))
                if intron_end in range(region_start, region_start + 1 + jitter*2) and \
                    (chrm, intron_start + 1 + jitter, intron_end - 2  - jitter + 1, 'A') not in known_splices and \
                    feature.start + 1 != db[transcript_id].start:
                    # intron matches a known acceptor
                    known_splices.add((chrm, intron_start + 1 + jitter, intron_end - 2  - jitter + 1, 'A'))

            elif region_type == 'CDS' and region_start < intron_start + 1 + jitter \
                and region_end > intron_end - 1 - jitter:
                    transcript_id = feature.transcript_id if not arabidopsis else feature.attrs['Parent'].split(',')[0]
                    if gene_id in blacklist: continue
                    if (intron_start, intron_end, region_type) not in exitrons_added:
                        exitrons.append({'chrom':chrm,
                                        'start':intron_start + 1 + jitter,
                                        'end':intron_end - 2 - jitter + 1, #plus 1 because of bedtools conventions,
                                        'name':f'{gene_name}d{intron_start + 1 + jitter}-{intron_end - 2 - jitter + 1}',
                                        'region':region_type,
                                        'ao':intron_witnesses,
                                        'strand':feature.strand,
                                        'gene_name':gene_name,
                                        'gene_id':gene_id,
                                        'length':intron_end - 2 - jitter + 1 - (intron_start + 1 + jitter) - 1,
                                        'splice_site':'splice_site',
                                        'transcript_id':transcript_id})
                        exitrons_added.append((intron_start, intron_end, region_type))
                    elif 'tag' in feature.keys():
                        if ('basic' in feature.tag or
                              'CCDS' in feature.tag or
                              'ccdsid' in feature.keys()):
                            #TODO pybedtools only grants access to ONE tag, we need to look at all the tags
                            #submit an issue with pybedtools
                            for e in exitrons:
                                if e['name'] == f'{gene_name}d{intron_start + 1 + jitter}-{intron_end - 2 - jitter + 1}':
                                    e['transcript_id'] = transcript_id
                                    break

        # del intersection

    return ([exitron for exitron in exitrons if (((exitron['chrom'], exitron['start'], exitron['end'], 'D') not in known_splices and
                                                 ((exitron['chrom'], exitron['start'], exitron['end'], 'A') not in known_splices)))],
            reads,
            meta_data)



def filter_exitrons(exitrons, reads, bamfile, genome, meta_data, db, skip_realign, mapq = 50, pso_min = 0.005, ao_min = 1, jitter = 10):
    """
    Parameters
    ----------
    exitrons : list
        list of unfiltered exitrons from exitron_caller.
    reads : dict
        Each intron is a key, and the value is list of (read_seq, ref_seq, left_anchor, right_anchor)
    bamfile : pysam.AlignmentFile
    genome : pysam.libcfaidx.FastaFile
        Random access to reference genome.
    meta_data : dict
    verbose : bool
        If true, we report optional statistics
    mapq : TYPE, optional
        Only considers reads from bamfile with quality >= mapq. The default is 50.
        This is needed to calculate pso.
    pso_min : float, optional
        Number reads that witness the exitron over the number of reads within the
        spliced out region. The default is 0.05.
    ao_min : int, optional
        Minimum number of reads witnessing the exitron. The default is 3.

    pso_ufmin : float, optional
        Number reads that witness the exitron over the number of reads within the
        spliced out region, before filtering (used for backwards compatibility). The default is 0.05.
    ao_ufmin : int, optional
        Minimum number of reads witnessing the exitron, before filtering (used for
        backwards compatibility). The default is 3.
    anchor_min : int, optional
        Minimum anchor length.  The default is 5.

    Returns
    -------
    filtered exitrons and meta data
    """

    """
    The Plan:

        iterate through exitrons, finding clusters
        i.e.: collect intervals that are +/- 10 away


    """
    if not exitrons:
        return [], meta_data

    jitter = jitter*2
    #Need to compute reverse complement for splice junctions.
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    res = []
    groups = {'+':[],
              '-':[]}
    collection = {'+':[],
                  '-':[]}

    # need to sort exitrons for the clustering algorithm
    exitrons.sort(key = lambda x: (x['start'], x['end']))

    # DEBUG: [(e['start'], e['end']) for e in exitrons]
    # filter one exitron at a time
    for exitron in exitrons:
        ao = exitron['ao']
        strand = exitron['strand']

        if collection[strand]:
            if (collection[strand][-1]['start'] - jitter <= exitron['start'] <= collection[strand][-1]['start'] + jitter and
                collection[strand][-1]['end'] - jitter <= exitron['end'] <= collection[strand][-1]['end'] + jitter):
                collection[strand].append(exitron)
            else:
                groups[strand].append(collection[strand])
                collection[strand] = [exitron]
        else:
            collection[strand] = [exitron]

    groups['+'].append(collection['+'])
    groups['-'].append(collection['-'])

    for strand in groups:
        for group in groups[strand]:
            if not group:
                continue # no exitrons found in this strand
            # calculate canonical spice sites and append read names.
            for e in group:
                start = e['start']
                end = e['end']
                # pos = db[e['transcript_id']].end - start if strand == '-' else start - db[e['transcript_id']].start
                # genome_seq = genome[e['chrom']][start:end - 1].upper()
                if strand == '+':
                    pos = start - db[e['transcript_id']].start + 1
                    genome_seq = db[e['transcript_id']].sequence(genome)[pos:pos + e['length']]
                    e['splice_site'] = genome_seq[:2] + '-' + genome_seq[-2:]
                elif strand == '-':
                    pos = db[e['transcript_id']].end - start
                    genome_seq = db[e['transcript_id']].sequence(genome)[pos - e['length']:pos]
                    right = genome_seq[:2]
                    left = genome_seq[-2:]
                    e['splice_site'] = genome_seq[:2] + '-' + genome_seq[-2:]
            try:
                consensus_e = max([e for e in group if e['splice_site'] in ['GT-AG','GC-AG','AT-AC']],
                                  key = lambda e: e['ao'])
                tot_ao = sum(e['ao'] for e in group)
                consensus_e['consensus_prop'] = round(consensus_e['ao']/tot_ao, ndigits = 2)
                consensus_e['ao'] = tot_ao
                consensus_reads = ''
                for e in group:
                    consensus_reads += ',' + ','.join(reads[(e['start'],
                                                           e['end'] - 1,
                                                           e['strand'])])
                consensus_e['reads'] = consensus_reads
            except ValueError:
                continue # this occurs when there are no cannonical splice sites within the group

            ao = consensus_e['ao']
            start = consensus_e['start']
            end = consensus_e['end']
            chrm = consensus_e['chrom']

            if ao < ao_min:
                meta_data['low_ao'].append(consensus_e)
                continue

            # We subtract 1 because these coords are in BED format.
            mid = (start+end)/2
            a = bamfile.count(chrm, start = start - 1, stop = start, read_callback = lambda x: x.mapq > mapq)
            b = bamfile.count(chrm, start = end - 1, stop = end, read_callback = lambda x: x.mapq > mapq)
            c = bamfile.count(chrm, start = mid - 1, stop = mid, read_callback = lambda x: x.mapq > mapq)

            pso = ao/((a + b + c - ao*3)/3.0 + ao)
            dp = int(ao/pso) if pso > 0 else 0

            # Check whether attributes exceed minimum values
            if pso >= pso_min:
                consensus_e['pso'] = round(pso, ndigits = 4)
                consensus_e['dp'] = dp
                consensus_e['a'] = a
                consensus_e['b'] = b
                consensus_e['c'] = c
                res.append(consensus_e)
            else:
                meta_data['low_pso'].append(consensus_e)


    # realignmnet may double count
    if not skip_realign and res:
        # realign
        print(f'Realigning exitrons in {res[0]["chrom"]}')
        called_reads = ''.join(e['reads'] for e in res)
        for exitron in res:
            e_start = exitron['start']
            e_end = exitron['end']
            chrm = exitron['chrom']
            e_length = exitron['length']

            for exon in db.children(db[exitron['transcript_id']], featuretype = 'exon', order_by = 'start'):
                if exon.start <= e_start <= exon.end:
                    exon_start = exon.start
                    exon_end = exon.end
                    break


            for read in bamfile.fetch(chrm, e_start, e_end):
                if read.query_name in called_reads or read.mapping_quality < mapq:
                    continue
                if (any([max(0, min(exon_end, i[1]) - max(exon_start, i[0])) > 0
                        for i in bamfile.find_introns([read]).keys()])):
                    pos = [p for p in read.get_aligned_pairs() if (p[1] != None and exon_start <= p[1] <= exon_end)]
                    try:
                        start = min(p[0] for p in pos if p[0] != None)
                        end = max(p[0] for p in pos if p[0] != None)
                    except:
                        continue # no nt overlap with exon

                    try:
                        r_seq = read.seq[start:end].upper()
                    except:
                        continue # strangely, sometimes pysam returns None from read.seq
                    if not r_seq: continue
                    # exon.sequence(-) is much faster than using pysam Fast
                    g_seq = exon.sequence(genome).upper() if exon.strand == '+' else rc(exon.sequence(genome).upper())

                    exitron_seq = g_seq[e_start - exon.start + 1 - 2: e_end - exon.start + 2]

                    left = exitron_seq[:2]
                    right = exitron_seq[-2:]
                    exitron_seq = exitron_seq[2:-2]
                    similarity = pairwise2.align.localms(read.seq[start-e_length:end+e_length], exitron_seq, 2, -1, -2, -1, score_only = True)
                    similarity = similarity/(e_length*2) if similarity else 1 # if no alignment, just continue
                    if similarity <= 0.7:
                        alignment = pairwise2.align.localms(g_seq, r_seq, 4, -2, -2, 0)[:10]
                        if (any(re.findall(f'{left}--*', aln.seqB) and (e_length - 10 - len(r_seq)*0.05 <= aln.seqB.count('-') <= e_length + 10 + len(r_seq)*0.05) for aln in alignment[:10]) or
                             any(re.findall(f'--*{right}', aln.seqB) and (e_length - 10 - len(r_seq)*0.05 <= aln.seqB.count('-') <= e_length + 10 + len(r_seq)*0.05) for aln in alignment[:10])):
                                exitron['reads'] += f',REALIGNED_{read.query_name}'
                                exitron['ao'] += 1
                                exitron['pso'] = exitron['ao']/((exitron['a'] + exitron['b'] + exitron['c'] - exitron['ao']*3)/3.0 + exitron['ao'])
                                called_reads += f',{read.query_name}'

    # for group in groups['-']:
    #     if not group:
    #         break
    #     # same but use reverse complement for splice-site
    #     for e in group:
    #         start = e['start']
    #         end = e['end']
    #         genome_seq = genome[e['chrom']][start:end - 1].upper()
    #         right = "".join(complement.get(base, base) for base in reversed(genome_seq[:2]))
    #         left = "".join(complement.get(base, base) for base in reversed(genome_seq[-2:]))
    #         e['splice_site'] = left + '-' + right
    #     try:
    #         consensus_e = max([e for e in group if e['splice_site'] in ['GT-AG','GC-AG','AT-AC']],
    #                           key = lambda e: e['ao'])

    #         consensus_e['ao'] = sum(e['ao'] for e in group)
    #         processed_exitrons.append(consensus_e)
    #     except:
    #         pass # this occurs when there are no cannonical splice sites within group
    # for exitron in processed_exitrons:
    #     ao = exitron['ao']
    #     start = exitron['start']
    #     end = exitron['end']
    #     chrm = exitron['chrom']

    #     if ao < ao_min:
    #         meta_data['low_ao'].append(exitron)
    #         continue

    #     # We subtract 1 because these coords are in BED format.
    #     mid = (start+end)/2
    #     a = bamfile.count(chrm, start = start - 1, stop = start, read_callback = lambda x: x.mapq > mapq)
    #     b = bamfile.count(chrm, start = end - 1, stop = end, read_callback = lambda x: x.mapq > mapq)
    #     c = bamfile.count(chrm, start = mid - 1, stop = mid, read_callback = lambda x: x.mapq > mapq)

    #     pso = ao/((a + b + c - ao*3)/3.0 + ao)
    #     psi = 1 - pso
    #     dp = int(ao/pso) if pso > 0 else 0

    #     # Check whether attributes exceed minimum values
    #     if pso >= pso_min:
    #         exitron['pso'] = pso
    #         exitron['psi'] = psi
    #         exitron['dp'] = dp
    #         res.append(exitron)
    #     else:
    #         meta_data['low_pso'].append(exitron)
    return res, meta_data


def identify_transcripts(exitrons, db, bamfilename, tmp_path, save_abundance, out_fn, arabidopsis, cores):
    bamfile = pysam.AlignmentFile(bamfilename, 'rb', require_index = True)

    # construct new bamfile
    tmp_bamfile_exitrons = pysam.AlignmentFile(tmp_path + '/e_tmp.bam', 'wb', template = bamfile)
    if save_abundance:
        tmp_bamfile_normals = pysam.AlignmentFile(tmp_path + '/n_tmp.bam', 'wb', template = bamfile)

    for e in exitrons:
        e_reads = e['reads'].split(',')
        # fetch reads at exitron junction
        for read in bamfile.fetch(e['chrom'], start = int(e['start']), stop = int(e['end'])):
            if read.query_name in e_reads:
                tmp_bamfile_exitrons.write(read)
            elif save_abundance:
                tmp_bamfile_normals.write(read)
    tmp_bamfile_exitrons.close()
    if save_abundance: tmp_bamfile_normals.close()
    bamfile.close()

    pysam.sort('-@', str(cores), '-o', tmp_path + '/e_tmp_sorted.bam', tmp_path + '/e_tmp.bam')
    if save_abundance: pysam.sort('-o', tmp_path + '/n_tmp_sorted.bam', tmp_path + '/n_tmp.bam')
    pysam.index(tmp_path + '/e_tmp_sorted.bam')
    if save_abundance: pysam.index(tmp_path + '/n_tmp_sorted.bam')

    # build a small gtf file of only thoes exitron spliced genes
    with open(tmp_path + '/tmp.gtf', 'w') as f:
        for e in exitrons:
            for region in db.children(db[e['gene_id']], order_by = 'start'):
                if arabidopsis and region.featuretype == 'exon':
                    gene_id = e['gene_id']
                    transcript_id = region.attributes['Parent'][0]
                    region = str(region).split('\t')
                    region[8] = f'gene_id "{gene_id}"; transcript_id "{transcript_id}";'
                    f.write('\t'.join(region) + '\n')
                elif not arabidopsis:
                    f.write(str(region) + '\n')


    # run liqa to create refgene
    subprocess.run(['liqa',
                    '-task',
                    'refgene',
                    '-ref',
                    f'{tmp_path + "/tmp.gtf"}',
                    '-format',
                    'gtf',
                    '-out',
                    f'{tmp_path + "/tmp.refgene"}'])

    # run liqa to quantify isoform expression
    jitter = 10 #TODO make this an argument
    subprocess.run(['liqa',
                    '-task',
                    'quantify',
                    '-refgene',
                    f'{tmp_path + "/tmp.refgene"}',
                    '-bam',
                    f'{tmp_path + "/e_tmp_sorted.bam"}',
                    '-out',
                    f'{tmp_path + "/isoform_estimates.out"}',
                    '-max_distance',
                    f'{jitter}',
                    '-f_weight',
                    '0'])
    if save_abundance:
        # run liqa to quantify isoform expression

        jitter = 10 #TODO make this an argument
        subprocess.run(['liqa',
                        '-task',
                        'quantify',
                        '-refgene',
                        f'{tmp_path + "/tmp.refgene"}',
                        '-bam',
                        f'{tmp_path + "/n_tmp_sorted.bam"}',
                        '-out',
                        f'{os.path.splitext(out_fn)[0] + ".isoform.normals"}',
                        '-max_distance',
                        f'{jitter}',
                        '-f_weight',
                        '0'])

    ie = pd.read_csv(f'{tmp_path}/isoform_estimates.out', sep='\t')
    if save_abundance: ie.to_csv(os.path.splitext(out_fn)[0] + ".isoform.exitrons", sep = '\t', index = False)
    for e in exitrons:
        gene = e['gene_name']
        ie_slice = ie[ie['GeneName'] == gene].sort_values(ascending = False, by = 'RelativeAbundance')
        transcripts = list(ie_slice[ie_slice['RelativeAbundance'] > 0.1]['IsoformName'])
        try:
            if not transcripts: transcripts = list(ie_slice[ie_slice['RelativeAbundance'] == max(ie_slice['RelativeAbundance'])]['IsoformName'])
        except ValueError:
            # transcript could not be measured by liqa
            e['transcript_id'] += ',NA' # revert back to default transcript_id
            continue
        t_str = ''
        for t in transcripts:
            t_str += f'{t},{round(float(ie_slice[ie_slice["IsoformName"] == t]["RelativeAbundance"]), 4)};'
        e['transcript_id'] = t_str[:-1] # leave off trailing ;
    return exitrons


# for reference, here is the short read version:

# def filter_exitrons(exitrons, reads, bamfile, genome, meta_data, verbose, mapq = 50, pso_min = 0.005, ao_min = 1, pso_ufmin = 0, ao_ufmin = 0, anchor_min = 5):
#     """
#     Parameters
#     ----------
#     exitrons : list
#         list of unfiltered exitrons from exitron_caller.
#     reads : dict
#         Each intron is a key, and the value is list of (read_seq, ref_seq, left_anchor, right_anchor)
#     bamfile : pysam.AlignmentFile
#     genome : pysam.libcfaidx.FastaFile
#         Random access to reference genome.
#     meta_data : dict
#     verbose : bool
#         If true, we report optional statistics
#     mapq : TYPE, optional
#         Only considers reads from bamfile with quality >= mapq. The default is 50.
#         This is needed to calculate pso.
#     pso_min : float, optional
#         Number reads that witness the exitron over the number of reads within the
#         spliced out region. The default is 0.05.
#     ao_min : int, optional
#         Minimum number of reads witnessing the exitron. The default is 3.

#     pso_ufmin : float, optional
#         Number reads that witness the exitron over the number of reads within the
#         spliced out region, before filtering (used for backwards compatibility). The default is 0.05.
#     ao_ufmin : int, optional
#         Minimum number of reads witnessing the exitron, before filtering (used for
#         backwards compatibility). The default is 3.
#     anchor_min : int, optional
#         Minimum anchor length.  The default is 5.

#     Returns
#     -------
#     filtered exitrons and meta data
#     """
#     #Need to compute reverse complement for splice junctions.
#     complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
#     res = []

#     # filter one exitron at a time
#     for exitron in exitrons:
#         chrm = exitron['chrom']
#         ao = exitron['ao']
#         start = exitron['start']
#         end = exitron['end']
#         strand = exitron['strand']

#         if ao < ao_ufmin:
#             continue #not enough unfiltered reads support this exitron
#         read_data = reads[(start, end - 1, strand)]

#         genome_seq = genome[chrm][start-10:end - 1+10].upper()
#         intron_seq = genome_seq[10:-10]


#         if strand == '+':
#             splice_site = intron_seq[:2] + '-' + intron_seq[-2:]
#             motif = genome_seq[:10] + '|' + splice_site[:2] + '-' + splice_site[-2:] + '|' + genome_seq[-10:]
#         elif strand == '-':
#             left = "".join(complement.get(base, base) for base in reversed(intron_seq[-2:]))
#             right = "".join(complement.get(base, base) for base in reversed(intron_seq[:2]))
#             splice_site = left + '-' + right
#             rgs = "".join(complement.get(base, base) for base in reversed(genome_seq[:10]))
#             lgs = "".join(complement.get(base, base) for base in reversed(genome_seq[-10:]))
#             motif = lgs + '|' + splice_site[:2] + '-' + splice_site[-2:] + '|' + rgs

#         exitron['splice_site'] = splice_site
#         if verbose:
#             exitron['splice_motif'] = motif

#         if splice_site.upper() not in ['GT-AG','GC-AG','AT-AC']:
#             meta_data['non_canonical'].append(exitron)
#             continue

#         right_anchor = (0,'')
#         left_anchor = (0,'')
#         ao_true = 0
#         for read in read_data:

#             l_anchor_len = read[1]
#             r_anchor_len = read[2]

#             # check 5' intron and 3' exon similarity
#             three_prime_exon_r = read[0][read[3]:read[3] + r_anchor_len]
#             five_prime_intron_r = intron_seq[:len(three_prime_exon_r)]
#             # check 5' exon and 3' intron similarity
#             five_prime_exon_r =  read[0][read[3] - l_anchor_len:read[3]]
#             three_prime_intron_r = intron_seq[-len(five_prime_exon_r):]

#             if three_prime_exon_r != five_prime_intron_r and five_prime_exon_r != three_prime_intron_r:
#                 ao_true += 1


#         if ao_true == 0 and ao_min > 0:
#             exitron['ao_unfiltered'] = ao
#             meta_data['no_anchors_after_filtering'].append(exitron)
#             continue

#         # We subtract 1 because these coords are in BED format.
#         mid = (start+end)/2
#         a = bamfile.count(chrm, start = start - 1, stop = start, read_callback = lambda x: x.mapq > mapq)
#         b = bamfile.count(chrm, start = end - 1, stop = end, read_callback = lambda x: x.mapq > mapq)
#         c = bamfile.count(chrm, start = mid - 1, stop = mid, read_callback = lambda x: x.mapq > mapq)


#         pso = ao/((a + b + c - ao*3)/3.0 + ao)
#         psi = 1 - pso
#         dp = int(ao/pso) if pso > 0 else 0

#         pso_true = ao_true/((a + b + c - ao_true*3)/3.0 + ao_true)
#         psi_true = 1 - pso_true
#         dp_true = int(ao_true/pso_true) if pso_true > 0 else 0

#         # Check whether attributes exceed minimum values
#         if (pso >= pso_ufmin and
#             pso_true >= pso_min and
#             ao_true >= ao_min):

#             exitron['ao'] = ao_true
#             exitron['pso'] = pso_true
#             exitron['psi'] = psi_true
#             exitron['dp'] = dp_true

#             exitron['ao_unfiltered'] = ao
#             exitron['pso_unfiltered'] = pso
#             exitron['psi_unfiltered'] = psi
#             exitron['dp_unfiltered'] = dp

#             exitron['delta_ao'] = ao - ao_true


#             left_anchor = max(read_data, key = lambda x: x[1])
#             right_anchor = max(read_data, key = lambda x: x[2])
#             if left_anchor[2] < anchor_min or right_anchor[3] < anchor_min:
#                 meta_data['low_anchor_exitron'].append(exitron)
#             else:
#                 res.append(exitron)
#     return res, meta_data


#===============================================================================
# Main
#===============================================================================


def exitrons_in_chrm(bamfilename, referencename, genomename, chrm, mapq, pso_min, ao_min, jitter, stranded, skip_realign, arabidopsis):
    """
    Wrapper that calls main functions *per chromosome*.
    """
    print(f'Finding exitrons in {chrm}')
    sys.stdout.flush()
    bamfile = pysam.AlignmentFile(bamfilename, 'rb', require_index = True)
    db = gffutils.FeatureDB(referencename + '.db')
    exitrons, reads, meta_data  = exitron_caller(bamfile,
                              referencename,
                              chrm,
                              db,
                              arabidopsis,
                              stranded,
                              mapq,
                              jitter)
    exitrons, meta_data = filter_exitrons(exitrons,
                    reads,
                    bamfile,
                    genomename,
                    meta_data,
                    db,
                    skip_realign,
                    mapq,
                    pso_min,
                    ao_min)
    bamfile.close()
    del db

    return exitrons, meta_data, chrm

def main(tmp_path):
    # Get arguments
    args = parse_args()
    # Define chrms
    chrms = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5',
        'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
        'chr11','chr12', 'chr13', 'chr14', 'chr15',
        'chr16','chr17', 'chr18', 'chr19', 'chr20',
        'chr21', 'chr22', 'chrX', 'chrY'] if not args.arabidopsis else \
        ['Chr1', 'Chr2', 'Chr3', 'Chr4',
         'Chr5', 'ChrM', 'Chrchloroplast']

    # genome_ref, annotation_ref = config_getter('config.ini') # no longer using config.ini

    # Check if bamfile can be opened and there is an index.
    try:
        bamfile = pysam.AlignmentFile(args.input, 'rb', require_index = True)
    except FileNotFoundError:
        try:
            print('Building bam index file')
            pysam.index(args.input)
            bamfile = pysam.AlignmentFile(args.input, 'rb', require_index = True)
        except FileNotFoundError:
            print(f'ERROR: There is a problem opening bam file at: {args.input}')
            sys.exit(1)
    bamfile.close()

    # Check if LIQA is available
    try:
        subprocess.run(['liqa'], capture_output = True)
    except FileNotFoundError:
        print('ERROR: Unable to locate LIQA. Please install with "pip install liqa"')
        sys.exit(1)

    # Check for gziped annotation
    if args.annotation_ref[-2:] != 'gz':
        print('ERROR: Annotation file is required to be zipped in .gz format. Please compress your GTF/GFF file with a command such as: gzip -c in.gtf > out.gtf.gz')
        sys.exit(1)

    # Prepage gffutils database
    try:
        db = gffutils.FeatureDB(args.annotation_ref + '.db')
        print(f'Using annotation databse {args.annotation_ref + ".db"}')
    except ValueError:
        print('Preparing annotation database... This may take awhile ... ')
        db = gffutils.create_db(args.annotation_ref,
                                dbfn = args.annotation_ref + '.db',
                                disable_infer_genes = True,
                                disable_infer_transcripts = True)
        db = gffutils.FeatureDB(args.annotation_ref + '.db')
        print(f'Using annotation databse {args.annotation_ref + ".db"}')

    # check for blacklist
    this_dir = os.path.dirname(os.path.realpath(__file__))
    blacklist_path = os.path.join(this_dir, 'blacklist.tsv')

    try:
        b = open(blacklist_path, 'r')
        b.close()
    except FileNotFoundError:
        print('No blacklist found -- continuing without it, beware of false positives')


    # Begin exitron calling
    global results
    global meta_data_out
    results = {}
    meta_data_out = defaultdict(list)

    def collect_result(output):
        exitrons = output[0]
        meta_data_result = output[1]
        chrm = output[2]
        print(f'Collecting data from {chrm}')
        sys.stdout.flush()
        results[chrm] = exitrons
        for key in meta_data_result:
            meta_data_out[key] = [item for item in meta_data_result[key]]

    if int(args.cores) > 1:
        pool = mp.Pool(int(args.cores))
        threads = []
        for chrm in chrms:
            threads.append(pool.apply_async(exitrons_in_chrm, args = (args.input,
                                                        args.annotation_ref,
                                                        args.genome_ref,
                                                        chrm,
                                                        args.mapq,
                                                        args.pso_min,
                                                        args.ao_min,
                                                        args.jitter,
                                                        args.stranded,
                                                        args.skip_realign,
                                                        args.arabidopsis), callback = collect_result))
        pool.close()
        for t in threads:
            t.get() # this gets any exceptions raised
        pool.join()
    else:
        for chrm in chrms:
            sys.stdout.flush()
            output = exitrons_in_chrm(args.input,
                                        args.annotation_ref,
                                        args.genome_ref,
                                        chrm,
                                        args.mapq,
                                        args.pso_min,
                                        args.ao_min,
                                        args.jitter,
                                        args.stranded,
                                        args.skip_realign,
                                        args.arabidopsis)
            collect_result(output)
    exitrons = []
    for chrm in results:
        exitrons.extend(results[chrm])

    out_file_name = args.out
    if not out_file_name:
        out_file_name = f'{os.path.splitext(args.input)[0]}.exitron'

    print('Quantifying transcripts.')
    sys.stdout.flush()
    # update transcripts
    identify_transcripts(exitrons,
                         db,
                         args.input,
                         tmp_path,
                         args.save_abundance,
                         out_file_name,
                         args.arabidopsis,
                         args.cores)
    print(f'Finished exitron calling and filtering. Printing to {out_file_name}')
    sys.stdout.flush()
    with open(out_file_name, 'w') as out:
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
                  'reads']
        #write header
        for column in header:
            out.write(column + '\t')
        out.write('\n')
        for exitron in exitrons:
            for column in header:
                out.write(str(exitron[column]) + '\t')
            out.write('\n')

    if args.meta_data:
        with open(args.meta_data, 'w') as out:
            out.write('Non-canonical exitrons:\n')
            for exitron in meta_data_out['non_canonical']:
                for key in exitron:
                    out.write(str(exitron[key])+ '\t')
                out.write('\n')
            out.write('\nExitrons that had no anchors after filtering:\n')
            for exitron in meta_data_out['no_anchors_after_filtering']:
                for key in exitron:
                    out.write(str(exitron[key])+ '\t')
                out.write('\n')
            out.write('\nExitrons that were filtered due to low anchor length:\n')
            for exitron in meta_data_out['low_anchor_exitron']:
                for key in exitron:
                    out.write(str(exitron[key])+ '\t')
                out.write('\n')

    # Clear tmp directory
    pybedtools.helpers.cleanup(remove_all=True)
    rmtree(tmp_path)


if __name__ == '__main__':
    # Set tmp directory
    this_dir = os.path.dirname(os.path.realpath(__file__))
    tmp_path = os.path.join(this_dir, f'scanexitron_tmp{os.getpid()}')
    try: os.mkdir(tmp_path)
    except FileExistsError: pass
    pybedtools.helpers.set_tempdir(tmp_path)
    main(tmp_path)
    # try:
    #     main(tmp_path)
    # except Exception as e:
    #     if e.__class__.__name__ == 'InterruptedError':
    #         sys.stderr.write("User interrupt!")
    #     else:
    #         traceback.print_tb(e.__traceback__)
    #     pybedtools.helpers.cleanup(remove_all=True)
    #     rmtree(tmp_path)
    #     sys.exit(1)