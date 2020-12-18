#!/usr/bin/env python

import os
import pysam
import argparse
import subprocess
import logging

from uuid import uuid4
from collections import defaultdict as dd

from bx.intervals.intersection import Intersecter, Interval

FORMAT = '%(asctime)s %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def load_falib(infa):
    seqdict = {}

    with open(infa, 'r') as fa:
        seqid = ''
        seq   = ''
        for line in fa:
            if line.startswith('>'):
                if seq != '':
                    seqdict[seqid] = seq
                seqid = line.lstrip('>').strip().split()[0]

                if ':' not in seqid:
                    logger.error("TE reference sequence .fasta headers must be in the format superfamily:subfamily e.g. >L1:L1Ta")
                    return None

                seq = ''
            else:
                assert seqid != ''
                seq = seq + line.strip()

    if seqid not in seqdict and seq != '':
        seqdict[seqid] = seq

    return seqdict


def align(qryseq, refseq, elt='PAIR', minmatch=85.0):
    rnd = str(uuid4())
    tgtfa = 'tmp.' + rnd + '.tgt.fa'
    qryfa = 'tmp.' + rnd + '.qry.fa'

    tgt = open(tgtfa, 'w')
    qry = open(qryfa, 'w')

    tgt.write('>ref' + '\n' + refseq + '\n')
    qry.write('>qry' + '\n' + qryseq + '\n')

    tgt.close()
    qry.close()

    cmd = ['exonerate', '--bestn', '1', '-m', 'ungapped', '--showalignment','0', '--ryo', elt + '\t%s\t%qab\t%qae\t%tab\t%tae\t%pi\t%qS\t%tS\n', qryfa, tgtfa]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    best = []
    topscore = 0

    for pline in p.stdout.readlines():
        pline = pline.decode()
        if pline.startswith(elt):
            c = pline.strip().split()
            if int(c[1]) > topscore and float(c[6]) >= minmatch:
                topscore = int(c[1])
                best = c

    os.remove(tgtfa)
    os.remove(qryfa)

    return best


def mm2_search(ref_fa, trd_fa, ref_elts, nr_elts, tldr_forest, window=100):
    FNULL = open(os.devnull, 'w')

    ref_tbx = pysam.Tabixfile(ref_elts)
    nr_tbx = pysam.Tabixfile(nr_elts)

    result = {}

    mm2_cmd  = ['minimap2', '-x', 'map-ont', '-Y', '-a', ref_fa, trd_fa]
    view_cmd = ['samtools', 'view', '-b', '-']
    aln  = subprocess.Popen(mm2_cmd, stdout=subprocess.PIPE, stderr=FNULL)
    view = subprocess.Popen(view_cmd, stdin=aln.stdout, stdout=subprocess.PIPE, stderr=FNULL)

    bamstream = pysam.AlignmentFile(view.stdout, 'rb')

    for read in bamstream:
        if None in (read.reference_name, read.reference_start, read.reference_end, read.mapq):
            continue

        align_info = ['%s:%d:%d:%d' % (read.reference_name, read.reference_start, read.reference_end, read.mapq)]
        if read.reference_name in ref_tbx.contigs:
            for ref_elt in ref_tbx.fetch(read.reference_name, read.reference_start-window, read.reference_end+window):
                align_info.append(':'.join(ref_elt.strip().split()) + ':REF')

        if read.reference_name in nr_tbx.contigs:
            for nr_elt in nr_tbx.fetch(read.reference_name, read.reference_start-window, read.reference_end+window):
                align_info.append(':'.join(nr_elt.strip().split()) + ':NR')

        if read.reference_name in tldr_forest:
            for tldr_elt in tldr_forest[read.reference_name].find(read.reference_start-window, read.reference_end+window):
                rec = tldr_elt.value
                align_info.append(':'.join((rec['Chrom'], rec['Start'], rec['End'], rec['Subfamily'], rec['UUID'], 'TLDR')))

        result[read.query_name] = align_info

    return result


def filter(rec, telocs, maplocs, mask, window=10000):
    filters = []

    rec_chrom = rec['Chrom']
    rec_start = int(rec['Start'])-window
    rec_end   = int(rec['End'])+window

    for maploc in maplocs:
        map_chrom, map_start, map_end, map_qual = maploc.split(':')
        map_start = int(map_start)-window
        map_end   = int(map_end)+window
        map_qual  = int(map_qual)

        if map_chrom == rec_chrom:
            if min(rec_end, map_end) - max(rec_start, map_start) > 0:
                filters.append('InsClose')

        if mask:
            map_chrom, orig_map_start, orig_map_end, map_qual = maploc.split(':')
            orig_map_start = int(orig_map_start)
            orig_map_end   = int(orig_map_end)
            map_qual       = int(map_qual)

            if map_chrom in mask.contigs:
                for mask_rec in mask.fetch(map_chrom, orig_map_start, orig_map_end):
                    mask_start, mask_end = map(int, mask_rec.split()[1:3])
                    if orig_map_start >= mask_start and orig_map_end <= mask_end:
                        filters.append('Masked')

        if map_qual == 0:
            filters.append('ZeroMapQ')

    if len(maplocs) > 1:
        filters.append('Multimap')

    if len(filters) == 0:
        filters.append('PASS')

    return filters


def double_trd_filter(map_5p, map_3p, window=10000):
    map_5p_chrom, map_5p_start, map_5p_end, map_5p_qual = map_5p.split(':')
    map_3p_chrom, map_3p_start, map_3p_end, map_3p_qual = map_3p.split(':')

    map_5p_start = int(map_5p_start)-window
    map_5p_end   = int(map_5p_end)+window
    map_3p_start = int(map_3p_start)-window
    map_3p_end   = int(map_3p_end)+window

    if map_5p_chrom == map_3p_chrom:
        if min(map_5p_end, map_3p_end) - max(map_5p_start, map_3p_start) > 0:
            return False

    return True


def main(args):
    tldr_forest = dd(Intersecter)
    telib = load_falib(args.elts)

    mask = None

    if args.mask:
        mask = pysam.Tabixfile(args.mask)

    header = []

    trd_fa = 'trd.' + str(uuid4()) + '.fa'

    with open(trd_fa, 'w') as trd_out:
        with open(args.table, 'r') as table:
            for i, line in enumerate(table):
                if i == 0:
                    header = line.strip().split('\t')
                    if 'Transduction_5p' not in header or 'Transduction_3p' not in header:
                        logger.error('Transduction columns not present, please run tldr with --trdcol option.')
                        break

                    continue

                rec = {}

                for n, field in enumerate(line.strip().split('\t')):
                    rec[header[n]] = field

                if rec['Transduction_5p'] != 'NA' and len(rec['Transduction_5p']) > int(args.mintrd):
                    trd_out.write('>%s\n%s\n' % (rec['UUID']+'_5p', rec['Transduction_5p']))

                if rec['Transduction_3p'] != 'NA' and len(rec['Transduction_3p']) > int(args.mintrd):
                    trd_out.write('>%s\n%s\n' % (rec['UUID']+'_3p', rec['Transduction_3p']))

                tldr_forest[rec['Chrom']].add_interval(Interval(int(rec['Start']), int(rec['End']), value=rec))

    result = mm2_search(args.ref, trd_fa, args.refelts, args.nonrefelts, tldr_forest, window=int(args.window))

    with open(args.table, 'r') as table:
        for i, line in enumerate(table):
            if i == 0:
                header = line.strip().split('\t')
                header_out = line.strip()

                for s in ['5p','3p']:
                    header_out += '\tTransduction_%s_Maploc' % s
                    header_out += '\tTransduction_%s_TE' % s
                    header_out += '\tTransduction_%s_Filter' % s

                print(header_out)
                continue

            rec = {}

            for n, field in enumerate(line.strip().split('\t')):
                rec[header[n]] = field

            trd_map  = {'5p':'NA', '3p':'NA'}
            trd_te   = {'5p':'NA', '3p':'NA'}
            trd_filt = {'5p':['NA'], '3p':['NA']}
            
            for s in ['5p','3p']:
                if rec['UUID']+'_'+s in result:
                    maplocs = []
                    telocs = []

                    for trd in result[rec['UUID']+'_'+s]:
                        if trd.split(':')[-1] in ('REF', 'NR', 'TLDR'):
                            telocs.append(trd)
                        else:
                            maplocs.append(trd)

                        trd_map[s] = ','.join(maplocs)

                        if len(telocs) > 0:
                            trd_te[s]  = ','.join(telocs)

                        filters = filter(rec, telocs, maplocs, mask)

                        trd_filt[s] = filters

            if trd_filt['5p'][0] == trd_filt['3p'][0] == 'PASS':
                if double_trd_filter(trd_map['5p'], trd_map['3p']):
                    trd_filt['5p'] = ['LocDisagree']
                    trd_filt['3p'] = ['LocDisagree']

            match = 'NA'

            for s in ['5p', '3p']:
                if trd_filt[s][0] == 'PASS':
                    tr_seq = rec['Transduction_'+s].upper()
                    te_seq = telib[rec['Family'] + ':' + rec['Subfamily']].upper().rstrip('A')
                    match  = align(tr_seq, te_seq, elt=rec['UUID'])

            rec_out = line.strip()

            for s in ['5p', '3p']:
                rec_out += '\t%s' % trd_map[s]
                rec_out += '\t%s' % trd_te[s]
                rec_out += '\t%s' % ','.join(trd_filt[s])

            print(rec_out)

    os.remove(trd_fa)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='foo')
    parser.add_argument('-t', '--table', required=True, help='tldr table')
    parser.add_argument('-f', '--ref', required=True, help='reference genome (fasta or minimap2 index)')
    parser.add_argument('-e', '--elts', required=True, help='reference elements .fa')
    parser.add_argument('-r', '--refelts', required=True, help='reference TE locations (tabix-indexed)')
    parser.add_argument('-n', '--nonrefelts', default=None, help='known non-reference elements (tabix-indexed)')
    parser.add_argument('-w', '--window', default=100, help='search window (default = 100)')
    parser.add_argument('-m', '--mintrd', default=30, help='minimum transduction size (default=30)')
    parser.add_argument('--mask', default=None, help='tabix')
    args = parser.parse_args()
    main(args)

