#!/usr/bin/env python3

from collections import defaultdict as dd
from collections import Counter

import os
import pysam
import argparse

import pandas as pd
import numpy as np
import scipy.stats as ss

import gzip


class Read:
    def __init__(self, read_name, cpg_loc, llr, cutoff=2.5, phase=None):
        self.read_name = read_name
        self.llrs = {}
        self.meth_calls = {}
        self.phase = phase

        self.add_cpg(cpg_loc, llr, cutoff=cutoff)

    def add_cpg(self, cpg_loc, llr, cutoff = 2.5):
        #assert abs(llr) > cutoff

        self.llrs[cpg_loc] = llr

        if llr > cutoff:
            self.meth_calls[cpg_loc] = 1
        elif llr < -1*cutoff:
            self.meth_calls[cpg_loc] = -1
        else:
            self.meth_calls[cpg_loc] = 0


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
                seq   = ''
            else:
                assert seqid != ''
                seq = seq + line.strip()

    if seqid not in seqdict and seq != '':
        seqdict[seqid] = seq

    return seqdict


def get_reads(fn, min_mapq=10, tag_untagged=False, ignore_tags=False):
    reads = {}

    bam = pysam.AlignmentFile(fn)
    for read in bam.fetch():
        if read.mapq >= min_mapq:    
            phase = None

            if tag_untagged or ignore_tags:
                phase = 'unphased'

            HP = None
            PS = None

            if not ignore_tags:
                for tag in read.get_tags():
                    if tag[0] == 'HP':
                        HP = tag[1]
                    if tag[0] == 'PS':
                        PS = tag[1]

            if None not in (HP, PS):
                phase = str(PS) + ':' + str(HP)

            reads[read.query_name] = phase

    return reads


def sorted_unmapped_segments(seq):
    segs = []

    current_seg_start = 0
    current_seg_end = 0

    in_seg = False

    for i, b in enumerate(list(seq)):
        if b.islower():
            if not in_seg:
                in_seg = True
                current_seg_start = i

        if b.isupper():
            if in_seg:
                in_seg = False
                current_seg_end = i

                segs.append([current_seg_start, current_seg_end])

    if list(seq)[-1].islower():
        current_seg_end = len(seq)

        segs.append([current_seg_start, current_seg_end])

    if len(segs) == 0:
        return [[0,0]]

    return sorted(segs, key=lambda x: x[1]-x[0], reverse=True)


def main(args):
    tldr_table = pd.read_csv(args.table, sep='\t', header=0, index_col=0)

    tldr_dir = args.tldr_dir

    if tldr_dir is None:
        tldr_dir = '.'.join(args.table.split('.')[:-2])

    samples = []

    for sample_calls in tldr_table['SampleReads']:
        for sample_c in sample_calls.split(','):
            samples.append(sample_c.split('|')[0])

    samples = sorted(list(set(samples)))

    meth_table = dd(dict)
    meth_output = dd(dict)

    assert os.path.exists(tldr_dir)


    for uuid in tldr_table.index:
        ins = tldr_table.loc[uuid]

        cons_fa = tldr_dir + '/' + uuid + '.cons.ref.fa'

        if not os.path.exists(cons_fa):
            continue

        cons_dict = load_falib(cons_fa)
        cons_seq = cons_dict[uuid]

        for sample in samples:
            bam_fn = tldr_dir + '/' + sample + '.' + uuid + '.te.bam'
            meth_fn = tldr_dir + '/' + sample + '.' + uuid + '.te.meth.tsv.gz'

            meth_table[uuid][sample] = [0,0,0] # meth, unmeth, no_call

            if not os.path.exists(bam_fn):
                continue

            if not os.path.exists(meth_fn):
                continue

            chrom = uuid

            h_start, h_end = sorted_unmapped_segments(cons_seq)[0]  # defines TE start / end positions in contig

            # get relevant genome chunk to tmp tsv

            meth_tbx = pysam.Tabixfile(meth_fn)

            tmp_methdata = uuid + '.' + sample + '.tmp.methdata.tsv'

            with open(tmp_methdata, 'w') as meth_out:
                # header
                with gzip.open(meth_fn, 'rt') as _:
                    for line in _:
                        assert line.startswith('chromosome')
                        meth_out.write(line)
                        break

                assert chrom in meth_tbx.contigs

                for rec in meth_tbx.fetch(chrom, h_start, h_end):
                    meth_out.write(str(rec)+'\n')

            # index by read_name
            methdata = pd.read_csv(tmp_methdata, sep='\t', header=0, index_col=4)

            if not args.keep_tmp_table:
                os.remove(tmp_methdata)

            # get list of relevant reads
            reads = get_reads(bam_fn, tag_untagged=args.tag_untagged, ignore_tags=args.ignore_tags)

            readnames = []
            for r in reads.keys():
                if r in methdata.index:
                    readnames.append(r)

            methdata = methdata.loc[readnames]

            methreads = {}

            for index, row in methdata.iterrows():
                r_start = row['start']
                r_end   = row['end']
                llr     = row['log_lik_ratio']
                seq     = row['sequence']

                # get per-CG position (nanopolish/calculate_methylation_frequency.py)
                cg_pos = seq.find("CG")
                first_cg_pos = cg_pos
                while cg_pos != -1:
                    cg_start = r_start + cg_pos - first_cg_pos
                    cg_pos = seq.find("CG", cg_pos + 1)

                    cg_elt_start = cg_start - h_start

                    if cg_start >= h_start and cg_start <= h_end:
                        #print (cg_start, cg_elt_start, llr, index)
                        if index not in methreads:
                            methreads[index] = Read(index, cg_elt_start, llr, phase=reads[index], cutoff=float(args.cutoff))
                        else:
                            methreads[index].add_cpg(cg_elt_start, llr, cutoff=float(args.cutoff))

            for name, read in methreads.items():
                for loc in read.llrs.keys():
                    if read.meth_calls[loc] == 1:
                        meth_table[uuid][sample][0] += 1

                    if read.meth_calls[loc] == -1:
                        meth_table[uuid][sample][1] += 1

                    if read.meth_calls[loc] == 0:
                        meth_table[uuid][sample][2] += 1


        meth_output[uuid]['seg_id']     = uuid
        meth_output[uuid]['seg_chrom']  = ins['Chrom']
        meth_output[uuid]['seg_start']  = ins['Start']
        meth_output[uuid]['seg_end']    = ins['End']
        meth_output[uuid]['seg_name']   = ins['Subfamily']
        meth_output[uuid]['seg_strand'] = ins['Strand']

        for sample in samples:
            meth_output[uuid][sample + '_meth_calls']   = meth_table[uuid][sample][0]
            meth_output[uuid][sample + '_unmeth_calls'] = meth_table[uuid][sample][1]
            meth_output[uuid][sample + '_no_calls']     = meth_table[uuid][sample][2]

    col_order = ['seg_id', 'seg_chrom', 'seg_start', 'seg_end', 'seg_name', 'seg_strand']

    for sample in samples:
        col_order += [sample + '_meth_calls', sample + '_unmeth_calls', sample + '_no_calls']

    meth_output = pd.DataFrame.from_dict(meth_output).T
    meth_output = meth_output[col_order]

    out_fn = '.'.join(args.table.split('.')[:-1]) + '.nr.segmeth.table.txt'

    meth_output.to_csv(out_fn, sep='\t', index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='giant bucket')
    parser.add_argument('-t', '--table', required=True, help='tldr table')
    parser.add_argument('-d', '--tldr_dir', default=None)
    parser.add_argument('-c', '--cutoff', default=2.5, help='llr cutoff (absolute value), default=2.5')
    parser.add_argument('--keep_tmp_table', action='store_true', default=False)
    parser.add_argument('--ignore_tags', action='store_true', default=True)
    parser.add_argument('--tag_untagged', action='store_true', default=False)

    args = parser.parse_args()
    main(args)
