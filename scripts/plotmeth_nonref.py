#!/usr/bin/env python3

from __future__ import print_function

from collections import defaultdict as dd
from collections import Counter

import os
import pysam
import argparse

from operator import itemgetter

import pandas as pd
import numpy as np
import scipy.stats as ss

import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

#sns.set_palette('viridis', n_colors=2)

from matplotlib import gridspec
from matplotlib.patches import ConnectionPatch

from uuid import uuid4
import gzip

from statsmodels.nonparametric.smoothers_lowess import lowess


class Gene:
    def __init__(self, ensg, name):
        self.ensg = ensg
        self.name = name
        self.tx_start = None
        self.tx_end = None
        self.cds_start = None
        self.cds_end = None
        self.exons = []

    def add_exon(self, block):
        assert len(block) == 2
        assert block[0] < block[1]
        self.exons.append(block)
        self.exons = sorted(self.exons, key=itemgetter(0))

    def add_tx(self, block):
        assert len(block) == 2
        assert block[0] < block[1]
        if self.tx_start is None or self.tx_start > block[0]:
            self.tx_start = block[0]

        if self.tx_end is None or self.tx_end < block[1]:
            self.tx_end = block[1]

    def add_cds(self, block):
        assert len(block) == 2
        assert block[0] < block[1]

        if self.cds_start is None or self.cds_start > block[0]:
            self.cds_start = block[0]

        if self.cds_end is None or self.cds_end < block[1]:
            self.cds_end = block[1]

    def has_tx(self):
        return None not in (self.tx_start, self.tx_end)

    def has_cds(self):
        return None not in (self.cds_start, self.cds_end)

    def merge_exons(self):
        new_exons = []
        if len(self.exons) == 0:
            return

        last_block = self.exons[0]

        for block in self.exons[1:]:
            if min(block[1], last_block[1]) - max(block[0], last_block[0]) > 0: # overlap
                last_block = [min(block[0], last_block[0]), max(block[1], last_block[1])]

            else:
                new_exons.append(last_block)

            last_block = block

        new_exons.append(last_block)

        self.exons = new_exons


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


def slide_window(meth_table, sample, width=20, slide=2):
    midpt_min = min(meth_table['loc'])
    midpt_max = max(meth_table['loc'])

    sample_table = meth_table.loc[meth_table['sample'] == sample]

    win_start = int(midpt_min - width/2)
    win_end = win_start + width

    meth_frac = {}
    meth_n = {}

    while int((win_start+win_end)/2) < midpt_max:
        win_start += slide
        win_end += slide

        meth_count = len(meth_table.loc[(meth_table['sample'] == sample) & (meth_table['loc'] > win_start) & (meth_table['loc'] < win_end) & (meth_table['call'] == 1)])
        unmeth_count = len(meth_table.loc[(meth_table['sample'] == sample) & (meth_table['loc'] > win_start) & (meth_table['loc'] < win_end) & (meth_table['call'] == -1)])

        midpt = int((win_start+win_end)/2)

        if meth_count + unmeth_count > 0:
            meth_frac[midpt] = meth_count/(meth_count+unmeth_count)
            meth_n[midpt] = meth_count+unmeth_count

    return meth_frac, meth_n


def smooth(x, window_len=8, window='hanning'):
    ''' modified from scipy cookbook: https://scipy-cookbook.readthedocs.io/items/SignalSmooth.html '''

    assert window_len % 2 == 0, '--smoothwindowsize must be an even number'
    assert x.ndim == 1
    assert x.size > window_len

    if window_len<3:
        return x

    assert window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']

    s=np.r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]]
    
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='valid')

    return y[(int(window_len/2)-1):-(int(window_len/2))]


def mask_methfrac(data, cutoff=20):
    data = np.asarray(data)
    data = data > int(cutoff)

    segs = []

    in_seg = False
    seg_start = 0

    for i in range(len(data)):
        if data[i]:
            if in_seg:
                segs.append(list(range(seg_start, i)))

            in_seg = False

        else:
            if not in_seg:
                seg_start = i

            in_seg = True

    if in_seg:
        segs.append(list(range(seg_start, len(data))))

    return segs


def build_genes(gtf, chrom, start, end):
    genes = {}

    for line in gtf.fetch(chrom, start, end):

        chrom, source, feature, start, end, score, strand, frame, attribs = line.split('\t')

        block = [int(start), int(end)]

        attribs = attribs.strip()

        attr_dict = {}

        for attrib in attribs.split(';'):
            if attrib:
                key, val = attrib.strip().split()[:2]
                key = key.strip()
                val = val.strip().strip('"')
                attr_dict[key] = val

        if 'gene_id' not in attr_dict:
            continue

        if 'gene_name' not in attr_dict:
            continue

        ensg = attr_dict['gene_id']
        name = attr_dict['gene_name']

        if ensg not in genes:
            genes[ensg] = Gene(ensg, name)

        if feature == 'exon':
            genes[ensg].add_exon(block)

        if feature == 'CDS':
            genes[ensg].add_cds(block)

        if feature == 'transcript':
            genes[ensg].add_tx(block)

    return genes


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
    teont_table = pd.read_csv(args.table, sep='\t', header=0, index_col=0)

    teont_dir = None
    if args.teont_dir is None:
        teont_dir = '.'.join(args.table.split('.')[:-2])

    cons_fa = teont_dir + '/' + args.uuid + '.cons.ref.fa'
    cons_dict = load_falib(cons_fa)
    cons_seq = cons_dict[args.uuid]

    meth_table = dd(dict)
    sample_order = []


    for sample in args.sample.split(','):

        sample_order.append(sample)

        bam_fn = teont_dir + '/' + sample + '.' + args.uuid + '.te.bam'
        meth_fn = teont_dir + '/' + sample + '.' + args.uuid + '.te.meth.tsv.gz'

        assert os.path.exists(teont_dir)
        assert os.path.exists(bam_fn), 'not found: %s' % bam_fn
        assert os.path.exists(meth_fn), 'not found: %s' % meth_fn
        assert args.uuid in teont_table.index

        ins = teont_table.loc[args.uuid]

        chrom = args.uuid
        elt_start = 0
        elt_end = len(cons_seq)

        fn_prefix = '.'.join(('_'.join(args.sample.split(',')), args.uuid, ins['Chrom'], str(ins['Start']), ins['Subfamily']))

        h_start, h_end = sorted_unmapped_segments(cons_seq)[0]  # defines TE start / end positions in contig
        h_cpg_start = None
        h_cpg_end = None

        # get relevant genome chunk to tmp tsv

        meth_tbx = pysam.Tabixfile(meth_fn)

        tmp_methdata = fn_prefix+'.tmp.methdata.tsv'

        with open(tmp_methdata, 'w') as meth_out:
            # header
            with gzip.open(meth_fn, 'rt') as _:
                for line in _:
                    assert line.startswith('chromosome')
                    meth_out.write(line)
                    break

            assert chrom in meth_tbx.contigs

            for rec in meth_tbx.fetch(chrom, elt_start, elt_end):
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

                cg_elt_start = cg_start - elt_start

                if cg_start >= elt_start and cg_start <= elt_end:
                    #print (cg_start, cg_elt_start, llr, index)
                    if index not in methreads:
                        methreads[index] = Read(index, cg_elt_start, llr, phase=reads[index], cutoff=float(args.cutoff))
                    else:
                        methreads[index].add_cpg(cg_elt_start, llr, cutoff=float(args.cutoff))

        
        #meth_table = dd(dict)

        #sample_order = []

        for name, read in methreads.items():

            for loc in read.llrs.keys():
                uuid = str(uuid4())
                meth_table[uuid]['loc'] = loc
                meth_table[uuid]['llr'] = read.llrs[loc]
                meth_table[uuid]['read'] = name
                meth_table[uuid]['sample'] = sample
                meth_table[uuid]['call'] = read.meth_calls[loc]



    # table for plotting
    meth_table = pd.DataFrame.from_dict(meth_table).T
    meth_table['loc'] = pd.to_numeric(meth_table['loc'])
    meth_table['llr'] = pd.to_numeric(meth_table['llr'])

    # cpg space
    meth_table['orig_loc'] = meth_table['loc']
    meth_table['loc'] = ss.rankdata(meth_table['loc'], method='dense')

    coord_to_cpg = {}
    for orig_loc, new_loc in zip(meth_table['orig_loc'], meth_table['loc']):
        coord_to_cpg[orig_loc] = new_loc


    h_cpg_start = coord_to_cpg[min(meth_table['orig_loc'], key=lambda x:abs(x-h_start))]
    h_cpg_end = coord_to_cpg[min(meth_table['orig_loc'], key=lambda x:abs(x-h_end))]

    fig = plt.figure()

    gs_ratio = [0.01,1,3,3]

    if args.gtf is not None:
        gs_ratio = [1,1,3,3]

    gs = gridspec.GridSpec(4,1,height_ratios=gs_ratio)


    # plot genes
    ax0 = plt.subplot(gs[0])

    ax0.spines['bottom'].set_visible(False)
    ax0.spines['left'].set_visible(False)
    ax0.spines['right'].set_visible(False)
    ax0.xaxis.set_ticks_position('top')

    gtf = None

    if args.gtf is not None:
        gtf = pysam.Tabixfile(args.gtf)

    genes = []
    if gtf is not None:
        genes = build_genes(gtf, ins['Chrom'], int(ins['Start'])-int(args.windowsize), int(ins['End'])+int(args.windowsize))

    exon_patches = []
    tx_lines = []

    genes_of_interest = []

    if args.genes is not None:
        genes_of_interest = args.genes.strip().split(',')

    i = 0
    for ensg in genes:
        window_offset = int(ins['Start'])-int(args.windowsize)

        if genes_of_interest:
            if genes[ensg].name not in genes_of_interest:
                continue

        if genes[ensg].has_tx():
            tx_lines.append(matplotlib.lines.Line2D([genes[ensg].tx_start-window_offset, genes[ensg].tx_end-window_offset], [0.4+i, 0.4+i], zorder=1))

            print('transcript: %d-%d %s' % (genes[ensg].tx_start, genes[ensg].tx_end, genes[ensg].name))

        genes[ensg].merge_exons()
        for exon_start, exon_end in genes[ensg].exons:
            exon_len = exon_end - exon_start

            exon_patches.append(matplotlib.patches.Rectangle([exon_start-window_offset, i], exon_len, 0.8, edgecolor='#777777', facecolor='#ff4500', zorder=2))


        blocks_str = ','.join(['%d-%d' % (s,e) for s, e in genes[ensg].exons])
        print('%s exons: %s' % (genes[ensg].name, blocks_str))

        i += 1

    if i < 3:
        i = 3

    ax0.set_ylim(0,i)
    ax0.set_yticks([])

    for p in exon_patches:
        ax0.add_patch(p)

    for tx in tx_lines:
        ax0.add_line(tx)


    # plot correspondence between genome space and cpg space
    ax1 = plt.subplot(gs[1])
    ax2 = ax1.twiny()

    # set window size a close to window as possible (nearest CpG locations)

    view_orig_start = 0
    view_orig_end = 0
    view_cpg_start = 0
    view_cpg_end = 0

    for oc in coord_to_cpg.keys():
        if oc > h_start - int(args.windowsize):
            view_orig_start = oc
            break

    for oi, oc in enumerate(coord_to_cpg.keys()):
        if oc > h_end + int(args.windowsize):
            view_orig_end = list(coord_to_cpg.keys())[oi-1]
            break

    view_cpg_start = coord_to_cpg[view_orig_start]
    view_cpg_end = coord_to_cpg[view_orig_end]

    ax2.set_xlim(view_orig_start, view_orig_end)
    ax1.set_xlim(view_cpg_start, view_cpg_end)

    ax1.set_ylim(0,10)
    ax1.set_yticklabels([])


    x1 = []
    x2 = []

    step = int(args.topspacing)

    for i, x in enumerate(meth_table['orig_loc']):
        if i in (0, len(meth_table['orig_loc'])-1):
            x2.append(x)
            x1.append(coord_to_cpg[x])

        elif i % step == 0:
            x2.append(x)
            x1.append(coord_to_cpg[x])

    
    ax1.vlines(x1, 0, 1, color='#777777', zorder=1)
    ax2.vlines(x2, 9, 10, color='#777777', zorder=1)


    orig_highlight_box = matplotlib.patches.Rectangle((h_start,9), h_end-h_start, 1.0, lw=1, edgecolor='#777777', facecolor='#a8caff', zorder=2)
    cpg_highlight_box = matplotlib.patches.Rectangle((h_cpg_start,0), h_cpg_end-h_cpg_start, 1.0, lw=1, edgecolor='#777777', facecolor='#a8caff', zorder=3)

    # arrows

    orig_arr_x = h_start
    orig_arr_y = 9.5
    orig_arr_dx = (h_end-h_start) - ((h_end-h_start)*0.03)
    orig_arr_dy = 0

    cpg_arr_x = h_cpg_start
    cpg_arr_y = 0.5
    cpg_arr_dx = (h_cpg_end-h_cpg_start) - ((h_cpg_end-h_cpg_start)*0.03)
    cpg_arr_dy = 0

    if ins['Strand'] == '-':
        orig_arr_x = h_end
        orig_arr_dx = 0 - orig_arr_dx

        cpg_arr_x = h_cpg_end
        cpg_arr_dx = 0- cpg_arr_dx

    orig_arrow = ax2.arrow(orig_arr_x, orig_arr_y, orig_arr_dx, orig_arr_dy, head_width=1, head_length=(h_end-h_start)*0.03, zorder=4)
    cpg_arrow = ax1.arrow(cpg_arr_x, cpg_arr_y, cpg_arr_dx, cpg_arr_dy, head_width=1, head_length=(h_cpg_end-h_cpg_start)*0.03, zorder=4)


    ax2.add_patch(orig_highlight_box)
    ax1.add_patch(cpg_highlight_box)

    for x1_x, x2_x in zip(x1, x2):
        link_end1 = (x1_x, 1)
        link_end2 = (x2_x, 9)

        l_col = '#777777'

        if x2_x >= h_start and x2_x <= h_end:
            l_col = '#3080ff'

        con = ConnectionPatch(xyA=link_end1, xyB=link_end2, coordsA="data", coordsB="data", axesA=ax1, axesB=ax2, color=l_col)
        ax2.add_artist(con)

    ax0.set_xlim(ax2.get_xlim()) # sync axes between orig coords and gtf plot
    ax2.set_xticks([])

    xt_labels = [str(int(t+elt_start)) for t in list(ax0.get_xticks())]
    ax0.set_xticklabels(xt_labels)


    # llr plot

    ax3 = plt.subplot(gs[2])

    ax3.axhline(y=float(args.cutoff), c='k', linestyle='--',lw=1)
    ax3.axhline(y=0, c='#bbbbbb', linestyle='--',lw=1)
    ax3.axhline(y=0-float(args.cutoff), c='k', linestyle='--',lw=1)

    ax3 = sns.lineplot(x='loc', y='llr', hue='sample', data=meth_table)

    sample_color = {}
    for i, sample in enumerate(sample_order):
        sample_color[sample] = sns.color_palette(n_colors=len(sample_order))[i]

    ax3.set_xlim(ax1.get_xlim())

    # meth frac plot

    ax5 = plt.subplot(gs[3])

    for sample in sample_order:
        windowed_methfrac, meth_n = slide_window(meth_table, sample, width=int(args.slidingwindowsize), slide=int(args.slidingwindowstep))
        
        smoothed_methfrac = smooth(np.asarray(list(windowed_methfrac.values())), window_len=int(args.smoothwindowsize))

        #print(sample, ','.join(map(str, list(meth_n.values()))))

        masked_segs = mask_methfrac(list(meth_n.values()))
        #print(sample, masked_segs)

        ax5.plot(list(windowed_methfrac.keys()), smoothed_methfrac, marker='', color=sample_color[sample], zorder=1)

        for seg in masked_segs:
            if len(seg) > 2:
                mf_seg = np.asarray(smoothed_methfrac)[seg]
                pos_seg = np.asarray(list(windowed_methfrac.keys()))[seg]
            
                ax5.plot(pos_seg, mf_seg, marker='', color='#ffffff', zorder=2)


    ax5.set_xlim(ax1.get_xlim())
    ax5.set_ylim((-0.05,1.05))

    fig.set_size_inches(16, 8)

    imgtype = 'png'

    if args.svg:
        imgtype = 'svg'

    if args.ignore_tags:
        plt.savefig('%s.unphased.meth.%s' % (fn_prefix, imgtype), bbox_inches='tight')
    else:
        plt.savefig('%s.phased.meth.%s' % (fn_prefix, imgtype), bbox_inches='tight')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='giant bucket')
    parser.add_argument('-g', '--gtf', default=None, help='genes or intervals to display in gtf format')
    parser.add_argument('-t', '--table', required=True, help='teont table')
    parser.add_argument('-d', '--teont_dir', default=None)
    parser.add_argument('-u', '--uuid', required=True, help='target UUID')
    parser.add_argument('-s', '--sample', required=True)
    parser.add_argument('-c', '--cutoff', default=2.5, help='llr cutoff (absolute value), default=2.5')
    parser.add_argument('-w', '--windowsize', required=True)
    parser.add_argument('--slidingwindowsize', default=20, help='size of sliding window for meth frac (default=20)')
    parser.add_argument('--slidingwindowstep', default=2, help='step size for meth frac (default=2)')
    parser.add_argument('--smoothwindowsize', default=8, help='size of window for smoothing (default=8)')
    parser.add_argument('--methcall_ymax', default=None)
    parser.add_argument('--topspacing', default=10, help='spacing between links in top panel (default=10)')
    parser.add_argument('--genes', default=None, help='genes of interest (comma delimited)')
    parser.add_argument('--keep_tmp_table', action='store_true', default=False)
    parser.add_argument('--excl_ambig', action='store_true', default=False)
    parser.add_argument('--ignore_tags', action='store_true', default=True)
    parser.add_argument('--tag_untagged', action='store_true', default=False)
    parser.add_argument('--skip_callplot', action='store_true', default=False)
    parser.add_argument('--skip_wiggle', action='store_true', default=False)
    parser.add_argument('--svg', action='store_true', default=False)


    args = parser.parse_args()
    main(args)
