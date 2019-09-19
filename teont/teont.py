#!/usr/bin/env python3

import os
import sys
import pysam
import subprocess
import argparse
import pickle
import multiprocessing as mp
import numpy as np
#import scipy.stats as ss

from uuid import uuid4
from collections import defaultdict as dd
from collections import OrderedDict as od
from collections import Counter
from bx.intervals.intersection import Intersecter, Interval

from sklearn import metrics
from sklearn.mixture import GaussianMixture

import logging
FORMAT = '%(asctime)s %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


class AlignedColumn:
    def __init__(self):
        self.bases = od() # species name --> base
        self.annotations = od() # metadata

    def gap(self):
        if '-' in self.bases.values():
            return True
        return False

    def subst(self):
        if self.gap():
            return False

        if len(set(self.bases.values())) > 1:
            return True
        return False

    def cons(self, min_override_gap=2):
        ''' weighted consensus '''
        most_common = Counter(map(str.upper, self.bases.values())).most_common()[0]

        if most_common[0] != '-':
            return most_common

        if len(Counter(map(str.upper, self.bases.values())).most_common()) > 1:
            if Counter(map(str.upper, self.bases.values())).most_common()[1][1] >= min_override_gap:
                return Counter(map(str.upper, self.bases.values())).most_common()[1]

        return most_common

    def __str__(self):
        return str(self.bases)


class MSA:
    ''' multiple sequence alignment class '''
    def __init__(self, infile=None):
        self.columns = []
        self.ids     = []
        self.seqs    = od()

        if infile is not None:
            self.readFastaMSA(infile)

    def __len__(self):
        return len(self.columns)

    def readFastaMSA(self, infile):
        id   = None
        seq  = ''

        with open(infile, 'r') as fasta:
            for line in fasta:
                line = line.strip()
                if line.startswith('>'):
                    if id is not None:
                        self.seqs[id] = seq
                    seq = ''
                    id = line.lstrip('>')
                    self.ids.append(id)
                else:
                    seq += line
            self.seqs[id] = seq

        first = True
        colen = 0
        for ID, seq in self.seqs.items():
            if first:
                colen = len(seq)
                for base in list(seq):
                    ac = AlignedColumn()
                    ac.bases[ID] = base
                    self.columns.append(ac)
                first = False
            else:
                assert len(seq) == colen
                pos = 0
                for base in list(seq):
                    ac = self.columns[pos]
                    ac.bases[ID] = base
                    pos += 1

    def consensus(self, gapped=False):
        ''' compute consensus '''
        bases = []
        for column in self.columns:
            b = column.cons()[0]
            if gapped or b != '-':
                bases.append(b)
        return ''.join(bases)


class InsRead:
    def __init__(self, bamfile, q_chrom, q_start, q_end, r_start, r_end, r_name, r_seq, r_qual, is_ins, is_clip, clip_end):
        self.bamname  = '.'.join(os.path.basename(bamfile).split('.')[:-1])
        self.q_chrom  = q_chrom
        self.q_start  = q_start
        self.q_end    = q_end
        self.r_start  = r_start
        self.r_end    = r_end
        self.r_name   = r_name
        self.r_seq    = r_seq
        self.r_qual   = r_qual
        self.is_ins   = is_ins
        self.is_clip  = is_clip
        self.clip_end = clip_end

        self.r_seq_trimmed   = None
        self.r_qual_trimmed  = None
        self.q_start_trimmed = None
        self.q_end_trimmed   = None

        # te alignment info
        self.te             = None
        self.te_fam         = None
        self.te_subfam      = None
        self.te_elt_start   = None
        self.te_elt_end     = None
        self.te_r_seq_start = None
        self.te_r_seq_end   = None
        self.te_read_orient = None

        self.useable = False # set by process_cluster

        self.r_pos = None
        if self.r_start:
            self.r_pos = self.r_start
        else:
            self.r_pos = self.r_end

    def te_info(self, te):
        self.te = te
        self.te_fam, self.te_subfam = te.subfamily().split(':')
        self.te_elt_start   = te.min_qry_coord()
        self.te_elt_end     = te.max_qry_coord()
        self.te_r_seq_start = te.min_tgt_coord()
        self.te_r_seq_end   = te.max_tgt_coord()
        self.te_read_orient = te.max_qry_orient()

        #if self.te_read_orient == '-':
        #    self.te_r_seq_start, self.te_r_seq_end = self.te_r_seq_end, self.te_r_seq_start

    def fastq(self):
        return '@%s\n%s\n+\n%s\n' % (self.r_name, self.r_seq_trimmed, self.r_qual_trimmed)

    def fasta(self):
        type = 'ins'
        if self.is_clip:
            type = 'clip_%s' % self.clip_end

        return '>%s_%s\n%s\n' % (self.r_name, type, self.r_seq_trimmed) 

    def trim(self, flanksize=200):

        trim_start = self.q_start
        trim_end = self.q_end

        if trim_start < flanksize:
            trim_start = 0
            self.q_start_trimmed = self.q_start

        else:
            trim_start -= flanksize
            self.q_start_trimmed = flanksize

        if trim_end + flanksize > len(self.r_seq):
            trim_end = len(self.r_seq)

        else:
            trim_end += flanksize

        self.q_end_trimmed = self.q_start_trimmed + (self.q_end - self.q_start)
        self.r_seq_trimmed = self.r_seq[trim_start:trim_end]
        self.r_qual_trimmed = self.r_qual[trim_start:trim_end]

    def te_overlap_frac(self):
        if None in (self.q_end_trimmed, self.te_r_seq_end, self.q_start_trimmed, self.te_r_seq_start):
            return 0.0

        return float(self.te_r_seq_end - self.te_r_seq_start) / float(self.q_end_trimmed - self.q_start_trimmed)

    def __lt__(self, other):
        return self.r_pos < other.r_pos


class InsCluster:
    def __init__(self, uuid):
        self.reads  = []
        self.uuid   = uuid
        self.cons   = None

        self.cons_te_align = None
        self.cons_align_map = None

        self.reset_major_family = None
        self.reset_major_subfam = None

        self.breakpoints = []
        self.breakpoints_remappable = False
        self.tsd_seq = 'NA'

    def add(self, read, flanksize=200):
        self.reads.append(read)

    def chrom(self):
        return self.reads[0].q_chrom

    def interval(self):
        p = sorted([read.r_pos for read in self.reads if read.useable])
        return p[0], p[-1]

    def dump_fastq(self, te_fam=None, sample=None):
        out_fn = '%s.%s.cluster.fq' % (str(sample), self.uuid)

        with open(out_fn, 'w') as out:
            for read in self.reads:
                if not read.useable:
                    continue

                if sample is not None and read.bamname != sample:
                    continue

                if read.te: # require TE alignment
                    if te_fam:
                        if read.te_fam == te_fam:
                            out.write(read.fastq())
                    else:
                        out.write(read.fastq())
        return out_fn

    def dump_fasta(self, te_fam=None):
        out_fn = '%s.cluster.fa' % self.uuid

        with open(out_fn, 'w') as out:
            for read in self.reads:
                if not read.useable:
                    continue

                if read.te: # require TE alignment
                    if te_fam:
                        if read.te_fam == te_fam:
                            out.write(read.fasta())
                    else:
                        out.write(read.fasta())

        return out_fn

    def te_align_count(self, family=None):
        if family:
            return len([r for r in self.reads if r.te_fam == family])
        else:
            return len([r for r in self.reads if r.te])

    def te_major_family(self):
        if self.reset_major_family is None:
            return Counter([read.te_fam for read in self.reads]).most_common()[0][0]

        return self.reset_major_family

    def te_major_subfam(self):
        if self.reset_major_subfam is None:
            return Counter([read.te_subfam for read in self.reads]).most_common()[0][0]

        return self.reset_major_subfam

    def te_orientation(self):
        return Counter([read.te_read_orient for read in self.reads]).most_common()[0][0]

    def te_fam_frac(self, family):
        return float(len([r for r in self.reads if r.te_fam == family])) / float(len(self))

    def te_median_overlap(self, span_only=False):
        if span_only:
            return np.median([read.te_overlap_frac() for read in self.reads if read.useable and read.is_ins])

        return np.median([read.te_overlap_frac() for read in self.reads if read.useable])

    def te_median_length(self):
        l = []
        for read in self.reads:
            if read.useable and read.is_ins:
                l.append(read.te_elt_end - read.te_elt_start)

        return np.median(l)

    def te_useable_count(self):
        return len([read for read in self.reads if read.useable])

    def te_embedded_count(self):
        return len([read for read in self.reads if read.useable and read.is_ins])

    def te_samples(self):
        return len(set([read.bamname for read in self.reads if read.useable]))

    def te_sample_count(self):
        sample_count = Counter([read.bamname for read in self.reads if read.useable])
        return ','.join(['|'.join(map(str, k)) for k in sample_count.items()])

    def detect_inversion(self):
        if self.cons_te_align is None:
            return False

        return self.cons_te_align.max_qry_orient() != self.cons_te_align.min_qry_orient()

    def trim_reads(self, flanksize=200):
        for read in self.reads:
            read.trim(flanksize=flanksize)

    def adjust_cons(self, left_corr, right_corr):
        left_orig, right_orig = self.sorted_unmapped_segments()[0]

        cons = list(self.cons)

        if left_corr < 0:
            for b in range(left_orig, left_orig+abs(left_corr)):
                cons[b] = cons[b].upper()

        if left_corr > 0:
            for b in range(left_orig-left_corr, left_orig):
                cons[b] = cons[b].lower()

        if right_corr < 0:
            for b in range(right_orig, right_orig+abs(right_corr)):
                cons[b] = cons[b].lower()

        if right_corr > 0:
            for b in range(right_orig-right_corr, right_orig):
                cons[b] = cons[b].upper()

        if cons != self.cons:
            self.cons = ''.join(cons)
            return True

        return False

    def extend_tsd(self):
        if self.tsd_seq == 'NA':
            return 0,0

        tsd_len = len(self.tsd_seq)

        left_bp, right_bp = self.sorted_unmapped_segments()[0]

        tsd_left = self.cons[left_bp-tsd_len:left_bp]
        tsd_right = self.cons[right_bp:right_bp+tsd_len]

        cons = list(self.cons)

        # extend left
        l = left_bp-tsd_len-1
        r = right_bp-1

        l_ext = 0
        while l > 0 and cons[l].upper() == cons[r].upper():
            #print(cons[l], cons[r])
            cons[l] = cons[l].upper()
            cons[r] = cons[r].upper()
            self.tsd_seq = cons[l].upper() + self.tsd_seq
            l -= 1
            r -= 1
            l_ext += 1


        #print('left', l_ext)


        # extend right
        l = left_bp
        r = right_bp+tsd_len

        r_ext = 0
        while r < len(cons) and cons[l].upper() == cons[r].upper():
            #print(cons[l], cons[r])
            cons[l] = cons[l].upper()
            cons[r] = cons[r].upper()
            self.tsd_seq = self.tsd_seq + cons[l].upper()
            l += 1
            r += 1
            r_ext += 1


        #print('right', r_ext)

        #print(tsd_left, tsd_right)

        self.cons = ''.join(cons)

        return l_ext, r_ext


    def sorted_unmapped_segments(self):
        if self.cons_align_map is None:
            return [[0,0]]

        segs = []

        current_seg_start = 0
        current_seg_end = 0

        in_seg = False

        for i, b in enumerate(list(self.cons)):
            if b.islower():
                if not in_seg:
                    in_seg = True
                    current_seg_start = i

            if b.isupper():
                if in_seg:
                    in_seg = False
                    current_seg_end = i

                    segs.append([current_seg_start, current_seg_end])

        if list(self.cons)[-1].islower():
            current_seg_end = len(self.cons)

            segs.append([current_seg_start, current_seg_end])

        if len(segs) == 0:
            return [[0,0]]

        return sorted(segs, key=lambda x: x[1]-x[0], reverse=True)


    def join_unmapped_segments(self, exp_flank_len, telib_seq, minlen=100, minfrac=0.7):
        segs = self.sorted_unmapped_segments()
        assert len(segs) > 1

        if segs[1][1] - segs[1][0] < minlen:
            return None

        # close enough to the expected flank length?
        if segs[1][0] > exp_flank_len*1.5 and segs[1][1] < len(self.cons)-(exp_flank_len*1.5):
            return None

        seg_seq = self.cons[segs[1][0]:segs[1][1]]

        subalign = te_align(telib_seq, seg_seq, te_fa_is_seq=True)

        if subalign is None:
            return None

        if (subalign.max_qry_coord() - subalign.min_qry_coord()) / len(seg_seq) < minfrac:
            return None
        
        cons = list(self.cons)

        for p in range(min(segs[0][0], segs[1][0]), max(segs[0][1], segs[1][1])):
            cons[p] = cons[p].lower()

        self.cons = ''.join(cons)

        return self.cons


    def reset_breakpoints(self, map_start, map_end):
        ''' reset breakpoints based on refined te mapping information '''
        assert len(self.breakpoints) == 2
        if cons_align_map is None:
            return None


    def mask_cons(self, trim=True):
        if self.cons_align_map is None:
            return self.cons

        assert len(self.cons_align_map) == len(self.cons)

        masked_cons = []

        for i in range(len(self.cons_align_map)):
            if self.cons_align_map[i] == 0:
                masked_cons.append(self.cons[i].lower())

            if self.cons_align_map[i] == 1:
                masked_cons.append(self.cons[i].upper())

        if trim:
            # trim forward
            trimmed_cons = []
            trimming = True

            for b in masked_cons:
                if b.isupper():
                    trimming = False

                if not trimming:
                    trimmed_cons.append(b)

            masked_cons = trimmed_cons

            # trim reverse
            trimmed_cons = []
            trimming = True

            for b in masked_cons[::-1]:
                if b.isupper():
                    trimming = False

                if not trimming:
                    trimmed_cons.append(b)

            masked_cons = trimmed_cons[::-1]

        self.cons = ''.join(masked_cons)

        return self.cons


    def te_overlap_zscores(self):
        z = {}
        for read in self.reads:
            if read.te_overlap_frac() > 0.0:
                z[read.r_name] = read.te_overlap_frac()

        o_mean = np.mean(list(z.values()))
        o_std   = np.std(list(z.values()))

        for name, s in z.items():
            z[name] = abs(s-o_mean)/o_std

        return z

    def prelim_filter(self, minreads=4): # TODO passthrough params
        passed = True

        # minimum cluster size
        if len(self.reads) < minreads:
            passed = False

        # require at least one fully embedded insertion
        if len([r for r in self.reads if r.is_ins]) < 1:
            passed = False

        return passed

    def __len__(self):
        return len(self.reads)


class ExonerateAlignment:
    def __init__(self, alignment):
        c = alignment.strip().split('\t')
        self.original  = c
        self.subfamily = c[1]
        self.score     = int(c[2])
        self.qry_start = int(c[3])
        self.qry_end   = int(c[4])
        self.tgt_start = int(c[5])
        self.tgt_end   = int(c[6])
        self.match     = float(c[7])
        self.orient    = c[9]

        if self.orient == '-':
            self.tgt_start, self.tgt_end = self.tgt_end, self.tgt_start

    def __lt__(self, other):
        return self.qry_start < other.qry_start


class AlignmentGroup:
    def __init__(self, alignment):
        self.alignments = []
        self.add_alignment(alignment)

    def add_alignment(self, alignment):
        if len(self.alignments) == 0:
            self.alignments.append(alignment)

        else:
            non_redundant = True

            for a in self.alignments:
                if alignment.qry_start > a.qry_start and alignment.qry_end < a.qry_end:
                    non_redundant = False

            if non_redundant:    
                self.alignments.append(alignment)

    def subfamily(self):
        return self.alignments[0].subfamily

    def best_score(self):
        return max([a.score for a in self.alignments])

    def best_match(self):
        return max([a.match for a in self.alignments])

    def min_qry_coord(self):
        return min([a.qry_start for a in self.alignments])

    def max_qry_coord(self):
        return max([a.qry_end for a in self.alignments])

    def min_tgt_coord(self):
        return min([a.tgt_start for a in self.alignments])

    def max_tgt_coord(self):
        return max([a.tgt_end for a in self.alignments])

    def min_qry_orient(self):
        ''' return orientation of leftmost alignment '''
        return sorted(self.alignments)[0].orient

    def max_qry_orient(self):
        ''' return orientation of rightmost alignment '''
        return sorted(self.alignments)[-1].orient


def rc(dna):
    ''' reverse complement '''
    complements = str.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
    return dna.translate(complements)[::-1]


def mafft(in_fn, iterations=1000, threads=1):
    ''' use MAFFT to create MSA '''

    out_fn = '.'.join(in_fn.split('.')[:-1]) + '.msa.fa'

    args = ['mafft', '--maxiterate', str(iterations), '--thread', str(threads), in_fn]

    FNULL = open(os.devnull, 'w')

    p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=FNULL)

    with open(out_fn, 'w') as out_fa:
        for line in p.stdout:
            line = line.decode()
            out_fa.write(line)

    return out_fn


def te_align(te_fa, refseq, minmatch=80.0, te_fa_is_seq=False, max_segs=3):
    #print('te_align debug:')
    rnd = str(uuid4())
    tgt_fa = 'tmp.' + rnd + '.tgt.fa'

    with open(tgt_fa, 'w') as tgt:
        tgt.write('>ref\n%s\n' % refseq)

    if te_fa_is_seq:
        te_seq = te_fa
        te_fa = 'tmp.' + rnd + '.te.fa'
        with open(te_fa, 'w') as te:
            te.write('>te\n%s\n' % te_seq)

    cmd = ['exonerate', '--bestn', str(max_segs), '--model', 'affine:local', '--showalignment','0', '--ryo', 'TE\t%qi\t%s\t%qab\t%qae\t%tab\t%tae\t%pi\t%qS\t%tS\n', te_fa, tgt_fa]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    align_groups = {}

    for pline in p.stdout.readlines():
        pline = pline.decode()
        if pline.startswith('TE'):
            a = ExonerateAlignment(pline.strip())
            if a.match >= minmatch:
                if a.subfamily in align_groups:
                    align_groups[a.subfamily].add_alignment(a)

                else:
                    align_groups[a.subfamily] = AlignmentGroup(a)

    os.remove(tgt_fa)

    if te_fa_is_seq:
        os.remove(te_fa)

    best_group = None

    for subfam in align_groups:
        if best_group is None:
            best_group = best_group = align_groups[subfam]

        elif align_groups[subfam].best_score() > best_group.best_score():
            best_group = align_groups[subfam]

    return best_group


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


def minimap2_pileup(cluster, dump_fq, outbase):
    out_name = cluster.uuid
    ins_ref = cluster.cons

    tmp_ref = str(uuid4()) + '.fa'
    with open(tmp_ref, 'w') as out:
        out.write('>ins_ref\n%s\n' % ins_ref)

    filter_start, filter_end = cluster.sorted_unmapped_segments()[0]

    FNULL = open(os.devnull, 'w')

    mm2_cmd  = ['minimap2', '-x', 'map-ont', '-a', tmp_ref, dump_fq]
    view_cmd = ['samtools', 'view', '-b', '-']
    sort_cmd = ['samtools', 'sort', '-']
    pile_cmd = ['samtools', 'mpileup', '-B', '-Q', '1', '-f', tmp_ref, '-']
    aln  = subprocess.Popen(mm2_cmd, stdout=subprocess.PIPE, stderr=FNULL)
    view = subprocess.Popen(view_cmd, stdin=aln.stdout, stdout=subprocess.PIPE, stderr=FNULL)
    sort = subprocess.Popen(sort_cmd, stdin=view.stdout, stdout=subprocess.PIPE, stderr=FNULL)
    pile = subprocess.Popen(pile_cmd, stdin=sort.stdout, stdout=subprocess.PIPE, stderr=FNULL)

    pileup = {}
    depth = {}

    for line in pile.stdout:
        c = line.decode().split()
        pileup[int(c[1])] = parse_pileup(c[4])
        depth[int(c[1])] = int(c[3])

    altered_ref = alter_ref(ins_ref, pileup, depth)#, filter_start, filter_end)

    os.remove(tmp_ref)

    if os.path.exists(tmp_ref+'.fai'):
        os.remove(tmp_ref+'.fai')

    return altered_ref


def minimap2_bam(ins_ref, dump_fq, outbam):

    FNULL = open(os.devnull, 'w')

    mm2_cmd  = ['minimap2', '-x', 'map-ont', '-a', ins_ref, dump_fq]
    view_cmd = ['samtools', 'view', '-b', '-']
    sort_cmd = ['samtools', 'sort', '-']
    aln  = subprocess.Popen(mm2_cmd, stdout=subprocess.PIPE, stderr=FNULL)
    view = subprocess.Popen(view_cmd, stdin=aln.stdout, stdout=subprocess.PIPE, stderr=FNULL)
    sort = subprocess.Popen(sort_cmd, stdin=view.stdout, stdout=subprocess.PIPE, stderr=FNULL)

    with open(outbam, 'wb') as out:
        for line in sort.stdout:
            out.write(line)

    subprocess.call(['samtools', 'index', outbam])

    return outbam


def minimap2_breakpoint(ref, cons, side):
    assert side in ('L', 'R')

    tmp_ref = str(uuid4()) + '.fa'
    with open(tmp_ref, 'w') as out:
        out.write('>ref\n%s\n' % ref.upper())

    tmp_cons = str(uuid4()) + '.fa'
    with open(tmp_cons, 'w') as out:
        out.write('>cons\n%s\n' % cons.upper())

    FNULL = open(os.devnull, 'w')

    mm2_cmd  = ['minimap2', '-x', 'map-ont', '-Y', '-a', tmp_ref, tmp_cons]
    view_cmd = ['samtools', 'view', '-b', '-']
    aln  = subprocess.Popen(mm2_cmd, stdout=subprocess.PIPE, stderr=FNULL)
    view = subprocess.Popen(view_cmd, stdin=aln.stdout, stdout=subprocess.PIPE, stderr=FNULL)

    bamstream = pysam.AlignmentFile(view.stdout, 'rb')

    for read in bamstream:
        if not read.is_unmapped and not read.is_secondary and not read.is_supplementary:
            if side == 'L':
                #print(read.get_aligned_pairs(matches_only=True)[-1])
                os.remove(tmp_ref)
                os.remove(tmp_cons)
                return read.get_aligned_pairs(matches_only=True)[-1]

            if side == 'R':
                #print(read.get_aligned_pairs(matches_only=True)[0])
                os.remove(tmp_ref)
                os.remove(tmp_cons)
                return read.get_aligned_pairs(matches_only=True)[0]

    os.remove(tmp_ref)
    os.remove(tmp_cons)

    return None


def minimap2_extend(ref, cons):

    tmp_ref = str(uuid4()) + '.fa'
    with open(tmp_ref, 'w') as out:
        out.write('>ref\n%s\n' % ref.upper())

    tmp_cons = str(uuid4()) + '.fa'
    with open(tmp_cons, 'w') as out:
        out.write('>cons\n%s\n' % cons.upper())

    FNULL = open(os.devnull, 'w')

    mm2_cmd  = ['minimap2', '-x', 'map-ont', '-Y', '-a', tmp_ref, tmp_cons]
    view_cmd = ['samtools', 'view', '-b', '-']
    aln  = subprocess.Popen(mm2_cmd, stdout=subprocess.PIPE, stderr=FNULL)
    view = subprocess.Popen(view_cmd, stdin=aln.stdout, stdout=subprocess.PIPE, stderr=FNULL)

    bamstream = pysam.AlignmentFile(view.stdout, 'rb')

    map_ref_start  = None
    map_ref_end    = None
    map_cons_start = None
    map_cons_end   = None

    for read in bamstream:
        if not read.is_unmapped:
            for q, r in read.get_aligned_pairs(matches_only=True):
                if map_ref_start is None:
                    map_ref_start = r
                    map_cons_start = q

                elif map_ref_start > r:
                    map_ref_start = r
                    map_cons_start = q

                if map_ref_end is None:
                    map_ref_end = r
                    map_cons_end = q

                elif map_ref_end < r:
                    map_ref_end = r
                    map_cons_end = q


    os.remove(tmp_ref)
    os.remove(tmp_cons)

    if None in (map_ref_start, map_ref_end, map_cons_start, map_cons_end):
        return cons

    extended_seq = ref[:map_ref_start] + cons[map_cons_start:map_cons_end] + ref[map_ref_end:]

    if len(extended_seq) <= len(cons):
        return cons

    return extended_seq


def minimap2_profile(ref, cons, max_de = 0.12):
    tmp_ref = str(uuid4()) + '.fa'
    with open(tmp_ref, 'w') as out:
        out.write('>ref\n%s\n' % ref.upper())

    tmp_cons = str(uuid4()) + '.fa'
    with open(tmp_cons, 'w') as out:
        out.write('>cons\n%s\n' % cons)

    FNULL = open(os.devnull, 'w')

    mm2_cmd  = ['minimap2', '-x', 'map-ont', '-Y', '-a', tmp_ref, tmp_cons]
    view_cmd = ['samtools', 'view', '-b', '-']
    aln  = subprocess.Popen(mm2_cmd, stdout=subprocess.PIPE, stderr=FNULL)
    view = subprocess.Popen(view_cmd, stdin=aln.stdout, stdout=subprocess.PIPE, stderr=FNULL)

    bamstream = pysam.AlignmentFile(view.stdout, 'rb')
    profiles = []

    for read in bamstream:
        if not read.is_unmapped:
            if read.has_tag('de'): # gap-compressed per-base divergence
                if read.get_tag('de') > max_de:
                    continue

            #print(read)
            profile = np.zeros(len(cons))-1

            b = 0
            for cig in read.cigartuples:
                if cig[0] in (0,1,4): # matches, insertions, and soft clips add bases to reads
                    profile[b:b+cig[1]] = cig[0]
                    b += cig[1]

            profiles.append(profile)

    if len(profiles) == 0:
        os.remove(tmp_ref)
        os.remove(tmp_cons)
        return None

    merged_profile = []

    for i in range(len(profiles[0])):
        matched = 0
        for p in profiles:
            assert int(p[i]) in (0,1,4)

            if p[i] == 0: # match
                matched = 1

        merged_profile.append(matched)

    #print(','.join(map(str, merged_profile)))

    os.remove(tmp_ref)
    os.remove(tmp_cons)

    return merged_profile


def parse_pileup(pile):
    pile = pile.upper()
    parse = []

    i = 0
    while i < len(pile):
        if pile[i] in (',','.','A','T','C','G','*'):
            parse.append(pile[i])
            i += 1

        elif pile[i] in ('+','-'):
            op = [pile[i]]
            i += 1

            num_bases = [pile[i]]
            i += 1

            while pile[i] in map(str, list(range(10))):
                num_bases.append(pile[i])
                i += 1

            num_bases = int(''.join(num_bases))

            for _ in range(num_bases):
                op.append(pile[i])
                i += 1

            parse.append(''.join(op))

        elif pile[i] == '^':
            i += 2 # char after "^" is mapq

        else:
            i += 1

    return parse


def mixmodel(pos, reg_covar=1e-4, max_iter=1000, max_pos=2):
    pos = np.asarray(pos).reshape(-1,1)

    N = np.arange(1,max_pos+1)
    models = [None for i in range(len(N))]

    for i in range(len(N)):
        models[i] = GaussianMixture(N[i], reg_covar=reg_covar, max_iter=max_iter, n_init=max_pos, random_state=1).fit(pos)

    AIC = [m.aic(pos) for m in models]

    return models[np.argmin(AIC)] # best-fit mixture


def alter_ref(ins_ref, pileup, depth):#, filter_start, filter_end):
    altered = []

    # mpileup coords are 1-based
    for i, b in enumerate(ins_ref, 1):

        if i in pileup:
            snv_p = [p for p in pileup[i] if p[0] not in ('+', '-', '.', ',')]
            ins_p = [p for p in pileup[i] if p[0] == '+']

            alt_count = 0
            alt_base = ''

            if snv_p:
                alt_base, alt_count = Counter(snv_p).most_common()[0]

            ins_count = 0
            ins_seq = ''

            if ins_p:
                ins_seq, ins_count = Counter(ins_p).most_common()[0]

            if alt_count > depth[i]*.5 and alt_count >= 3:
                #print('%d: ref %s, alt %s, pileup: %s' % (i, b, alt_base, str(pileup[i])))
                b = alt_base

            if ins_count > depth[i]*.5 and ins_count >= 3:
                #print('%d: ref %s, ins %s, pileup: %s' % (i, b, ins_seq, str(pileup[i])))
                b += ins_seq.lstrip('+-')

            if b != '*': # deletion
                altered.append(b)

    return ''.join(altered)


def process_cluster(cluster, inslib, outbase, args):
    minreads = int(args.minreads)
    mafft_threads = int(args.mafft_threads)
    te_lib = args.elts
    ref = pysam.Fastafile(args.ref)

    if cluster.prelim_filter(minreads=minreads):
        logger.debug('processing cluster %s, length: %d' % (cluster.uuid, len(cluster)))
        cluster.trim_reads(int(args.flanksize))

        for insread in cluster.reads:
            te = te_align(te_lib, insread.r_seq_trimmed)
            if te:
                insread.te_info(te)

        if cluster.te_align_count() >= minreads:
            major_family = cluster.te_major_family()
            major_subfam = cluster.te_major_subfam()
            family_frac  = cluster.te_fam_frac(major_family)

            use_reads = []

            for r_name, z in cluster.te_overlap_zscores().items():
                if abs(z) < 2:
                    use_reads.append(r_name)

            for read in cluster.reads:
                if read.r_name in use_reads and read.te:
                    read.useable = True

            if family_frac < 0.5 or None in (major_family, major_subfam):
                return None

            te_type = major_family+':'+major_subfam

            # set initial breakpoints from GMM
            r_pos_list = sorted([read.r_pos for read in cluster.reads if read.useable])
            model = mixmodel(r_pos_list)

            means = list(model.means_.reshape(1,-1)[0])
            if len(means) == 1:
                means.append(means[0])

            cluster.breakpoints = sorted(map(round, means))

            dump_fa = cluster.dump_fasta(te_fam=major_family)
            dump_fq = cluster.dump_fastq(te_fam=major_family)
            msa_fa = mafft(dump_fa, threads=mafft_threads)
            msa = MSA(msa_fa)

            cluster.cons = msa.consensus(gapped=False)
            cluster.cons = minimap2_pileup(cluster, dump_fq, outbase)

            per_bam_fq = {}

            if te_type in inslib:
                samples = list(set([read.bamname for read in cluster.reads if read.useable]))

                for sample in samples:
                    per_bam_fq[sample] = cluster.dump_fastq(te_fam=major_family, sample=sample)

                ref_region = ref.fetch(cluster.chrom(), int(cluster.breakpoints[0]) - int(args.flanksize), int(cluster.breakpoints[1]) + int(args.flanksize))
                
                cluster.cons_align_map = minimap2_profile(ref_region, cluster.cons)
                cluster.mask_cons()

                cons_seg_coords = cluster.sorted_unmapped_segments()[0]
                cons_seg_seq = cluster.cons[cons_seg_coords[0]:cons_seg_coords[1]]

                # this alignment can overrule previous subfamily classification
                cluster.cons_te_align = te_align(te_lib, cons_seg_seq)

                if cluster.cons_te_align is None:
                    cluster.reset_major_family, cluster.reset_major_subfam = 'NA', 'NA'

                else:
                    cluster.reset_major_family, cluster.reset_major_subfam = cluster.cons_te_align.subfamily().split(':')

                # handle case where te alignment is split across > 1 unaligned segment
                if cluster.cons_te_align is not None and len(cluster.sorted_unmapped_segments()) > 1:
                    joined_cons = cluster.join_unmapped_segments(int(args.flanksize), inslib[cluster.cons_te_align.subfamily()])
                    if joined_cons is not None:
                        # repeat te alignment, subfam reassignment if unmapped segments have been joined
                        cons_seg_coords = cluster.sorted_unmapped_segments()[0]
                        cons_seg_seq = cluster.cons[cons_seg_coords[0]:cons_seg_coords[1]]

                        cluster.cons_te_align = te_align(te_lib, cons_seg_seq)

                        if cluster.cons_te_align is None:
                            cluster.reset_major_family, cluster.reset_major_subfam = 'NA', 'NA'

                        else:
                            cluster.reset_major_family, cluster.reset_major_subfam = cluster.cons_te_align.subfamily().split(':')

                # clean up breakpoints
                bp_left, bp_right = cluster.sorted_unmapped_segments()[0]

                bp_left_seq = cluster.cons[bp_left-150:bp_left+50]
                bp_right_seq = cluster.cons[bp_right-50:bp_right+150]

                break_left = minimap2_breakpoint(ref_region, bp_left_seq, 'L')
                break_right = minimap2_breakpoint(ref_region, bp_right_seq, 'R')

                if None not in (break_left, break_right):
                    # fine-tune breakpoints in consensus seq
                    left_corr = 150-break_left[0]
                    right_corr = 50-break_right[0]

                    adjusted = cluster.adjust_cons(left_corr, right_corr)

                    if adjusted:
                        # repeat te alignment, subfam reassignment if breakpoints have been tweaked
                        cons_seg_coords = cluster.sorted_unmapped_segments()[0]
                        cons_seg_seq = cluster.cons[cons_seg_coords[0]:cons_seg_coords[1]]

                        cluster.cons_te_align = te_align(te_lib, cons_seg_seq)

                        if cluster.cons_te_align is None:
                            cluster.reset_major_family, cluster.reset_major_subfam = 'NA', 'NA'

                        else:
                            cluster.reset_major_family, cluster.reset_major_subfam = cluster.cons_te_align.subfamily().split(':')


                    # index from start of reference segment
                    cluster.breakpoints = sorted([(int(cluster.breakpoints[0]) - int(args.flanksize)) + break_left[1], 
                                                  (int(cluster.breakpoints[0]) - int(args.flanksize)) + break_right[1]])

                    cluster.breakpoints_remappable = True

                    #TSD
                    if break_right[1] < break_left[1]:
                        cluster.tsd_seq = ref_region[break_right[1]:break_left[1]]

                    l_ext, r_ext = cluster.extend_tsd()

                    if l_ext + r_ext > 0:
                        # repeat te alignment, subfam reassignment if breakpoints have been tweaked
                        cons_seg_coords = cluster.sorted_unmapped_segments()[0]
                        cons_seg_seq = cluster.cons[cons_seg_coords[0]:cons_seg_coords[1]]

                        cluster.cons_te_align = te_align(te_lib, cons_seg_seq)

                        if cluster.cons_te_align is None:
                            cluster.reset_major_family, cluster.reset_major_subfam = 'NA', 'NA'

                        else:
                            cluster.reset_major_family, cluster.reset_major_subfam = cluster.cons_te_align.subfamily().split(':')


                if args.detail_output:

                    longest_read = max([len(read.r_seq) for read in cluster.reads if read.useable])

                    ext_ref = ref.fetch(cluster.chrom(), (cluster.breakpoints[0] - longest_read) - int(args.flanksize), cluster.breakpoints[1] + longest_read + int(args.flanksize))
                    ext_cons = minimap2_extend(ext_ref, cluster.cons)

                    cons_out_fa = outbase + '/' + cluster.uuid + '.cons.ref.fa'

                    with open(cons_out_fa, 'w') as out:
                        out.write('>%s\n%s\n' % (cluster.uuid, ext_cons))

                    for sample in samples:
                        te_outbam = minimap2_bam(cons_out_fa, per_bam_fq[sample], '%s/%s.%s.te.bam' % (outbase, sample, cluster.uuid))

            else:
                pass
                # TODO: add "warnings" to cluster

            os.remove(msa_fa)
            os.remove(dump_fa)
            os.remove(dump_fq)

            for fq in per_bam_fq.values():
                os.remove(fq)

    return cluster


def build_clusters(args, outbase, chrom):
    min_te_len = int(args.min_te_len)
    max_te_len = int(args.max_te_len)
    wiggle = int(args.wiggle)

    bams = [pysam.AlignmentFile(bam) for bam in args.bams.split(',')]

    ins_forest = dd(Intersecter)
    clusters = dict()

    for bam in bams:

        for read in bam.fetch(contig=chrom):
            if read.is_secondary or read.is_supplementary:
                continue

            b = 0

            for cig in read.cigartuples:
                if cig[1] > min_te_len and cig[0] in (1,4): # insertion or soft clip in read
                    is_ins  = cig[0] == 1
                    is_clip = cig[0] == 4

                    if is_clip and cig[1] > max_te_len:
                        continue

                    clip_end = None
                    if is_clip:
                        clip_end = 'R'


                    q_start = b-1
                    if q_start < 0:
                        q_start = 0
                        clip_end = 'L'

                    q_end = b+cig[1]

                    ap = dict(read.get_aligned_pairs())

                    r_pos = None

                    if q_start in ap:
                        r_pos = ap[q_start]

                    if r_pos is None and q_end in ap:
                        r_pos = ap[q_end]

                    if r_pos is None:
                        logger.debug('skip read warning: %s' % read.qname)
                        continue

                    r_start = None
                    r_end = None

                    if q_start in ap:
                        r_start = ap[q_start]

                    if q_end in ap:
                        r_end = ap[q_end]

                    #ins_qual = np.asarray(list(map(ord, read.qual[q_start:q_end+1]))) - 34.
                    #flank_qual = np.asarray(list(map(ord, read.qual[:q_start])) + list(map(ord, read.qual[q_end+1:]))) - 34.

                    ins_read = InsRead(bam.filename.decode(), read.reference_name, q_start, q_end, r_start, r_end, read.qname, read.seq, read.qual, is_ins, is_clip, clip_end)

                    cluster_id = None

                    if read.reference_name in ins_forest:
                        for i, cluster in enumerate(ins_forest[read.reference_name].find(r_pos-1,r_pos+1)):
                            cluster_id = cluster.value
                            clusters[cluster_id].add(ins_read)

                            if i > 0:
                                logger.debug('r_pos overlaps multiple clusters on read: %s' % read.qname)

                    if cluster_id is None:
                        cluster_id = str(uuid4())
                        ins_forest[read.reference_name].add_interval(Interval(r_pos-wiggle, r_pos+wiggle, value=cluster_id))
                        clusters[cluster_id] = InsCluster(cluster_id)
                        clusters[cluster_id].add(ins_read)


                if cig[0] in (0,1,4): # matches, insertions, and soft clips add bases to reads
                    b += cig[1] 

    pickle_fn = '%s/%s.pickle' % (outbase, chrom)

    if len(clusters) > 0:
        logger.info('writing clusters to %s' % pickle_fn)

        with open(pickle_fn, 'wb') as p_out:
            pickle.dump(clusters, p_out)

        return pickle_fn

    return None


def assign_filter(output):
    filters = []

    if output['UnmapCover'] == 'NA':
        filters.append('UnmapCoverNA')

    elif output['UnmapCover'] < 0.5:
        filters.append('UnmapCover<0.5')

    if output['Family'] == 'NA':
        filters.append('NoFamily')

    if output['StartTE'] == 'NA':
        filters.append('NoTEAlignment')

    if not output['Remappable']:
        filters.append('NonRemappable')

    if 'NA' not in (output['StartTE'], output['EndTE']):
        if (output['EndTE'] - output['StartTE'])*2 < output['LengthIns']:
            filters.append('TEMapTooLong')

    if len(filters) > 0:
        return ','.join(filters)

    
    return 'PASS'


header = [
'UUID',
'Chrom',
'Start',
'End',
'Strand',
'Family',
'Subfamily',
'StartTE',
'EndTE',
'LengthIns',
'Inversion',
'UnmapCover',
'UsedReads',
'SpanReads',
'NumSamples',
'SampleReads',
'NonRef',
'TSD',
'Consensus',
'Remappable',
'Filter'
]

detail_header = [
'Cluster',
'BamName',
'ReadName',
'IsSpanRead',
'TEAlign',
'TEOverlap',
'Useable',
'RefPos'
]
                    

def main(args):

    if args.debug:
        logger.setLevel(logging.DEBUG)

    nonref = None
    if args.nonref:
        nonref = pysam.Tabixfile(args.nonref)

    bams = args.bams.split(',')

    outbase = '_'.join(['.'.join(os.path.basename(bam).split('.')[:-1]) for bam in bams])

    if args.outbase is not None:
        outbase = args.outbase

    logger.info('te-ont started with command: %s' % ' '.join(sys.argv))
    logger.info('output basename: %s' % outbase)

    if not os.path.exists(outbase):
        os.mkdir(outbase)

    bam = pysam.AlignmentFile(bams[0])
    chroms = bam.references

    if args.chroms is not None:
        chroms = []
        with open(args.chroms) as _:
            for line in _:
                chroms.append(line.strip())


    inslib = load_falib(args.elts)


    # parallelise cluster building across chromosomes
    pool = mp.Pool(processes=int(args.procs))

    reslist = []
    for chrom in chroms:
        res = pool.apply_async(build_clusters, [args, outbase, chrom])
        reslist.append(res)

    clusters = []

    pickles = []

    for res in reslist:
        pickle_fn = res.get()
        if pickle_fn is not None:
            pickles.append(pickle_fn)

    for pickle_fn in pickles:
        with open(pickle_fn, 'rb') as p:
            p_clusters = pickle.load(p)
            
            if len(p_clusters) > 0:
                logger.info('loaded %d clusters from %s' % (len(p_clusters), pickle_fn))

                for cluster in p_clusters.values():
                    clusters.append(cluster)

        if not args.keep_pickles:
            os.remove(pickle_fn)

    # parallelise cluster processing per cluster
    reslist = []
    for cluster in clusters:
        res = pool.apply_async(process_cluster, [cluster, inslib, outbase, args])
        reslist.append(res)

    processed_clusters = []
    for res in reslist:
        processed_clusters.append(res.get()) 

    table_out = open(outbase+'.table.txt', 'w')

    table_out.write('%s\n' % '\t'.join(header))

    c_count = 0

    for cluster in processed_clusters:
        if cluster and cluster.cons:

            nr = 'NA'

            if nonref:
                ins_start, ins_end = cluster.interval()

                if cluster.chrom() in nonref.contigs:
                    for nr in nonref.fetch(cluster.chrom(), ins_start, ins_end):
                        nr = '_'.join(nr.strip().split())

            c_unmap = cluster.sorted_unmapped_segments()[0]

            start_te    = 'NA'
            end_te      = 'NA'
            start_unmap = 'NA'
            end_unmap   = 'NA'
            unmap_cover = 'NA'

            if cluster.cons_te_align is not None:
                start_te = cluster.cons_te_align.min_qry_coord()
                end_te = cluster.cons_te_align.max_qry_coord()
                start_unmap = cluster.cons_te_align.min_tgt_coord()
                end_unmap = cluster.cons_te_align.max_tgt_coord()
                unmap_cover = (end_unmap-start_unmap)/(c_unmap[1] - c_unmap[0])

            output = od()

            output['UUID']        = cluster.uuid
            output['Chrom']       = cluster.chrom()
            output['Start']       = int(cluster.breakpoints[0])
            output['End']         = int(cluster.breakpoints[1])
            output['Strand']      = cluster.te_orientation()
            output['Family']      = cluster.te_major_family()
            output['Subfamily']   = cluster.te_major_subfam()
            output['StartTE']     = start_te
            output['EndTE']       = end_te
            output['LengthIns']   = c_unmap[1] - c_unmap[0]
            output['Inversion']   = cluster.detect_inversion()
            output['UnmapCover']  = unmap_cover
            output['UsedReads']   = cluster.te_useable_count()
            output['SpanReads']   = cluster.te_embedded_count()
            output['NumSamples']  = cluster.te_samples()
            output['SampleReads'] = cluster.te_sample_count()
            output['NonRef']      = nr
            output['TSD']         = cluster.tsd_seq
            output['Consensus']   = cluster.cons
            output['Remappable']  = cluster.breakpoints_remappable
            output['Filter']      = assign_filter(output)

            detail = od()

            if args.detail_output:
                with open(outbase+'/'+cluster.uuid+'.detail.out', 'w') as detail_out:
                    detail_out.write('%s\n' % '\t'.join(detail_header))

                    for insread in cluster.reads:
                        detail['Cluster']    = cluster.uuid
                        detail['BamName']    = insread.bamname
                        detail['ReadName']   = insread.r_name
                        detail['IsSpanRead'] = insread.is_ins
                        detail['TEAlign'] = 'NA'

                        if insread.te:
                            detail['TEAlign'] = '|'.join([','.join(a.original) for a in insread.te.alignments])

                        detail['TEOverlap'] = insread.te_overlap_frac()
                        detail['Useable'] = insread.useable
                        detail['RefPos'] = insread.r_pos

                        detail_out.write('%s\n' % '\t'.join(map(str, [detail[col] for col in detail_header])))

            if args.color_consensus:
                C_YELLOW  = "\033[33m"
                C_BLUE    = "\033[34m"
                C_PINK    = "\033[91m"
                C_END     = "\033[0m"

                col_seg = cluster.cons[c_unmap[0]:c_unmap[1]]

                if cluster.cons_te_align is not None:
                    col_seg=col_seg[:start_unmap] + C_BLUE + col_seg[start_unmap:end_unmap] + C_YELLOW + col_seg[end_unmap:]

                tsd_len = 0
                if output['TSD'] != 'NA':
                    tsd_len = len(output['TSD'])

                tsd_1 = C_PINK + cluster.cons[c_unmap[0]-tsd_len:c_unmap[0]] + C_END
                tsd_2 = C_PINK + cluster.cons[c_unmap[1]:c_unmap[1]+tsd_len] + C_END

                output['Consensus'] = cluster.cons[:c_unmap[0]-tsd_len]+tsd_1+C_YELLOW+col_seg+C_END+tsd_2+cluster.cons[c_unmap[1]+tsd_len:]

            table_out.write('%s\n' % '\t'.join(map(str, [output[col] for col in header])))
            c_count += 1

    table_out.close()

    logger.info('finished. wrote %d records to %s' % (c_count, outbase+'.table.txt'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Find TE insertions in ONT data')
    parser.add_argument('-b', '--bams', required=True, help='comma seperated bam list')
    parser.add_argument('-e', '--elts', help='reference elements .fa', required=True)
    parser.add_argument('-r', '--ref', help='reference fasta (samtools faidx indexed)', required=True)
    parser.add_argument('-p', '--procs', help='process count', default=1)
    parser.add_argument('-m', '--minreads', default=3, help='minimum read count to form cluster (default = 3)')
    parser.add_argument('-o', '--outbase', default=None, help='base name for output files (defaults to bam name(s)')
    parser.add_argument('-c', '--chroms', default=None, help='limit to chromsomes in file')
    parser.add_argument('--max_te_len', default=7000, help='maximum insertion size (default = 7000)')
    parser.add_argument('--min_te_len', default=200, help='minimum insertion size (default = 200)')
    parser.add_argument('--wiggle', default=50, help='interval padding for breakpoint search (default = 50)')
    parser.add_argument('--flanksize', default=500, help='flank size for read trimming (default = 500)')
    parser.add_argument('--mafft_threads', default=1)
    parser.add_argument('-n', '--nonref', help='known nonref tabix', default=None)
    parser.add_argument('--debug', action='store_true', default=False)
    parser.add_argument('--color_consensus', action='store_true', default=False)
    parser.add_argument('--detail_output', action='store_true', default=False, help='outputs per-insertion .bams and per-read mapping information')
    parser.add_argument('--keep_pickles', action='store_true', default=False, help='do not remove per-chromosome pickle files')
    args = parser.parse_args()
    main(args)
