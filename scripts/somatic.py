#!/usr/bin/env python

import os
import pysam
import argparse
import subprocess
from uuid import uuid4

from collections import defaultdict as dd


def strip_seq(seq):
     nt = ['A', 'T', 'C', 'G']
     return ''.join([b for b in list(seq.upper()) if b in nt])


def remap_full_check(lendict, seqdict, ref, flanksize=500):
    success = {}

    for id, seq in seqdict.items():
        seqdict[id] = strip_seq(seq)

    tmp_cons = str(uuid4()) + '.fa'

    with open(tmp_cons, 'w') as out:
        for id, seq in seqdict.items():
            out.write(f'>{id}\n{seq}\n')
            success[id] = True

    FNULL = open(os.devnull, 'w')

    mm2_cmd  = ['minimap2', '-x', 'map-ont', '-Y', '-a', ref, tmp_cons]
    view_cmd = ['samtools', 'view', '-b', '-']

    aln  = subprocess.Popen(mm2_cmd, stdout=subprocess.PIPE, stderr=FNULL)
    view = subprocess.Popen(view_cmd, stdin=aln.stdout, stdout=subprocess.PIPE, stderr=FNULL)

    save = pysam.set_verbosity(0)
    bamstream = pysam.AlignmentFile(view.stdout, 'rb')
    pysam.set_verbosity(save)

    n_rp_list = dd(list)

    for read in bamstream:
        #print(read.to_string())
        #if not read.is_unmapped and not read.is_secondary and not read.is_supplementary:
        if not read.is_unmapped:
            inslen = int(lendict[read.query_name])
            n_rp = len(read.get_reference_positions())
            exp_len = len(read.query_sequence)-inslen*.5

            #print (f'{read.query_name}\t{n_rp}\t{len(read.query_sequence)}\t{exp_len}\t{read.cigartuples}')

            if read.cigartuples[0][0]==4 and read.cigartuples[0][1] > flanksize/2:
                continue

            if read.cigartuples[-1][0]==4 and read.cigartuples[-1][1] > flanksize/2:
                continue

            n_rp_list[read.query_name].append(n_rp)

            if n_rp > exp_len:
                success[read.query_name] = False

    for uuid in n_rp_list:
        s0 = sorted(n_rp_list[uuid], reverse=True)[0]
        s1 = 0

        if len(n_rp_list[uuid]) > 1:
            s1 = sorted(n_rp_list[uuid], reverse=True)[1]
        
        #print(f'{uuid}\t{sorted(n_rp_list[uuid], reverse=True)}\t{s0}\t{s1}')

        if s1*1.1 > s0:
            success[uuid] = False

    os.remove(tmp_cons)

    return success


def predict_somatic(nonref, num_samples, filled_phase, empty_phase, max_vaf=0.4):
    if int(num_samples) > 1:
        return 'GERMLINE'
    
    if nonref != 'NA':
        return 'GERMLINE'

    if 'NA' in (filled_phase, empty_phase):
        return 'GERMLINE'

    fph = dd(list)
    eph = dd(list)

    fph_info = filled_phase.split("|")[1:]
    eph_info = empty_phase.split("|")[1:]

    shared_ps = [] # phase sets in filled and empty

    for ph in fph_info:
        ps, hp, rc = ph.split(':')
        if int(rc) >= 1:
            fph[ps].append((hp, int(rc)))
    
    for ph in eph_info:
        if ph.startswith('unphased'):
            continue

        ps, hp, rc = ph.split(':')
        if int(rc) >= 2:
            eph[ps].append((hp, int(rc)))

        if ps in fph:
            shared_ps.append(ps)
    
    shared_ps = list(set(shared_ps))

    if len(shared_ps) == 0:
        return 'GERMLINE'

    for ps in shared_ps:
        if not (len(fph[ps]) == 1 and len(eph[ps]) > 1):
            return 'GERMLINE'

    for ps in shared_ps:
        f_count = sum([ph[1] for ph in fph[ps]])
        e_count = sum([ph[1] for ph in eph[ps]])

        if f_count/e_count > max_vaf:
            return 'GERMLINE'

    return 'SOMATIC'


def main(args):

    header = []
    recs = []
    seqdict = {}
    lendict = {}

    with open(args.table, 'r') as table:
        for i, line in enumerate(table):
            if i == 0:
                header = line.strip().split('\t')
                print(line.strip()+"\tSomatic")

            else:
                rec = {}

                for n, field in enumerate(line.strip().split('\t')):
                    rec[header[n]] = field

                recs.append(rec)
                seqdict[rec['UUID']] = rec['Consensus']
                lendict[rec['UUID']] = rec['LengthIns']
    
    pass_remap = remap_full_check(lendict, seqdict, args.ref)

    for rec in recs:
        outcome = predict_somatic(rec['NonRef'], rec['NumSamples'], rec['Phasing'], rec['EmptyPhase'])
        
        if not pass_remap[rec['UUID']]:
            if rec['Filter'] == 'PASS':
                rec['Filter'] = 'AltIns'
            else:
                rec['Filter'] = rec['Filter'] + ',AltIns'
        
        print('\t'.join(rec.values()) + f'\t{outcome}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='foo')
    parser.add_argument('-t', '--table', required=True, help='table')
    parser.add_argument('-r', '--ref', required=True)
    parser.add_argument('-f', '--flanksize', default=500, help='flanksize parameter from tldr (default = 500)')
    args = parser.parse_args()
    main(args)

