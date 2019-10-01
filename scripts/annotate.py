#!/usr/bin/env python

import sys
import pysam
import argparse

def main(args):
    tbx = pysam.Tabixfile(args.tabix)
    
    with open(args.table) as l1seq:
        for line in l1seq:

            if line.startswith('UUID'):
                header = line.strip().split()
                header.append(args.name)
                print('\t'.join(header))
                continue

            c = line.strip().split()
            chrom = c[1]
            start = int(c[2]) - 100
            end   = int(c[3]) + 100

            annotations = []

            if chrom in tbx.contigs:
                for rec in tbx.fetch(chrom, start, end):
                    annotations.append('|'.join(rec.strip().split()))

            annotations = list(set(annotations)) # uniqify

            if len(annotations) == 0: annotations.append('NA')

            print(line.strip() + '\t' + ','.join(annotations))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='annotate TLDR table')
    parser.add_argument('-t', '--table', required=True)
    parser.add_argument('-x', '--tabix', required=True)
    parser.add_argument('-n', '--name', required=True)
    args = parser.parse_args()
    main(args)
