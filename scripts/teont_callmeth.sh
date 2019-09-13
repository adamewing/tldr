#!/bin/sh

if [ $# -ne 4 ]
then
    echo "USAGE: $0 <NANOPOLISH INDEXED READS FASTQ> <CALLDIR> <UUID> <SAMPLE>"
    echo "NOTE: must be executed from the directory where nanopolish index was run"
    exit 65
fi

READS=$1
CALLDIR=$2
UUID=$3
SAMPLE=$4

if [ ! -e $READS ]
then
    echo "reads fastq not found: $READS"
    exit 1
fi

if [ ! -e ${READS}.index ]
then
    echo "reads fastq not indexed: $READS"
    exit 1
fi

if [ ! -d $CALLDIR ]
then
    echo "call directory (--detail_out) not found: $CALLDIR"
    exit 1
fi

CONS=${CALLDIR}/${UUID}.cons.ref.fa
BAM=${CALLDIR}/${SAMPLE}.${UUID}.te.bam

if [ ! -e $CONS ]
then
    echo "expected consensus fasta not found: $CONS"
    exit 1
fi

if [ ! -e $BAM ]
then
    echo "expected BAM not found: $BAM"
    exit 1
fi

command -v nanopolish >/dev/null 2>&1 || { echo "nanopolish is not installed" >&2; exit 1; }
command -v bgzip >/dev/null 2>&1 || { echo "bgzip is not installed" >&2; exit 1; }
command -v tabix >/dev/null 2>&1 || { echo "tabix is not installed" >&2; exit 1; }

echo "fetching nanopolish methylation calls for $UUID sample $SAMPLE..."

TSV=${CALLDIR}/${UUID}.te.meth.tsv

nanopolish call-methylation -r $READS -g $CONS -b $BAM | sort -k3,3n > $TSV

bgzip $TSV

tabix -f -S 1 -s 1 -b 3 -e 4 ${TSV}.gz

echo "finished, methylation data indexed in ${TSV}.gz"
