## teONT

*TE insertion detection in long-read error-prone sequencing data*

# Installation

teONT requires python > 3.6 and has the following dependencies:

## HTSLIB / SAMtools
Easiest method is via conda:

```
conda install -c bioconda tabix
conda install -c bioconda samtools
```
Manual installation:
```
git clone https://github.com/samtools/htslib.git
git clone https://github.com/samtools/samtools.git

make -C htslib && sudo make install -C htslib
make -C samtools && sudo make install -C samtools
```

## minimap2
Via conda:
```
conda install -c bioconda minimap2
```
For manual installation see [minimap2 github](https://github.com/lh3/minimap2)


## MAFFT
Via conda:
```
conda install -c bioconda mafft
```

For manual installation see the [mafft website](https://mafft.cbrc.jp/alignment/software/linux.html)


## Exonerate
Via conda:

```
conda install -c bioconda exonerate
```
For manual installation see the [exonerate website](https://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate)

# Install

Install teONT package + python dependencies:

```
python setup.py install
```

# Running teONT

Synopsis (minimal input requirements), assuming reads aligned to hg38 using minimap2:
```
teont -b aligned_reads.bam -e /path/to/teont/ref/teref.human.fa -r /path/to/minimap2-indexed/reference/genome.fasta
```

# Additional options:

## multiple .bam files
Multiple .bam files can provided in a comma-delimited list.

## -o/--outbase
Specify a base name for output files. The default is to use the name of the input bam(s) without the .bam extension and joined with "_" if > 1 .bam file given

## -p/--procs
Spread work over _p_ processes. Uses python multiprocessing.

## -n/--nonref
Annotate insertion with known non-reference insertion sites (examples provided in `/path/to/teont/ref`)

## -c/--chroms
Specify a text file of chromosome names (one per line) and teont will focus only on these.

## -m/--minreads
Minimum supporting read count to trigger a consensus / insertion call (default = 3)

## --max_te_length
Maximum insertion size (default = 7000)

## --min_te_len
Minimum insertion size (default = 200)

## --wiggle
Allows for sloppy breakpoints in initial breakpoint search (default = 50)

## --flanksize
Trim reads to contain at most `--flanksize` bases on either side of the insertion. Setting too large makes consensus building slower and more error-prone.

## --mafft_threads
Number of threads to give each mafft job (consensus building)

## --color_consensus
This will annotate the consensus sequence with ANSI escape characters that yield coloured text on the command-line:
red = TSD, blue = TE insertion sequence, yellow = non-TE insertion sequence

## --detail_out
Creates a directory (name is the output base name) with extended consensus sequences, per-insertion read mapping information and per-insertion .bam files


