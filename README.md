[![DOI](https://zenodo.org/badge/205993555.svg)](https://zenodo.org/badge/latestdoi/205993555)

## tldr

*Transposons from Long DNA Reads*

# Installation

tldr requires python > 3.6 and has the following dependencies:

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

Install tldr package + python dependencies:

```
python setup.py install
```

# Running tldr

Synopsis (minimal input requirements), assuming reads aligned to hg38 using minimap2:
```
tldr -b aligned_reads.bam -e /path/to/tldr/ref/teref.ont.human.fa -r /path/to/minimap2-indexed/reference/genome.fasta --color_consensus
```

## Command-line Options

### multiple .bam files
Multiple .bam files can provided in a comma-delimited list.

### -o/--outbase
Specify a base name for output files. The default is to use the name of the input bam(s) without the .bam extension and joined with "_" if > 1 .bam file given

### -p/--procs
Spread work over _p_ processes. Uses python multiprocessing.

### -n/--nonref
Annotate insertion with known non-reference insertion sites (examples provided in `/path/to/tldr/ref`)

### -c/--chroms
Specify a text file of chromosome names (one per line) and tldr will focus only on these.

### -m/--minreads
Minimum supporting read count to trigger a consensus / insertion call (default = 3)

### --max_te_length
Maximum insertion size (default = 7000)

### --min_te_len
Minimum insertion size (default = 200)

### --wiggle
Allows for sloppy breakpoints in initial breakpoint search (default = 50)

### --flanksize
Trim reads to contain at most `--flanksize` bases on either side of the insertion. Setting too large makes consensus building slower and more error-prone.

### --mafft_threads
Number of threads to give each mafft job (consensus building)

### --color_consensus
This will annotate the consensus sequence with ANSI escape characters that yield coloured text on the command-line:
red = TSD, blue = TE insertion sequence, yellow = non-TE insertion sequence

### --detail_output
Creates a directory (name is the output base name) with extended consensus sequences, per-insertion read mapping information and per-insertion .bam files. Required for mCpG analysis.

### --extend_consensus
If --detail_output option is enabled, extend output per-sample consensus by n bases (default 0). This is useful in the analysis of CpG methylation to add context on either end of the insertion.

## Output

Some fields in the output table (basename.table.txt) may not be self-explainatory:

### StartTE / EndTE
Start / end position relative to TE consensus provided via `-e/--elts`

### LengthIns
Length of actual inserted sequence. _Not necessarily the same as EndTE-StartTE_

### Inversion
Internal inversion detected in TE

### UnmapCover
Fraction of inserted sequence covered by TE sequence

### UsedReads
Number of reads used in consensus generation

### SpanReads
Number of reads which completely embed the insertion

### NumSamples
Number of samples (.bam files) in which the insertion was detected

### SampleReads
Per-sample accounting of supporting reads

### NonRef
If `-n/--nonref` given, annotate whether insertion is a known non-reference insertion ("NA" otherwise)

### TSD
Target site duplication (based on reference genome)

### Consensus
Upper case bases = reference genome sequence, lower case bases = insertion sequence. If `--color_consensus` given TSD will be red, TE will be blue, other inserted sequence (e.g. transduction) will be yellow using ANSI terminal colours (may be affected by specific terminal config)

### Filter
Annotate whether an insertion call is problematic; "PASS" otherwise (similar to VCF filter column).

## Methylation

Non-reference methylation can be assessed through the use of scripts located in the `scripts/` directory:

| script                | description |
|-----------------------|-------------|
| tldr_callmeth.sh      | Must be run from within the diretory where `nanopolish index` was run to index a .fastq file against a set of ONT .fast5 files. Takes as input a .fastq (indexed via `nanopolish index`), an output directory created via the `--detail_output` option, a UUID and a sample name. Creates a tabix indexed table from the output of nanopolish call-methylation on the sample+uuid combination. Can be automated via xargs or GNU parallel. |
| tablemeth_nonref.py   | Creates a table with per-element mCpG summary data given a tldr output table and the directory created by `--detail_output`. Only considers element + sample combinations from the tldr table where `tldr_callmeth.sh` has been run. Requires pysam, pandas, numpy, and scipy. |
| plotmeth_nonref.py    | Makes a plot of a TE (requires running `tldr_callmeth.sh` first) plus the surrounding region if `--extend_consensus` is specified. Tracks include translation to CpG space, raw log-likelihood, and smoothed methylation fraction. Requires pysam, pandas, numpy, scipy, matplotlib, and seaborn. |

## Reference TEs

See https://github.com/adamewing/te-nanopore-tools

## References

Adam D. Ewing, Nathan Smits, Francisco J. Sanchez-Luque, Sandra R. Richardson, Seth W. Cheetham, Geoffrey J. Faulkner. Nanopore Sequencing Enables Comprehensive Transposable Element Epigenomic Profiling. 2020. Molecular Cell, Online ahead of print: https://doi.org/10.1016/j.molcel.2020.10.024

## Getting help

Reporting [issues](https://github.com/adamewing/tldr/issues) and questions through github is preferred versus e-mail.
