[![DOI](https://zenodo.org/badge/205993555.svg)](https://zenodo.org/badge/latestdoi/205993555)

## tldr

*Transposons from Long DNA Reads*

# Installation

tldr requires python > 3.6 and has the following dependencies:
- HTSLIB/Samtools
- minimap2
- MAFFT
- Exonerate
- some python dependecies in the background. 

## One-step Conda environment setup 
There is a pre-baked Conda (or [mamba](https://anaconda.org/conda-forge/mamba)) environment file provided (tldr.yml) that can be used to create a tldr Conda environment with all of the necessary dependencies. 

```
git clone https://github.com/adamewing/tldr.git
cd tldr
conda env create -f tldr.yml
conda activate tldr
python setup.py install
tldr -h
```
If you use the above method, make sure to activate the Conda environment first with `conda activate tldr` whenever using tldr.

## Installing dependencies seperately
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

Note: different versions of MAFFT may yield different results from tldr. We currently recommend MAFFT v7.480.

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

### -b/--bams
Multiple .bam files can provided in a comma-delimited list.

### -e/--elts
Reference elements in .fasta format. The header for each should be formatted as`>Superfamily:Subfamily` e.g. `>ALU:AluYb9`.
If `none` is specifed instead of a filename, tldr will run without a reference TE collection. This is useful for genomes where active mobile element content is not well understood or for unbiased identification of inserted sequenced and is also useful for identifying viral intregration and gene retrocopy insertions.

### -r/--ref
Reference genome .fasta, expects a samtools index i.e. `samtools faidx`.

### -p/--procs
Spread work over _p_ processes. Uses python multiprocessing.

### -m/--minreads
Minimum supporting read count to trigger a consensus / insertion call (default = 3)

### --embed_minreads
Minimum number of reads completely embeddeding the insertion (default = 1, requires at least 1).

### -o/--outbase
Specify a base name for output files. The default is to use the name of the input bam(s) without the .bam extension and joined with "_" if > 1 .bam file given

### -c/--chroms
Specify a text file of chromosome names (one per line) and tldr will focus only on these.

### --max_te_len
Maximum insertion size (default = 10000)

### --min_te_len
Minimum insertion size (default = 200)

### --min_alt_frac
Parameter for allowing base changes in consensus cleanup (default = 0.5)

### --min_alt_depth
Parameter for allowing base changes in consensus cleanup (default = 3)

### --min_total_depth_frac
Parameter for allowing base changes in consensus cleanup (default = 0.25)

### --max_cluster_size
Limit cluster size and downsample clusters larger than the cutoff (default = no limit). Downsampling is biased such that reads completely embedding the inserted sequence are preferred.

### --wiggle
Allows for sloppy breakpoints in initial breakpoint search (default = 50)

### --flanksize
Trim reads to contain at most `--flanksize` bases on either side of the insertion. Setting too large makes consensus building slower and more error-prone.

### -n/--nonref
Annotate insertion with known non-reference insertion sites (examples provided in `/path/to/tldr/ref`

### --color_consensus
This will annotate the consensus sequence with ANSI escape characters that yield coloured text on the command-line:
red = TSD, blue = TE insertion sequence, yellow = non-TE insertion sequence
While this looks nice on the command line (try `less -R`) and is helpful for evaluating insertion calls, the output may not translate well to other applications as the escape sequences for the ANSI colours will be embedded in the sequence.

### --detail_output
Creates a directory (name is the output base name) with extended consensus sequences, per-insertion read mapping information and per-insertion .bam files. Required for mCpG analysis.

### --extend_consensus
If --detail_output option is enabled, extend output per-sample consensus by n bases (default 0). This is useful in the analysis of CpG methylation to add context on either end of the insertion.

### --trdcol
Adds 5' and 3' transduction columns needed by the `call_transductions.py` script, if you're into that kind of thing.

### --keep_pickles
Saves pickles for later.

### --use_pickles <folder>
Search specified folder for .pickle files and use them instead of clustering reads. Faster for re-running with different options, requires `--keep_pickles`.

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

### MedianMapQ
Median mapping quality score from input .bam(s)

### TEMatch
Overall mean identity to TE in reference library (`-e/--elts`)

### UsedReads
Number of reads used in consensus generation

### SpanReads
Number of reads which completely embed the insertion

### NumSamples
Number of samples (.bam files) in which the insertion was detected

### SampleReads
Per-sample accounting of supporting reads

### EmptyReads
Number of reads spanning both TSDs +/- `--wiggle` parameter with no evidence for insertion, useful for inferring genotype

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
