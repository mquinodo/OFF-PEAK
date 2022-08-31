# OFF-PEAK
CNV detection tool for WES data

This software was written by Mathieu Quinodoz in the group of Prof. Rivolta from the IOB in Basel, Switzerland. It was developped on Ubuntu 22.04.1 LTS.

## Prerequisites
+ BEDTools [[Link](https://bedtools.readthedocs.io/en/latest/content/installation.html)] (>= v2.25.0)
+ mosdepth [[Link](https://github.com/brentp/mosdepth)] (>= v0.3.4)
+ R [[Link](https://cran.r-project.org/mirrors.html)] (>= v3.2.0)
+ R libraries: optparse, gplots, ExomeDepth, pROC and caTools

## Installation
The tool does not require compilation.

## Usage
This tool contains four modules: target processing, coverage computation, CNV discovery and CNV plot.

### 1) Target processing
The main script 01_targets-processing.sh takes as input a BED file which contains the target regions of the sequencing data as well as the reference genome (FASTA format). The output is a BED file containing processed on-target and off-target regions.
It is called with bash and its computation time for an exome BED file is few minutes:
```
bash 01_targets-processing.sh
  --genome [hg19|hg38]
  --targets target.bed
  --out output-name
  --ref ref_genome.fa
  [other options]
```
#### Required arguments
Option | Value | Description
--- | --- | ---
--targets | STRING | BED file containing the target regions of the exome or targeted sequencing
--genome | [hg19/hg38] | Genome build used
--out | STRING | Name of the output file (without extension)
--ref | STRING | Reference genome in FASTA format

#### Optional arguments
Option | Default | Value | Description
--- | --- | --- | ---
--minOntarget | 100 | 1-Inf | Minimal size of on-targets. Smaller on-target regions will be extended on each side to reach this value.
--maxOntarget | 300 | >minOntarget | Maximal size of on-targets. Larger on-target regions will be splitted in regions of equal size.
--minOfftarget | 1000 | 1-Inf | Minimal size of off-targets. Smaller on-target regions will be discarded.
--maxOfftarget | 50000 | >minOfftarget | Maximal size of off-targets. Larger on-target regions will be splitted in regions of equal size.
--paddingOfftarget | 300 | 0-Inf | Padding around on-target regions to avoid reads in off-target regions.


### 2) Coverage computation
The main script 02_coverage-count.sh takes as input a text file listing the sample BAM files and IDs as well as the processed BED file from step 1. It is outputing a file containing the coverage for each sample and for each target.
It is called with bash and its computation time is 2-3 minutes per BAM file:
```
bash 02_coverage-count.sh
  --listBAM list_BAM.txt
  --mosdepth mosdepth-executable
  --work output_directory
  --targetsBED processed-targets.bed
  [other options]
```
#### Required arguments
Option | Value | Description
--- | --- | ---
--listBAM | STRING | Text file containing two tab-delimited columns: BAM files and IDs. It can also have only BAM files and IDs will be deduced from the name of BAM files (not recommended).
--mosdepth | STRING | Mosdepth executable file that can be downloaded on the github page [[Link](https://github.com/brentp/mosdepth)]
--work | STRING | Output directory
--targetsBED | STRING | BED file with processed targets produced in step 1.



