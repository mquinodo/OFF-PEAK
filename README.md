# OFF-PEAK
CNV detection tool for WES data

This software was written by Mathieu Quinodoz in the group of Prof. Rivolta from the IOB in Basel, Switzerland. It was developped on Ubuntu 22.04.1 LTS.

It can be used under a Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. For commercial uses of OFF-PEAK, please contact mathieu.quinodoz[at]iob.ch.

## Prerequisites
+ BEDTools [[Link](https://bedtools.readthedocs.io/en/latest/content/installation.html)] (>= v2.25.0)
+ mosdepth [[Link](https://github.com/brentp/mosdepth)] (>= v0.3.3)
+ R [[Link](https://cran.r-project.org/mirrors.html)] (>= v3.2.0)
+ R libraries: optparse, gplots, ExomeDepth, pROC and caTools

## Installation
The tool does not require compilation.

## Usage
This tool contains four modules: 1° target processing, 2° coverage computation, 3° CNV discovery and 4° CNV plot.

### 1) Target processing
The main script 01_targets-processing.sh takes as input a BED file which contains the target regions of the sequencing data as well as the reference genome (FASTA format). It is outputing a BED file containing processed on-target and off-target regions.
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
--minOntarget | 100 | 1 - Inf | Minimal size of on-targets. Smaller on-target regions will be extended on each side to reach this value.
--maxOntarget | 300 | > minOntarget | Maximal size of on-targets. Larger on-target regions will be splitted into regions of equal size.
--minOfftarget | 1 | 1 - Inf | Minimal size of off-targets. Smaller on-target regions will be discarded.
--maxOfftarget | 50000 | > minOfftarget | Maximal size of off-targets. Larger on-target regions will be splitted into regions of equal size.
--paddingOfftarget | 300 | 0 - Inf | Padding around on-target regions to avoid targeted reads to be counted in off-target regions.


### 2) Coverage computation
The main script 02_coverage-count.sh takes as input a text file listing the sample BAM files and IDs as well as the processed BED file from step 1. It is outputing a file containing the coverage for each sample and for each target.
It is called with bash and its computation time is 2-3 minutes per BAM file:
```
bash 02_coverage-count.sh
  --listBAM list_BAM.txt
  --mosdepth mosdepth-executable
  --work output_directory
  --targetsBED processed-targets.bed
```
#### Required arguments
Option | Value | Description
--- | --- | ---
--listBAM | STRING | Text file containing two tab-delimited columns: BAM files and IDs. It can also have only one column with BAM files. In this case, IDs will be deduced from the name of BAM files (not recommended).
--mosdepth | STRING | Mosdepth executable file that can be downloaded from the github page here: [[Link](https://github.com/brentp/mosdepth/releases/download/v0.3.3/mosdepth)]
--work | STRING | Output directory
--targetsBED | STRING | BED file with processed targets produced in step 1


### 3) CNV discovery
The main script 03_OFF-PEAK.R takes as input the coverage data per target computed in step 2 and outputs multiple files including annotated CNVs, CNV plots and statistics. 
It is called with Rscript and its computation time is ~1 hour per 10 samples:
```
Rscript 03_OFF-PEAK.R
  --output output_directory
  --data target-coverage.tsv
  --databasefile data-hg19.RData
  [other options]
```
#### Required arguments
Option | Value | Description
--- | --- | ---
--output | STRING | BED file containing the target regions of the exome or targeted sequencing
--data | STRING | Genome build used
--databasefile | STRING | Absolute path to RData file containing various information found in data folder (data-hg19.RData or data-hg38.RData)

#### Optional arguments
Option | Default | Value | Description
--- | --- | --- | ---
--mincor | 0.9 | 0 - 0.99 | Minimal correlation of control samples compared to analyzed one.
--minsignal | 2500 | 1 - Inf | Minimal signal for a target to be analyzed.
--maxvar | -0.2 | -1 - 1 | Maximal variance for a target to be analyzed.
--leaveoneout | 1 | 0 or 1 | If 1, leave-one-out PCA (LOO-PCA) will be used. If 0, standard PCA will be used.
--downsample | 20000 | 100 - Inf | Number of downsampled targets for optimization of PC removal.
--nbFake | 500 | 10 - 10000 | Number of fake CNVs used for optimization of PC removal.
--stopPC | 0.0001 | 0 - 0.1 | Stopping criteria for optimization of PC removal.
--minZ | 4 | 2 - 10 | Minimal absolute Z-score for single target CNV processing.
--minOfftarget | 1000 | 1 - Inf | Minimal size of off-targets without exons.
--chromosome-plots | - | - | If present, coverage plots for each chromosome will be done.
--genome-plots | - | - | If present, genome-wide coverage plots will be done.

### 4) CNV plotting
In step 3 , the 20 best CNVs per sample will be automatically plotted. This script is used to further plot other CNVs. The main script 04_OFF-PEAK-plot.R takes as input intermediate data from step 3 in order to do additional plots.
It is called with Rscript and its computation time is few seconds per CNV:
```
Rscript 04_OFF-PEAK-plot.R
  --ID ID
  --chr chr
  --begin begin
  --end end
  --data 05_RData-files/data-ID.RData
  --out output-directory
  --databasefile data-hg19.RData
  [other options]
```
#### Required arguments
Option | Value | Description
--- | --- | ---
--ID | STRING | Sample ID
--chr | STRING | Chromosome
--begin | number | Begin position of CNV
--end | number | End position of CNV
--data | STRING | Intermediate RData file outputted in step 3 (in 05_RData-files folder)
--out | STRING | Output directory
--databasefile | STRING | Absolute path to RData file containing various information (data-hg19.RData or data-hg38.RData)

#### Optional arguments
Option | Default | Value | Description
--- | --- | --- | ---
--side | 10 | 1 - 100 | Number of target plotted on each side of the region.


## Output files

The output files are grouped in 6 general folders and one folder per sample.

### 01_general_stats with:

File name | Explanation | test
Pairwise-correlations-all.tsv | Text file with pairwise correlation between samples
Heatmap-correlations-all.pdf | PDF file with heatmap of the pairwise correlations
Maximal-correlation-per-sample.tsv | Text file with maximal correlation per sample
Maximal-correlation-per-sample.pdf | PDF file representing maximal correlation per sample
log_parameters.tsv | Text file with the log of parameters used

### 02_BED-files
CNVs-all-IGV.bed
CNVs-targets-only-IGV.bed
CNVs-all-AnnotSV.bed
CNVs-targets-only-AnnotSV.bed

### 03_Samples-info
Samples_quality_info.tsv
Samples_low-quality.tsv

### 04_CNVs-results
CNVs-all.tsv
CNVs-targets-only.tsv
The folder also includes filtered CNVs:
HQ: only high quality CNVs (PQ>2, CQ>2.5, QUAL>3)
unique: no other samples with overlapping CNVs
HQ-sample: only samples of high quality

### 05_RData-files
data-ID.RData
data-plot-ID.RData
data-targets-ID.RData

### 06_PC-plots
PC-removal-ID.pdf
Variance-explained-ID.pdf

### Individual folder
ID-merged.tsv
ID-merged_targets-only.tsv
#### Subfolders:
plots_on-targets_off-targets
plots_on-targets-only
plots_genome
plots_chromosomes
