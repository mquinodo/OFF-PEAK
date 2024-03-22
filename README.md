# OFF-PEAK
CNV detection tool for WES and targeted sequencing data.

This software was written by Mathieu Quinodoz in the group of Prof. Rivolta from the IOB in Basel, Switzerland. It was developped on Ubuntu 22.04.1 LTS. It is published in AJHG.

## Prerequisites
+ BEDTools [[Link](https://bedtools.readthedocs.io/en/latest/content/installation.html)] (>= v2.25.0)
+ mosdepth [[Link](https://github.com/brentp/mosdepth)] (>= v0.3.3)
+ R [[Link](https://cran.r-project.org/mirrors.html)] (>= v3.2.0) with following libraries: optparse, gplots, ExomeDepth, pROC and caTools

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
  --name target-panel
  --ref ref_genome.fa
  [other options]
```
The output file (target-panel.bed) will be placed in the data folder.

#### Required arguments
Option | Value | Description
--- | --- | ---
--targets | STRING | BED file containing the target regions of the exome or targeted sequencing
--genome | [hg19/hg38] | Genome build used
--name | STRING | Name of the output file (without extension)
--ref | STRING | Reference genome in FASTA format (same as used to map the reads in BAM files)

#### Optional arguments
Option | Default | Value | Description
--- | --- | --- | ---
--minOntarget | 100 | 1 - Inf | Minimal size of on-targets. Smaller on-target regions will be extended on each side to reach this value.
--maxOntarget | 300 | > minOntarget | Maximal size of on-targets. Larger on-target regions will be splitted into regions of equal size.
--minOfftarget | 1 | 1 - Inf | Minimal size of off-targets. Smaller on-target regions will be discarded.
--maxOfftarget | 50000 | > minOfftarget | Maximal size of off-targets. Larger on-target regions will be splitted into regions of equal size.
--paddingOfftarget | 300 | 0 - Inf | Padding around on-target regions to avoid targeted reads to be counted in off-target regions.


### 2) Coverage computation
The main script 02_coverage-count.sh takes as input a text file listing the sample BAM files (.bai index files are needed in the same folder, they can be generated with samtools index command) and IDs as well as the processed BED file from step 1. It is outputing a file containing the coverage for each sample and for each target.
It is called with bash and its computation time is 2-3 minutes per sample for WES:
```
bash 02_coverage-count.sh
  --listBAM list_BAM.txt
  --mosdepth mosdepth-executable
  --work output_directory
  --targetsBED data/target-panel.bed
```
#### Required arguments
Option | Value | Description
--- | --- | ---
--listBAM | STRING | Text file containing two tab-delimited columns: BAM files (with absolute location, for example /user/Download/Patient1.bam) and IDs (for example Patient1). It can also have only one column with BAM files. In this case, IDs will be deduced from the name of BAM files (not recommended).
--mosdepth | STRING | Mosdepth executable file that can be downloaded from the github page here: [[Link](https://github.com/brentp/mosdepth/releases/download/v0.3.3/mosdepth)]
--work | STRING | Output directory
--targetsBED | STRING | BED file with processed targets produced in step 1


### 3) CNV discovery
The main script 03_OFF-PEAK.R takes as input the coverage data per target computed in step 2 and outputs multiple files including annotated CNVs, CNV plots and statistics. 
It is called with Rscript and its computation time is ~1 hour per 10 samples:
```
Rscript 03_OFF-PEAK.R
  --output output_directory
  --data ALL.target.tsv
  --databasefile data/data-hg19.RData
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
--nb-plots | 10 | 0 - Inf | Number of plots for CNVs per individual

### 4) CNV plotting
In step 3 , the 20 best CNVs per sample will be automatically plotted. This script is used to further plot other CNVs. The main script 04_OFF-PEAK-plot.R takes as input intermediate data from step 3 in order to do additional plots.
It is called with Rscript and its computation time is few seconds per CNV:
```
Rscript 04_OFF-PEAK-plot.R
  --ID ID
  --chr chr
  --begin begin
  --end end
  --batch output_directory
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
--batch | STRING | Output directory of step 3 (output_directory)
--out | STRING | Output directory
--databasefile | STRING | Absolute path to RData file containing various information (data-hg19.RData or data-hg38.RData)

#### Optional arguments
Option | Default | Value | Description
--- | --- | --- | ---
--side | 20 | 1 - 100 | Number of target plotted on each side of the region.


## Output files

The output files are grouped in 6 general folders and one folder per sample.

#### 01_general_stats with:
File name | Explanation
--- | ---
Pairwise-correlations-all.tsv | Text file with pairwise correlation between samples
Heatmap-correlations-all.pdf | PDF file with heatmap of the pairwise correlations
Maximal-correlation-per-sample.tsv | Text file with maximal correlation per sample
Maximal-correlation-per-sample.pdf | PDF file representing maximal correlation per sample
log_parameters.tsv | Text file with the log of parameters used

#### 02_BED-files
File name | Explanation
--- | ---
CNVs-all-IGV.bed | BED file that be loaded on IGV for visualization of all CNVs
CNVs-targets-only-IGV.bed | BED file that be loaded on IGV for visualization of on-target CNVs
CNVs-all-AnnotSV.bed | BED file that be used for AnnotSV annotation of all CNVs
CNVs-targets-only-AnnotSV.bed | BED file that be used for AnnotSV annotation of on-target CNVs

#### 03_Samples-info
File name | Explanation
--- | ---
Samples_quality_info.tsv | Quality information about samples
Samples_low-quality.tsv | Samples with low quality

#### 04_CNVs-results
File name | Explanation
--- | ---
CNVs-all.tsv | All CNVs detected
CNVs-targets-only.tsv | On-target CNVs
The folder also includes filtered CNVs:
Type | Explanation
--- | ---
HQ | only high quality CNVs (PQ>2, CQ>2.5, QUAL>3)
unique | no other samples with overlapping CNVs
HQ-sample | only samples of high quality

#### 05_RData-files
File name | Explanation
--- | ---
data-ID.RData | Data for all targets
data-plot-ID.RData | Data used for plotting CNVs
data-targets-ID.RData | Data for on-targets only

#### 06_PC-plots
File name | Explanation
--- | ---
PC-removal-ID.pdf | PDF representing AUCs on fake CNVs depending on PC removal
Variance-explained-ID.pdf | PDF representing the variance explained by PCs

#### Individual folder
File name | Explanation
--- | ---
ID-merged.tsv | All CNVs found in the sample
ID-merged_targets-only.tsv | On-target CNVs found in the sample
##### Subfolders:
Sufolder name | Explanation
--- | ---
plots_on-targets_off-targets | Folder with graphical representation of best CNVs
plots_on-targets-only | Folder with graphical representation of best on-target CNVs
plots_genome | Folder with genome-wide representation of coverage
plots_chromosomes | Folder with chromosome representation of coverage

#### Explanation of fields in CNV text files:
Field | Explanation
--- | ---
ID | ID of the sample
Chromosome | Chromosome
Begin | Coordinate of the beginning of the CNV
End | Coordinate of the end of the CNV
Begin-min | Minimal begin position of detected CNV (including low quality targets)
End-max | Maximal end position of detected CNV (including low quality targets)
Type | Deletion or duplication
Ploidy | Most likely ploidy based on Z-scores
Nb_targets | Number of targets in the CNV
Targets | ID of targets in the CNV
Genes | Genes affected by the CNV
Exons | Exons affected by the CNV
Genes-possible | Genes affected by minimal to maximal positions (including low quality targets)
ncRNAs | non-coding RNAs affected by the CNV
RefSeq_Functional_Element | Function elements affected by the CNV
Nb_overlapping_samples | Number of samples with an overlapping CNV
gnomAD-CNV_1% | Frequent included CNVs from gnomAD with AF>1%
gnomAD-CNV_ALL | All included CNVs from gnomAD
Counts | Total coverage in the targets included in the CNV
Average_counts | Average total coverage for the controls
SD_counts | Standard deviation of the total coverage for the controls
Ratio | Ratio between sample coverage and control coverage
Z-score_cor | Corrected Z-score for the coverage
PQ | Ploidy quality
CQ | CNV quality (based on "shape")
QUAL | Overall quality
Rank_sample | Rank of the CNV in the sample based on absolute Z-score
Nb_CNV_HQ | Number of high quality CNVs in the sample
ClinVar-patho_included | All included CNVs from ClinVar
ClinVar-patho_overlap | All overlaping CNVs from ClinVar
Nb-exon-before-chr | Number of exons between the CNV and the beginning of the chromosome
Nb-exon-after-chr | Number of exons between the CNV and the end of the chromosome

