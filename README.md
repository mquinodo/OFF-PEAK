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



#### Test sample
A testing VCF with random variants can be found in the Test folder. The following command can be done to test AutoMap:
```
AUTOMAP_HOME=<path/to/automap/script>
bash $AUTOMAP_HOME/AutoMap_v1.0.sh
  --vcf $AUTOMAP_HOME/Test/TestSample.vcf
  --genome hg19
  --out $AUTOMAP_HOME/Results-test
```
And the output can be compared with the expected output files in the Test folder.
