#!/bin/bash

usage() { echo "## ERROR: Usage: $0 [--genome <hg19|hg38>] [--out <string>] [--minOntarget <0-Inf>] [--maxOntarget <minOntarget-Inf>] [--minOfftarget <0-Inf>] [--maxOfftarget <minOfftarget-Inf>] [--paddingOfftarget <0-Inf>]. Exit." 1>&2; exit 1; }
nogenome() { echo "## ERROR: You need to provide the genome version through --genome option (hg19 or hg38). Exit." 1>&2; exit 1; }
nooutput() { echo "## ERROR: You need to provide an output directory through --out option. Exit." 1>&2; exit 1; }
notargets() { echo "## ERROR: You need to provide a bed file with targets with --targets option. Exit." 1>&2; exit 1; }
noref() { echo "## ERROR: You need to provide a reference genome FASTA file with --ref option. Exit." 1>&2; exit 1; }

while getopts ":-:" o; do
    case "${o}" in
    	-)  
        	case $OPTARG in
                genome)
                    genome="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    [ "$genome" == "hg19" ] || [ "$genome" == "hg38" ] || usage
                    ;;
                out)
                    out="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ;;
                targets)
                    targets="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ;;
                minOntarget)
                    minOntarget=$(echo "${!OPTIND}" | bc); OPTIND=$(( $OPTIND + 1 ))
                    (($(echo "$minOntarget > 0" | bc -l))) || usage
                    ;;
                maxOntarget)
                    maxOntarget=$(echo "${!OPTIND}" | bc); OPTIND=$(( $OPTIND + 1 ))
                    (($(echo "$maxOntarget >= $minOntarget" | bc -l))) || usage
                    ;;
                minOfftarget)
                    minOfftarget=$(echo "${!OPTIND}" | bc); OPTIND=$(( $OPTIND + 1 ))
                    (($(echo "$minOfftarget > 0" | bc -l))) || usage
                    ;;
                maxOfftarget)
                    maxOfftarget=$(echo "${!OPTIND}" | bc); OPTIND=$(( $OPTIND + 1 ))
                    (($(echo "$maxOfftarget > $minOfftarget" | bc -l))) || usage
                    ;;
                paddingOfftarget)
                    paddingOfftarget=$(echo "${!OPTIND}" | bc); OPTIND=$(( $OPTIND + 1 ))
                    (($(echo "$paddingOfftarget >= 0" | bc -l))) || usage
                    ;;
                ref)
                    ref="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ;;
				*)
            		usage
            	;;
         esac ;;        
        *)
            usage
            ;;
    esac
done
shift $((OPTIND-1))

if [ -z "${genome}" ]; then
    nogenome
fi
if [ -z "${out}" ]; then
    nooutput
fi
if [ -z "${targets}" ]; then
    notargets
fi
if [ -z "${ref}" ]; then
    noref
fi
if [ -z "${minOntarget}" ]; then
    minOntarget=100
    echo " -> No use of --minOntarget option, value set as default: 100"
fi
if [ -z "${maxOntarget}" ]; then
    maxOntarget=300
    echo " -> No use of --maxOntarget option, value set as default: 300"
fi
if [ -z "${minOfftarget}" ]; then
    minOfftarget=1
    echo " -> No use of --minOfftarget option, value set as default: 1"
fi
if [ -z "${maxOfftarget}" ]; then
    maxOfftarget=50000
    echo " -> No use of --maxOfftarget option, value set as default: 50000"
fi
if [ -z "${paddingOfftarget}" ]; then
    paddingOfftarget=300
    echo " -> No use of --paddingOfftarget option, value set as default: 300"
fi


here="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

file=$here/data/${genome}_ncbiRefSeq
# NCBI genes (without .tsv)
# downloaded here: http://genome.ucsc.edu/cgi-bin/hgTables ncbiRefSeq "all fields from sleected table"
file2=$here/data/${genome}_ctgPos2
# Contigs (without.tsv)
# downloaded here: http://genome.ucsc.edu/cgi-bin/hgTables MappindandSequencing Map Contigs

if [ ! -f "$file.tsv" ]; then
    echo "Unzipping of the gene file."
    gunzip $file.tsv.gz
fi
if [ ! -f "$file2.tsv" ]; then
    echo "Unzipping of the contig file."
    gunzip $file2.tsv.gz
fi

echo "Step 1: extracting exons from gene file and selected the targeted ones"
# from genes to exons / parameter: size of side target 0bp
	# $10=exonStarts, $11=exonEnds -> chr1	67000041	67000051	SGIP1_NM_001308203.2_exon2
grep "NM_" $file.tsv | awk -F"\t" '{n=split($10,start,","); split($11,end,","); for(i=1; i<n; i=i+1) {ex=i; if($4=="-"){ex=(n-i)}; if($3 !~ /_/) {print $3 "\t" start[i]-00 "\t" end[i]+00 "\t" $13 "_" $2 "_exon" ex}  } }' | awk -F"\t" '{if($3>$2) print $0}' > $file.$out.exons.bed
bedtools sort -i $file.$out.exons.bed > $file.$out.exons.sort.bed
bedtools sort -i $file.$out.exons.bed > $file.$out.for-annotation.bed
# merge overlapping intervals
	# -> chr1	985282	985417	AGRN_NM_001305275.2_exon27,AGRN_NM_001364727.2_exon26,AGRN_NM_198576.4_exon27
bedtools merge -c 4 -o collapse -i $file.$out.exons.sort.bed > $file.$out.exons.sort.merge.bed
# merge targets file
bedtools sort -i $targets > $targets.$out.sort.bed
bedtools merge -i $targets.$out.sort.bed > $targets.$out.merged.bed
# script to separate exons covered by those which are not
Rscript $here/subscripts/annotate-targets.R $file.$out.exons.sort.merge.bed $targets.$out.merged.bed $file.$out.exons.sort.merge.ontar.bed $file.$out.exons.sort.merge.offtar.bed
bedtools sort -i ${file}.$out.exons.sort.merge.offtar.bed > $file.$out.exons.sort.merge.offtar.sort.bed

echo "Step 2: extend on-targets smaller than --minOntarget"
# extend if smaller than minOntarget 
awk -F"\t" -v minOntarget=$minOntarget '{s=$2; e=$3; if($3-$2<minOntarget) {s=$2-int((minOntarget-$3+$2)/2)-1; e=$3+int((minOntarget-$3+$2)/2)+1}; print $1 "\t" s "\t" e "\t" $4}' $file.$out.exons.sort.merge.ontar.bed > $file.$out.exons.sort.merge.ontar.ext.bed
# merge overlapping intervals
bedtools merge -c 4 -o collapse -i $file.$out.exons.sort.merge.ontar.ext.bed > $file.$out.exons.sort.merge.ontar.ext.merge.bed

echo "Step 3: split targets larger than --maxOntarget"
# split if big / paramater: target max size maxOntarget
awk -F"\t" -v maxOntarget=$maxOntarget '{if($3-$2>maxOntarget){n=int(($3-$2)/maxOntarget)+2; pos[1]=$2; pos[n]=$3; for (i=2; i<n; i=i+1){pos[i]=$2+(i-1)*int(($3-$2)/(n-1))}; for (i=1; i<n; i=i+1){print $1 "\t" pos[i] "\t" pos[i+1] "\t" $4 "-part" i} } else {print $0} }' $file.$out.exons.sort.merge.ontar.ext.merge.bed > $file.$out.exons.sort.merge.ontar.ext.merge.split.bed
bedtools sort -i $file.$out.exons.sort.merge.ontar.ext.merge.split.bed > $file.$out.exons.final.bed

echo "Step 4: padding off-targets with --paddingOfftarget"
# padding for genome-wide substraction: paddingOfftarget (not used for the actual on-targets)
awk -F"\t" -v paddingOfftarget=$paddingOfftarget '{print $1 "\t" $2-paddingOfftarget "\t" $3+paddingOfftarget "\t" $4}' $file.$out.exons.final.bed > $file.$out.exons.final.padded.bed

echo "Step 5: spliting off-targets larger than --maxOfftarget"
# removing alternative contigs
awk -F"\t" '{print $3 "\t" $4 "\t" $5}' $file2.tsv | grep -v -P "_|chrM" | grep -v "chrom" > $file2.$out.bed
bedtools sort -i $file2.$out.bed > $file2.$out.sort.bed
# substract padded exons from the genome
bedtools subtract -a $file2.$out.sort.bed -b $file.$out.exons.final.padded.bed > $file2.$out.sort.notar.bed
# split if bigger than maxOfftarget
awk -F"\t" -v maxOfftarget=$maxOfftarget '{if($3-$2>maxOfftarget){n=int(($3-$2)/maxOfftarget)+2; pos[1]=$2; pos[n]=$3; for (i=2; i<n; i=i+1){pos[i]=$2+(i-1)*int(($3-$2)/(n-1))}; for (i=1; i<n; i=i+1){print $1 "\t" pos[i] "\t" pos[i+1] "\t" "Off-target_" $1 ":" pos[i] "-" pos[i+1]} } else {print $1 "\t" $2 "\t" $3 "\t" "Off-target_" $1 ":" $2 "-" $3} }' $file2.$out.sort.notar.bed > $file2.$out.sort.notar.split.bed
bedtools sort -i $file2.$out.sort.notar.split.bed > $file2.$out.final.bed
#merging on- and off-targets
cat $file.$out.exons.final.bed $file2.$out.final.bed > $here/data/$out.hg19_both.bed
bedtools sort -i $here/data/$out.hg19_both.bed > $here/data/$out.hg19_both.sort.bed
# adding padding intervals as off-targets
bedtools subtract -a $file2.$out.sort.bed -b $here/data/$out.hg19_both.sort.bed | awk -F"\t" '{print $1 "\t" $2 "\t" $3 "\t" "Off-target_" $1 ":" $2 "-" $3}' > $here/data/$out.padded.bed
cat $here/data/$out.hg19_both.sort.bed $here/data/$out.padded.bed > $here/data/$out.hg19_both2.bed
bedtools sort -i $here/data/$out.hg19_both2.bed > $here/data/$out.hg19_both2.sort.bed

echo "Step 6: annotating non-coding exons"
# looking in which off-targets non-coding exons are belonging
bedtools intersect -a $here/data/$out.hg19_both2.sort.bed -b $file.$out.exons.sort.merge.offtar.sort.bed -wb | cut -f4,8 > $file.$out.exons-non-coding.targets.tsv
Rscript $here/subscripts/annotate-off-targets.R $here/data/$out.hg19_both2.sort.bed $file.$out.exons-non-coding.targets.tsv $here/data/$out.hg19_both.sort.noncoding.bed
grep "Off-target" $here/data/$out.hg19_both.sort.noncoding.bed > $here/data/$out.hg19_both.sort.noncoding.part1.bed
grep -v "Off-target" $here/data/$out.hg19_both.sort.noncoding.bed > $here/data/$out.hg19_both.sort.noncoding.part3.bed

echo "Step 7: removing off-targets smaller than --minOfftarget"
# remove too small than minOfftarget
awk -F"\t" -v minOfftarget=$minOfftarget '{if(($3-$2)>minOfftarget) print $0}' $here/data/$out.hg19_both.sort.noncoding.part1.bed > $here/data/$out.hg19_both.sort.noncoding.part1.filt.bed

echo "Step 8: writing output"
cat $here/data/$out.hg19_both.sort.noncoding.part1.filt.bed $here/data/$out.hg19_both.sort.noncoding.part3.bed > $here/data/$out.hg19_both.sort.noncoding.partall.bed
bedtools sort -i $here/data/$out.hg19_both.sort.noncoding.partall.bed > $here/data/$out.hg19_both.sort.noncoding.partall.sort.bed
grep -v -P "^chrY" $here/data/$out.hg19_both.sort.noncoding.partall.sort.bed > $here/data/$out.temp

echo "Step 9: annotating GC content"
awk -F"\t" '{a=int(($2+$3)/2); c1=a-10000; c2=a+10000; if(c1<1){c1=1}; print $1 "\t" c1 "\t" c2 "\t" $4}' $here/data/$out.temp > $here/data/$out.temp2
bedtools nuc -fi $ref -bed $here/data/$out.temp2 | tail -n+2 | cut -f6 > $here/data/$out.temp3
paste -d "\t" $here/data/$out.temp $here/data/$out.temp3 | sort -V -k1,1 -k2,2 > $here/data/$out.bed

echo "Step 10: cleaning temporary files"
# cleaning
rm -f $file.$out.exons* $file2.$out.sort* $file.$out.exons.final.bed $file2.$out.final.bed $here/data/$out.hg19_both.bed $file2.$out.bed $here/data/$out.hg19_both.sort.noncoding.bed $here/data/$out.hg19_NCBI_genes.exons* $here/data/$out.hg19_both*.sort.bed $here/data/*part* $here/data/$out.padded.bed $here/data/$out.hg19_both2.bed
rm -f $targets.$out.* $here/data/$out.temp*


