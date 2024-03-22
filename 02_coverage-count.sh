#!/bin/bash

nomosdepth() { echo "## ERROR: Mosdepth not found (--mosdepth). Exit." 1>&2; exit 1; }
nowork() { echo "## ERROR: You need to provide the working directory through --work option. Exit." 1>&2; exit 1; }
nolistBAM() { echo "## ERROR: You need to provide the text file containing the list of BAM files and samples ID through --listBAM option. Exit." 1>&2; exit 1; }
notargetsBED() { echo "## ERROR: You need to provide the BED file containing the on-targets and off-targets through --targetsBED option. Exit." 1>&2; exit 1; }
errorlistBAM() { echo "## ERROR: The list of BAM files and samples ID provided through --listBAM should have 2 columns for BAM files and sample IDs and be tab delimited. Exit." 1>&2; exit 1; }

while getopts ":-:" o; do
    case "${o}" in
    	-)  
        	case $OPTARG in
                work)
                    work="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ;;
                listBAM)
                    listBAM="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ;;
                mosdepth)
                    mosdepth="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ;;
                targetsBED)
                    targetsBED="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
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

if [ -z "${work}" ]; then
    nowork
fi
if [ -z "${listBAM}" ]; then
    nolistBAM
fi
if [ -z "${mosdepth}" ]; then
    nomosdepth
fi
if [ -z "${targetsBED}" ]; then
    notargetsBED
fi

sed '/^$/d' $listBAM > $listBAM.temp

mkdir -p $work
cut -f1-4 $targetsBED > $work/targetsBED.bed

echo "Step 1: Using mosdepth to extract coverage for each on-target and off-target"
if [[ `awk -F"\t" 'NF==1' $listBAM` ]]; then
    echo "   Only one field detected. ID will be deduced from the BAM file name"
fi

# check that all BAM files exists
while read file
do
	a=$(echo $file | awk -F" " '{print NF}')
	if [[ "$a" == 1 ]]; then
		pat=$(echo $file | cut -f1 -d" " | awk -F"\t" '{n=split($1,a,"/"); split(a[n],b,".bam"); print b[1]}')
		bam=$(echo $file | cut -f1 -d" ")
	fi
	if [ $a -gt 1 ]; then
		pat=$(echo $file | cut -f2 -d" ")
		bam=$(echo $file | cut -f1 -d" ")
	fi
	if [ ! -f $bam ]; then
		rm $listBAM.temp
		echo "## ERROR: BAM file $bam not found. Exit." 1>&2; exit 1;
	fi
done < $listBAM.temp

rm -f $work/list.txt
touch $work/list.txt
# coverage from the bam files with mosdepth 
while read file
do
	a=$(echo $file | awk -F" " '{print NF}')
	if [[ "$a" == 1 ]]; then
		pat=$(echo $file | cut -f1 -d" " | awk -F"\t" '{n=split($1,a,"/"); split(a[n],b,".bam"); print b[1]}')
		bam=$(echo $file | cut -f1 -d" ")
		echo $pat
	fi
	if [ $a -gt 1 ]; then
		echo "more than 1 field"
		pat=$(echo $file | cut -f2 -d" ")
		bam=$(echo $file | cut -f1 -d" ")
	fi
	echo $pat >> $work/list.txt
	if [ ! -f $work/coverage/$pat.cov.tsv ] || [ ! -s $work/coverage/$pat.cov.tsv ]; then
		echo "   Processing sample $pat"
		# downloaded here: https://github.com/brentp/mosdepth/releases/download/v0.3.2/mosdepth
		# -n = no-per-base; -t = threads; -Q mapq; precision is for rounding of numbers
		# gives the average per base coverage
		export MOSDEPTH_PRECISION=10
		$mosdepth  -n -t 2 -Q 50 -b $work/targetsBED.bed $work/$pat $bam
		zcat $work/$pat.regions.bed.gz > $work/$pat.regions.bed
		# remultiply the per-base average by the size to have the total coverage
		mkdir -p $work/coverage
		awk -F"\t" '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5*($3-$2)}' $work/$pat.regions.bed > $work/coverage/$pat.cov.tsv
		# cleaning
		rm $work/$pat.regions.bed.gz $work/$pat.regions.bed $work/$pat.mosdepth* $work/$pat.regions.bed.gz.csi
	fi
done < $listBAM.temp

rm $work/targetsBED.bed $listBAM.temp

echo "Step 2: Merging individual coverage files into one file"
# processing of coverage files to merge them and add header
printf "Chr\tbegin\tend\tname\tGC" > $work/header.target.tsv
cut -f1-5 $targetsBED  > $work/target.only.tsv
nb=$(wc -l $work/target.only.tsv | cut -f1 -d" ")

for pat in $(cat $work/list.txt)
do
	echo "   Processing sample $pat"
	nb2=$(wc -l $work/coverage/$pat.cov.tsv | cut -f1 -d" ")
	if [[ "$nb" != "$nb2" ]]; then
		echo "## ERROR: coverage file coverage/$pat.cov.tsv is incomplete. Please check the file and rerun the script. Exit." 1>&2; exit 1;
	fi
	cut -f5 $work/coverage/$pat.cov.tsv > $work/$pat.cov.only.tsv
	cat $work/target.only.tsv > $work/temp.tsv
	paste -d"\t" $work/temp.tsv $work/$pat.cov.only.tsv > $work/target.only.tsv
	printf "\t$pat" >> $work/header.target.tsv
done

echo "" >> $work/header.target.tsv
cat $work/header.target.tsv $work/target.only.tsv > $work/ALL.target.tsv

rm $work/temp.tsv $work/header.target.tsv $work/*.only.tsv $work/list.txt

