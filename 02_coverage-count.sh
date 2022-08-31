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
if [[ `awk -F"\t" 'NF==1' $listBAM` ]]; then
    errorlistBAM
fi

sed -i '/^$/d' $listBAM

cut -f1-4 $targetsBED > $work/targetsBED.bed

echo "Step 1: Using mosdepth to extract coverage for each on-target and off-target"
# coverage from the bam files with mosdepth 
while read file
do
	pat=$(echo $file | cut -f2 -d" ")
	bam=$(echo $file | cut -f1 -d" ")
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
done < $listBAM
rm $work/targetsBED.bed

echo "Step 2: Merging individual coverage files into one file"
# processing of coverage files to merge them and add header
printf "Chr\tbegin\tend\tname\tGC" > $work/header.target.tsv
cut -f1-5 $targetsBED  > $work/target.only.tsv

for pat in $(cut -f2 $listBAM)
do
	echo "   Processing sample $pat"
	cut -f5 $work/coverage/$pat.cov.tsv > $work/$pat.cov.only.tsv
	cat $work/target.only.tsv > $work/temp.tsv
	paste -d"\t" $work/temp.tsv $work/$pat.cov.only.tsv > $work/target.only.tsv
	printf "\t$pat" >> $work/header.target.tsv
done

echo "" >> $work/header.target.tsv
cat $work/header.target.tsv $work/target.only.tsv > $work/ALL.target.tsv

rm $work/temp.tsv $work/header.target.tsv $work/*.only.tsv

