
# cn.mops

here=/home/mquinodo/SYNO
folder=$here/scripts_NGS_analysis/OFF-PEAK
batch=unsolved-2020

here2=/home/mquinodo/$batch
a=$(ls $here2/*bam | cut -d"/" -f5 | cut -d"." -f1 | tr '\n' ',')
patients=${a::-1}
script=$folder/scripts-tools/cn.mops
output=$folder/$batch/cn.mops
filetype=".dup.BQSR.bam"
mkdir -p $output

Rscript $script/cn.mops.R $here2 $patients $filetype $output > $output/part1_log.tsv 2>&1

awk -F"," '{split($7,a,"."); split($10,b,"CN"); print $2 "\t" $3 "\t" $4"\t" a[1] "\t" $8 "\t" b[2]}' $output/cnmops.cnvs.csv | tr -d '"' | tail -n+2 > $output/$batch.cnmops.out.tsv


# CNVkit

here=/home/mquinodo/SYNO
batch=unsolved-2020
lab=WES/EXOMES/$batch
nametarget=Twist_Exome_Target_hg19
targetfile=$here/NGS/Twist/Twist_Exome_Target_hg19.bed
folder=$here/scripts_NGS_analysis/OFF-PEAK

output=$folder/$batch/cnvkit
mkdir -p $output
tools_dir=/opt/Tools/NGS_big_data
ref=$tools_dir/hg19.p13.plusMT.no_alt_analysis_set.fa

a=$(cat $folder/$batch/$batch.bam_list.txt | cut -f1 | cut -d"/" -f8 | tr '\n' ',')
b=${a::-1}
IFS=',' read -r -a patients <<< "$b"

a=$(cat $folder/$batch/$batch.bam_list.txt | cut -f1 | tr '\n' ',')
b=${a::-1}
IFS=',' read -r -a bams <<< "$b"

# autobin
rm $output/combine-command.sh
touch $output/combine-command.sh
printf '%s' "/home/mquinodo/cnvkit/cnvkit.py autobin " >> $output/combine-command.sh
for bam in "${bams[@]}"
do
printf "$bam " >> $output/combine-command.sh
done
printf '%s' "-t $targetfile -g /opt/Tools/NGS_big_data/access-5k-mappable.hg19.bed --annotate /opt/Tools/NGS_big_data/refFlat.txt" >> $output/combine-command.sh
cd $output
bash $output/combine-command.sh > $output/part1_log.tsv 2>&1
cd

# coverage
for bam in "${bams[@]}"
do
	pat=$(echo $bam | cut -d"/" -f8)
	/home/mquinodo/cnvkit/cnvkit.py coverage $bam $output/$nametarget.target.bed -o $output/$pat.targetcoverage.cnn &
	/home/mquinodo/cnvkit/cnvkit.py coverage $bam $output/$nametarget.antitarget.bed -o $output/$pat.antitargetcoverage.cnn
done

# reference
n=0
for i in "${patients[@]}"
do
	rm $output/combine-command.sh
	touch $output/combine-command.sh
	printf '%s' "/home/mquinodo/cnvkit/cnvkit.py reference " >> $output/combine-command.sh
	pat=${patients[$n]}
	controls=("${patients[@]}") 
	for con in "${controls[@]}"
	do
	if [ "$pat" != "$con" ]; then	
	printf "$output/$con.targetcoverage.cnn " >> $output/combine-command.sh
	printf "$output/$con.antitargetcoverage.cnn " >> $output/combine-command.sh
	fi
	done
	printf '%s' "--fasta $ref -o $output/$pat.ref.coverage.cnn" >> $output/combine-command.sh
	bash $output/combine-command.sh

	n=$((n+1))
done

# fix
for pat in "${patients[@]}"
do
	/home/mquinodo/cnvkit/cnvkit.py fix $output/$pat.targetcoverage.cnn $output/$pat.antitargetcoverage.cnn $output/$pat.ref.coverage.cnn -o $output/$pat.cnr
done

# segment
for pat in "${patients[@]}"
do
	/home/mquinodo/cnvkit/cnvkit.py segment $output/$pat.cnr -o $output/$pat.cns
done

# call
for pat in "${patients[@]}"
do
	/home/mquinodo/cnvkit/cnvkit.py call $output/$pat.cns -o $output/$pat.call.tsv
	awk -F"\t" -v pat=$pat '{print $1 "\t" $2 "\t" $3 "\t" pat "\t" $5 "\t" $6}' $output/$pat.call.tsv > $output/$pat.call2.tsv
done

cat $output/*.call2.tsv | awk -F"\t" '{if($6!=2) print $0}' | grep -v -P "chrY|chromosome|gl" > $output/$batch.cnvkit.out.tsv


# CODEX2

here=/home/mquinodo/SYNO
folder=$here/scripts_NGS_analysis/OFF-PEAK
batch=unsolved-2020
bedFile=$here/NGS/Twist/Twist_Exome_Target_hg19.bed
grep -v -P "chrUn|chrX|chrY" $bedFile > $bedFile.cano.bed

script=$folder/scripts-tools/codex2
here2=/home/mquinodo/$batch
a=$(ls $here2/*bam | cut -d"/" -f5 | cut -d"." -f1 | tr '\n' ',')
patients=${a::-1}
filetype=".dup.BQSR.bam"
output=$folder/$batch/codex2
mkdir -p $output

Rscript $script/codex2.R $here2 $patients $filetype $output $bedFile.cano.bed > $output/part1_log.tsv 2>&1

awk -F"\t" '{if($6>$5) print $3 "\t" $5 "\t" $6 "\t" $2 "\t" $13 "\t" $12}' $output/output.tsv | tail -n+2 > $output/$batch.codex2.out.tsv


# CoNIFER

here=/home/mquinodo/SYNO
folder=$here/scripts_NGS_analysis/OFF-PEAK
batch=unsolved-2020
targetfile=$here/NGS/Twist/Twist_Exome_Target_hg19.bed
here2=/home/mquinodo/$batch

script=$folder/scripts-tools/conifer
output=$folder/$batch/conifer
mkdir -p $output
sed 's/chr//g' $targetfile | grep -v -P "_|Y|Un" > $targetfile.temp
cat $script/header.txt $targetfile.temp > $output/probes.txt
rm $targetfile.temp

# From this point, pulled docker image
sudo docker run -v $output:/output -v $here2:/$batch -it ielis/conifer

batch=unsolved-2020

a=$(ls /${batch}/*bam | cut -d"/" -f3 | cut -d"." -f1 | tr '\n' ',')
b=${a::-1}
IFS=',' read -r -a patients <<< "$b"

# 1) Create RPKM files
mkdir -p /output/rpkm
for sam in "${patients[@]}"
do
	python conifer.py rpkm \
	--probes /output/probes.txt \
	--input /${batch}/${sam}.dup.BQSR.bam \
	--output /output/rpkm/${sam}.rpkm.txt
done

# 2) De-noising using SVD, while poorly responding probes are masked, and the Z- and SVD- transformations are applied. The final SVD-ZRPKM values for each exon are stored in a unified HDF5 file.
mkdir -p /output/results
svd=1 # Start point

python conifer.py analyze \
--probes /output/probes.txt \
--rpkm_dir /output/rpkm \
--svd ${svd} \
--write_svals /output/results/${batch}.svd${svd}.singular_values.txt \
--plot_scree /output/results/${batch}.svd${svd}.screeplot.png \
--write_sd /output/results/${batch}.svd${svd}.sd_values.txt \
--output /output/results/${batch}.svd${svd}.analysis.hdf5

# 3) CNV calls using an algorithm based on SVD-ZRPKM thresholds
python conifer.py call \
--input /output/results/${batch}.svd${svd}.analysis.hdf5 \
--output /output/results/${batch}.svd${svd}.calls.txt

awk -F"\t" '{split($1,a,"."); print $2 "\t" $3 "\t" $4 "\t" a[1] "\t" $5 "\t" $5}' /output/results/${batch}.svd${svd}.calls.txt | sed 's/del/1/g' | sed 's/dup/3/g' | tail -n+2 > /output/${batch}.conifer.out.tsv


# ExomeDepth

here=/home/mquinodo/SYNO
folder=$here/scripts_NGS_analysis/OFF-PEAK-
batch=unsolved-2020
here2=/home/mquinodo/$batch

script=$folder/scripts-tools/ExomeDepth
a=$(ls $here2/*bam | cut -d"/" -f5 | cut -d"." -f1 | tr '\n' ',')
patients=${a::-1}
filetype=".dup.BQSR.bam"
output=$folder/$batch/ExomeDepth
mkdir -p $output

Rscript $script/ExomeDepth.R $here2 $patients $filetype $output > $output/part1_log.tsv 2>&1

grep chr $output/exome_calls_*.csv > $output/all.csv
grep -v "reads.expected" $output/all.csv | $script/csv2tab.py | awk -F"\t" '{split($1,a,"/"); n=split(a[9],b,"_"); split(b[n],c,".csv"); print a[7] "\t" c[1] "\t" $3 "\t" $7 "\t" $5 "\t" $6 "\t" $9 "\t" $12 "\t" $4 "\t" $13 "\t" $14}' > $output/all.temp.tsv
cat $script/head.tsv $output/all.temp.tsv > $output/all.tsv
awk -F"\t" '{print "chr" $4 "\t" $5 "\t" $6 "\t" $2 "\t" $7 "\t" $8*2}' $output/all.temp.tsv | grep -v -P "chr\t" | tail -n+2 > $output/$batch.ExomeDepth.out.tsv
rm $output/all.temp*


# ControlFREEC

here=/home/mquinodo/SYNO
folder=$here/scripts_NGS_analysis/OFF-PEAK-
batch=unsolved-2020
here2=/home/mquinodo/$batch

script=$folder/scripts-tools/FREEC
output=$folder/$batch/FREEC
mkdir -p $output
a=$(ls $here2/*bam | cut -d"/" -f5 | cut -d"." -f1 | tr '\n' ',')
b=${a::-1}
IFS=',' read -r -a patients <<< "$b"

for patient in "${patients[@]}"
do
	sed "s|oouutt|$output|g" $script/config_default-2020.txt | sed "s|iinn|$here2\/$patient.dup.BQSR.bam|g" > $output/$patient.config.txt
	/home/mquinodo/FREEC-11.6/src/freec -conf $output/$patient.config.txt
done

for patient in "${patients[@]}"
do
	cat /home/mquinodo/FREEC-11.6/scripts/assess_significance.R | R --slave --args $output/$patient.dup.BQSR.bam_CNVs $output/$patient.dup.BQSR.bam_ratio.txt
	awk -F"\t" -v patient="$patient" '{print "chr" $1 "\t" $2 "\t" $3 "\t" patient "\t" $6 "\t" $4}' $output/$patient.dup.BQSR.bam_CNVs.p.value.txt > $output/$patient.dup.BQSR.bam_CNVs.p.value.pro.txt 
done

cat $output/*.dup.BQSR.bam_CNVs.p.value.pro.txt | grep -v "WilcoxonRankSumTestPvalue" > $output/$batch.controlFREEC.out.tsv


# GATK gCNV

here=/home/mquinodo/SYNO
folder=$here/scripts_NGS_analysis/OFF-PEAK-
batch=unsolved-2020
here2=/home/mquinodo/$batch
targetfile=$here/NGS/Twist/Twist_Exome_Target_hg19.bed
tools_dir=/opt/Tools/NGS_big_data

gatk=/home/mquinodo/gatk-4.2.6.1/gatk-package-4.2.6.1-local.jar
ref=$tools_dir/hg19.p13.plusMT.no_alt_analysis_set.fa
output=$folder/$batch/gatk
mkdir -p $output
a=$(ls $here2/*bam | cut -d"/" -f5 | cut -d"." -f1 | tr '\n' ',')
b=${a::-1}
IFS=',' read -r -a patients <<< "$b"

source activate gatk

# First to preprocess
# This only has to be run once per capture kit (here, Twist capture kit)

java -Xmx8G -jar $gatk PreprocessIntervals \
-R $ref \
--padding 0 --bin-length 0 \
-L $targetfile \
-imr OVERLAPPING_ONLY \
-O $output/targets.preprocessed.interval_list

### Collect read counts per target region
for pat in "${patients[@]}"
do
java -Xmx8G -jar $gatk CollectReadCounts \
-R $ref \
-imr OVERLAPPING_ONLY \
-L $output/targets.preprocessed.interval_list \
-I $here2/$pat.dup.BQSR.bam \
-O $output/$pat.hdf5
done 

### Annotate intervals with features and subset regions of interest with FilterIntervals
java -Xmx8G -jar $gatk AnnotateIntervals \
-R $ref  \
-L $output/targets.preprocessed.interval_list \
-imr OVERLAPPING_ONLY \
-O $output/gc.annotated.tsv

### Filter reads based on GC-content (remove extremes), here default parameters of GATK were used
rm $output/combine-command.sh
touch $output/combine-command.sh
printf "java -Xmx8G -jar $gatk FilterIntervals -L $output/targets.preprocessed.interval_list --annotated-intervals $output/gc.annotated.tsv --minimum-gc-content 0.1 --maximum-gc-content 0.9 --minimum-mappability 0.9 --maximum-mappability 1.0 --minimum-segmental-duplication-content 0.0 --maximum-segmental-duplication-content 0.5 --low-count-filter-count-threshold 5 --low-count-filter-percentage-of-samples 90.0 --extreme-count-filter-minimum-percentile 1.0 --extreme-count-filter-maximum-percentile 99.0 --extreme-count-filter-percentage-of-samples 90.0 --interval-merging-rule OVERLAPPING_ONLY -O $output/$batch.gcfiltered.interval_list " >> $output/combine-command.sh

for pat in "${patients[@]}"
do
	printf '%s' "-I $output/$pat.hdf5 " >> $output/combine-command.sh
done

bash $output/combine-command.sh
	
### Determine Contig Ploidy
rm $output/combine-command.sh
touch $output/combine-command.sh
printf "java -Xmx8G -jar $gatk DetermineGermlineContigPloidy -L $output/$batch.gcfiltered.interval_list --contig-ploidy-priors ${tools_dir}/ploidy_priors.tsv --interval-merging-rule OVERLAPPING_ONLY --output-prefix $batch --output $output " >> $output/combine-command.sh

for pat in "${patients[@]}"
do
	printf '%s' "-I $output/$pat.hdf5 " >> $output/combine-command.sh
done

bash $output/combine-command.sh


rm $output/combine-command.sh
touch $output/combine-command.sh
printf "java -Xmx8G -jar $gatk GermlineCNVCaller --run-mode COHORT -L $output/$batch.gcfiltered.interval_list --contig-ploidy-calls $output/${batch}-calls --annotated-intervals $output/gc.annotated.tsv --interval-merging-rule OVERLAPPING_ONLY --output-prefix ${batch} -O $output " >> $output/combine-command.sh

for pat in "${patients[@]}"
do
printf '%s' "-I $output/$pat.hdf5 " >> $output/combine-command.sh
done

bash $output/combine-command.sh

n=11
for ind in $(seq 0 1 $((n-1)))
do
	java -Xmx8G -jar $gatk PostprocessGermlineCNVCalls \
	--model-shard-path $output/${batch}-model \
	--calls-shard-path $output/${batch}-calls \
	--sequence-dictionary $tools_dir/hg19.p13.plusMT.no_alt_analysis_set.dict \
	--allosomal-contig chrX \
	--contig-ploidy-calls $output/${batch}-calls \
	--sample-index $ind \
	--output-genotyped-intervals $output/${patients[$ind]}_${batch}_cohort.vcf.gz \
	--output-genotyped-segments $output/${patients[$ind]}_${batch}_segments_cohort.vcf.gz \
	--output-denoised-copy-ratios $output/${patients[$ind]}_${batch}_denoised_copy_ratios.tsv
done

for ind in $(seq 0 1 $((n-1)))
do
	gunzip $output/${patients[$ind]}_${batch}_segments_cohort.vcf.gz
	grep -v -P "#|chrY|gl" $output/${patients[$ind]}_${batch}_segments_cohort.vcf | awk -F"\t" -v ind=${patients[$ind]} '{split($10,b,":"); if(b[2]!=2) {split($8,a,"="); print $1 "\t" $2 "\t" a[2] "\t" ind "\t" $6 "\t" b[2]}}' > $output/${patients[$ind]}.out.vcf
	grep -v -P "#|chrY|gl" $output/${patients[$ind]}_${batch}_segments_cohort.vcf | awk -F"\t" -v ind=${patients[$ind]} '{split($10,b,":"); if(b[2]!=2 && b[4]>100) {split($8,a,"="); print $1 "\t" $2 "\t" a[2] "\t" ind "\t" $6 "\t" b[2]}}' > $output/${patients[$ind]}.out2.vcf
done

cat $output/*.out.vcf > $output/$batch.gatk.out.tsv
cat $output/*.out2.vcf > $output/$batch.gatk2.out.tsv


# SavvyCNV

here=/home/mquinodo/SYNO
folder=$here/scripts_NGS_analysis/OFF-PEAK-
batch=unsolved-2020
here2=/home/mquinodo/$batch

script=$folder/scripts-tools/SavvyCNV
output=$folder/$batch/SavvyCNV
mkdir -p $output
a=$(ls $here2/*bam | cut -d"/" -f5 | cut -d"." -f1 | tr '\n' ',')
b=${a::-1}
IFS=',' read -r -a patients <<< "$b"

for pat in "${patients[@]}"
do
	java -cp /usr/local/bin/gatk-package-4.1.8.1-local.jar:/home/mquinodo/SavvySuite -XX:ParallelGCThreads=2 -XX:ConcGCThreads=2 -Xmx1g CoverageBinner $here2/$pat.dup.BQSR.bam > $output/$pat.coverageBinner
done

# -threshold For each 200bp chunk of the genome, the average number of reads  is calculated, and a chunk is counted as on-target if it is above this threshold, and off-target if it is below. The default is 5, and it can be a floating point number.
export CLASSPATH=/usr/local/bin/gatk-package-4.1.8.1-local.jar:/home/mquinodo/SavvySuite
java -Xmx30g CoverageOffTarget -readLength 100 -threshold 5 $output/*.coverageBinner > $output/CoverageOffTarget.threshold5.tsv 2> $output/log_messages-CoverageOffTarget-5.txt

##### ADJUST VALUES FROM LOG
ON=200
OFF=19800

java -cp /usr/local/bin/gatk-package-4.1.8.1-local.jar:/home/mquinodo/SavvySuite -Xmx30g SavvyCNV -d $OFF $output/*.coverageBinner > $output/cnv_list-off-$OFF.tsv 2> $output/log_messages-SavvyCNV-OFF.txt
java -cp /usr/local/bin/gatk-package-4.1.8.1-local.jar:/home/mquinodo/SavvySuite -Xmx30g SavvyCNV -d $ON $output/*.coverageBinner > $output/cnv_list-on-$ON.tsv 2> $output/log_messages-SavvyCNV-ON.txt

for pat in "${patients[@]}"
do
	grep $pat.coverageBinner $output/cnv_list-on-$ON.tsv | awk -F"\t" -v pat="$pat" '{print $1 "\t" $2 "\t" $3 "\t" pat "\t" $8 "\t" $9*2}' > $output/$pat.$ON.out.tsv
	grep $pat.coverageBinner $output/cnv_list-off-$OFF.tsv | awk -F"\t" -v pat="$pat" '{print $1 "\t" $2 "\t" $3 "\t" pat "\t" $8 "\t" $9*2}' > $output/$pat.$OFF.out.tsv
done

cat $output/*.$ON.out.tsv > $output/$batch.SavvyCNVON.out.tsv 
cat $output/*.$OFF.out.tsv > $output/$batch.SavvyCNVOFF.out.tsv 


# OFF-PEAK

###########################################################
###### PART 1: creation of targets and off-targets ########
###########################################################

here=/home/mquinodo/SYNO/scripts_NGS_analysis/OFF-PEAK-train5
script=$here/scripts/01_targets-offtargets.sh
targets=$here/refs/Twist_Exome_Target_hg19.bed
ref=/opt/Tools/NGS_big_data/hg19.p13.plusMT.no_alt_analysis_set.fa
out=Twist1
bash $script --genome hg19 --targets $targets --out $out --ref $ref

here=/home/mquinodo/SYNO/scripts_NGS_analysis/OFF-PEAK-train5
script=$here/scripts/01_targets-offtargets.sh
# Parameters
minOntarget=30
maxOntarget=75
targets=$here/refs/Twist_Exome_Target_hg19.bed
ref=/opt/Tools/NGS_big_data/hg19.p13.plusMT.no_alt_analysis_set.fa
out=Twist1-small
bash $script --genome hg19 --minOntarget $minOntarget --maxOntarget $maxOntarget --targets $targets --out $out --ref $ref

###########################################################
###### PART 2: take coverage from BAM files ###############
###########################################################

here=/home/mquinodo/SYNO/scripts_NGS_analysis/OFF-PEAK-train5
batch=unsolved-2020
script=$here/scripts/02_bam-count.sh
work=$here/$batch
listBAM=$work/$batch.bam_list.txt
targetsBED=$here/refs/Twist1.bed
mosdepth=/home/mquinodo/mosdepth
bash $script --work $work --listBAM $listBAM --mosdepth $mosdepth --targetsBED $targetsBED

here=/home/mquinodo/SYNO/scripts_NGS_analysis/OFF-PEAK-train5
batch=unsolved-2020-small
script=$here/scripts/02_bam-count.sh
work=$here/$batch
listBAM=$work/$batch.bam_list.txt
targetsBED=$here/refs/Twist1-small.bed
mosdepth=/home/mquinodo/mosdepth
bash $script --work $work --listBAM $listBAM --mosdepth $mosdepth --targetsBED $targetsBED

###########################################################
###### PART 3: testing on the batches #####################
###########################################################

here=/home/mquinodo/SYNO/scripts_NGS_analysis/OFF-PEAK-train5
batch=unsolved-2020
script=$here/scripts/02_R-PCA_testing14.R
data=$here/$batch/ALL.target.tsv
databasefile=$here/refs/data-hg19.RData
Rscript $script --output $here/$batch --data $data --databasefile $databasefile &> $here/$folder/log.tsv

here=/home/mquinodo/SYNO/scripts_NGS_analysis/OFF-PEAK-train5
batch=unsolved-2020-small
script=$here/scripts/02_R-PCA_testing14.R
data=$here/$batch/ALL.target.tsv
databasefile=$here/refs/data-hg19.RData
Rscript $script --output $here/$batch --data $data --databasefile $databasefile &> $here/$folder/log.tsv

###########################################################
###### PART 4: plots ######################################
###########################################################

here=/home/mquinodo/SYNO/scripts_NGS_analysis/OFF-PEAK-train5
script=$here/scripts/02_R-PCA_testing14-plot.R
pat="LL237"
chr="chr6"
begin=64694275
end=64791895
side=20
data=$here/unsolved-2020/05_RData-files/data-$pat.RData
out=$here/unsolved-2020/$pat
databasefile=$here/refs/data-hg19.RData
Rscript $script --ID $pat --chr $chr --begin $begin --end $end --side $side --data $data --out $out --databasefile $databasefile


