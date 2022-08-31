#!/bin/bash

### Made by Mathieu Quinodoz
### August 2022
### the same commands were used to analyse the samples from the retinal disease cohort

here=/home/mquinodo
batch=ICR96
folder=$here/SYNO/WES/$batch

# MAPPING AND VARIANT CALLING
for pat in {17296..17343} {17356..17403}
do
	ref=/opt/Tools/NGS_big_data/hg19.p13.plusMT.no_alt_analysis_set.fa
	fq1=$here/$batch/$pat/${pat}_1.fastq.gz
	fq2=$here/$batch/$pat/${pat}_2.fastq.gz
	cp $folder/$pat/original_data/${pat}_1.fastq.gz $here/$batch/$pat
	cp $folder/$pat/original_data/${pat}_2.fastq.gz $here/$batch/$pat
	# alignement
	export PATH="$PATH:/usr/local/bin/bwa-0.7.15"
	RG=$(echo "@RG\tID:$pat\tSM:$pat\tPL:Illumina")
	bwa mem -M -C -t 32 -R $RG $ref $fq1 $fq2 > $here/$batch/$pat/$pat.sam 2> $here/$batch/$pat/bwa.$pat.log
	rm $here/$batch/$pat/*fastq*
done


NAS=$here/SYNO/WES

tools_dir=/opt/Tools/NGS_big_data
gatk=/usr/local/bin/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar

target=$tools_dir/Twist_ComprehensiveExome_targets_hg19.1000pad.bed
ref=$tools_dir/hg19.p13.plusMT.no_alt_analysis_set.fa
db1=$tools_dir/dbsnp_138.hg19.vcf
omni=$tools_dir/hg19_v0_1000G_omni2.5.b37.chr.vcf.gz
G1000=$tools_dir/hg19_v0_1000G_phase1.snps.high_confidence.b37.chr.vcf.gz
hapmap=$tools_dir/hg19_v0_hapmap_3.3.b37.chr.vcf.gz
mills=$tools_dir/hg19_v0_Mills_and_1000G_gold_standard.indels.b37.chr.vcf.gz
indeldb=$tools_dir/hg19_v0_Homo_sapiens_assembly19.known_indels.chr.vcf.gz
axiom=$tools_dir/hg19_v0_Axiom_Exome_Plus.genotypes.all_populations.poly.chr.vcf.gz

# Combine all g.vcf by chromosome

for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrMT
do

	grep -P "$chr\t" $target > $here/$batch/target.$chr.bed

	rm $here/$batch/combine-command.$chr.sh
	touch $here/$batch/combine-command.$chr.sh
	printf "java -jar $gatk CombineGVCFs -G StandardAnnotation -R $ref -L $here/$batch/target.$chr.bed -O $here/$batch/$batch.merged.g.$chr.vcf " >> $here/$batch/combine-command.$chr.sh
	# -Xms32g -Xmx32g

	for file in $here/$batch/*/*.g.vcf; do
		printf '%s' "-V $file " >> $here/$batch/combine-command.$chr.sh
	done

	bash $here/$batch/combine-command.$chr.sh 2> $here/$batch/CombineGVCFs.$chr.$batch.log

	# Joint Genotyping by chromosome
	java -Xmx4g -jar $gatk GenotypeGVCFs \
	-R $ref \
	-V $here/$batch/$batch.merged.g.$chr.vcf \
	-O $here/$batch/$batch.merged.$chr.vcf \
	-D $db1 \
	-G StandardAnnotation \
	-L $here/$batch/target.$chr.bed 2> $here/$batch/GenotypeGVCFs.$chr.$batch.log

	rm $here/$batch/$batch.merged.g.$chr.vcf

done

# Merging VCFs by chromosome together

rm $here/$batch/merge.vcf.sh
touch $here/$batch/merge.vcf.sh
printf "bcftools concat --threads 10 " >> $here/$batch/merge.vcf.sh

for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrMT
do
printf '%s' "$here/$batch/$batch.merged.$chr.vcf " >> $here/$batch/merge.vcf.sh
done

bash $here/$batch/merge.vcf.sh > $here/$batch/$batch.merged.vcf 2> $here/$batch/MergeVCFs.$batch.log

rm $here/$batch/$batch.merged.g.chr*.vcf* $here/$batch/$batch.merged.chr*.vcf* $here/$batch/GenotypeGVCFs.chr*.$batch.log
rm $here/$batch/CombineGVCFs.chr*.$batch.log $here/$batch/combine-command.chr*.sh $here/$batch/target.*.bed


# VQSR

java -jar $gatk VariantRecalibrator \
-R $ref \
-V $here/$batch/$batch.merged.vcf \
-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
--resource:mills,known=false,training=true,truth=true,prior=12 $mills \
--resource:axiomPoly,known=false,training=true,truth=false,prior=10 $axiom \
--resource:dbsnp,known=true,training=false,truth=false,prior=2 $indeldb \
-O $here/$batch/$batch.merged.INDEL.recal \
--tranches-file $here/$batch/$batch.merged.INDEL.tranches \
-mode INDEL \
--output-model $here/$batch/$batch.merged.VQSR-INDEL.model 2> $here/$batch/VariantRecalibrator-INDEL.$batch.log &


java -jar $gatk VariantRecalibrator \
-R $ref \
-V $here/$batch/$batch.merged.vcf \
-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
--resource:hapmap,known=false,training=true,truth=true,prior=15 $hapmap \
--resource:omni,known=false,training=true,truth=true,prior=12 $omni \
--resource:1000G,known=false,training=true,truth=false,prior=10 $G1000 \
--resource:dbsnp,known=true,training=false,truth=false,prior=7 $db1 \
-O $here/$batch/$batch.merged.SNP.recal \
--tranches-file $here/$batch/$batch.merged.SNP.tranches \
-mode SNP \
--output-model $here/$batch/merged.VQSR-SNP.model 2> $here/$batch/VariantRecalibrator-SNP.$batch.log


java -jar $gatk ApplyVQSR \
-O $here/$batch/$batch.merged.SNP.VQSR.vcf \
-V $here/$batch/$batch.merged.vcf \
--recal-file $here/$batch/$batch.merged.SNP.recal \
--tranches-file $here/$batch/$batch.merged.SNP.tranches \
--truth-sensitivity-filter-level 99.0 \
--create-output-variant-index true \
-mode SNP 2> $here/$batch/ApplyVQSR-SNP.$batch.log


java -jar $gatk ApplyVQSR \
-O $here/$batch/$batch.merged.SNP-INDEL.VQSR.vcf \
-V $here/$batch/$batch.merged.SNP.VQSR.vcf \
--recal-file $here/$batch/$batch.merged.INDEL.recal \
--tranches-file $here/$batch/$batch.merged.INDEL.tranches \
--truth-sensitivity-filter-level 99.0 \
--create-output-variant-index true \
-mode INDEL 2> $here/$batch/ApplyVQSR-INDEL.$batch.log

# VQSRTrancheSNP99.90to100.00 means that the variant was in the range of VQSLODs corresponding to the remaining 0.1% of the truth set, which are considered false positives.


# Postprocessing of VCF file

 for pat in `bcftools query -l $here/$batch/$batch.merged.SNP-INDEL.VQSR.vcf`; do
   bcftools view -c1 -Ov -s $pat -o $here/$batch/$pat.VQSR.vcf $here/$batch/$batch.merged.SNP-INDEL.VQSR.vcf
   bcftools norm -m -both --fasta-ref $ref -o $here/$batch/$pat.VQSR.norm.vcf $here/$batch/$pat.VQSR.vcf
   mkdir -p $NAS/$batch/$pat/variants
   grep -P -v "AC=0;" $here/$batch/$pat.VQSR.norm.vcf | grep -P -v "\t\*\t" > $NAS/$batch/$pat/variants/$pat.VQSR.norm.noref.vcf
   rm $here/$batch/$pat.*.vcf 
 done


 # copy back to the NAS

 mkdir -p $NAS/$batch/10_joint-calling
 cp $here/$batch/*log $NAS/$batch/10_joint-calling
 cp $here/$batch/$batch.merged.SNP-INDEL.VQSR.vcf* $NAS/$batch/10_joint-calling

 rm $here/$batch/*log $here/$batch/$batch.*.vcf* $here/$batch/$batch.*SNP* $here/$batch/$batch.*INDEL* $here/$batch/merged.VQSR-SNP.model $here/$batch/merge.vcf.sh


