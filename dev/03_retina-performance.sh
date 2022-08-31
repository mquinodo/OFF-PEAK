

here=/home/mquinodo/SYNO/scripts_NGS_analysis/OFF-PEAK
here2=/home/mquinodo/SYNO/WES/EXOMES
batch=all

work=$here/$batch
mkdir -p $work

a=$(cat $here/unsolved-2020/unsolved-2020.bam_list.txt $here/unsolved-2021/unsolved-2021.bam_list.txt $here/unsolved-2022/unsolved-2022.bam_list.txt  | cut -f1 | cut -d"/" -f8 | tr '\n' ',')
b=${a::-1}
IFS=',' read -r -a patients <<< "$b"

# selecting phenotypes and then removing unsure ones
grep -P "Stargardt_disease|Cone-rod_dystrophy|retinitis_pigmentosa|Retinitis_pigmentosa|Retinal_dystrophy|Usher_syndrome|Choroideremia|CLN3|cone_dystrophy|retinal_dystrophy|Leber_congenital_amaurosis" /home/mquinodo/SYNO/NGS/OMIM/OMIM_2022.08.10/gene_to_disease.txt | grep -v -P "\?Leber|corneoretinal_dystrophy|\?Mic|\?Retinitis_pigmentosa|HK1\t|\?Retinal|\?Mac|\?Cone|\?Usher_|GDF6" > $work/retina-all.tsv
# selecting dominant and then removing 1° dominant but not retina 2° dominant unsure 3° dominant only missense
grep -P "RPGR|ominant" $work/retina-all.tsv | grep -v -P "DHDDS|ZNF408|CDH23|HARS1|MYO7A|PPP2R3C|RBP4|PDE6B|CRB1|PDE6H|ABCA4|ROM1|ADGRV1|PDZD7" | grep -v -P "ARL3|GUCA1A|IMPDH1|KIF3B|SEMA4A" | grep -v -P "AIPL1|IMPG2|KLHL7|NR2E3|PITPNM3|RLBP1|RP1L1|RPE65|PRPF6|SNRNP200|PRPF3\t|KCNJ13|RIMS1|IMPG1|RRM2B|TUBB4B|NRL|RDH12|GUCY2D" > $work/retina-AD.tsv
awk -F"\t" '{print "\t" $1 "\t"}' $work/retina-all.tsv > $work/retina-all-genes.tsv
awk -F"\t" '{print "\t" $1 "\t"}' $work/retina-AD.tsv > $work/retina-AD-genes.tsv
rm $work/retina-all.tsv $work/retina-AD.tsv

# interval with eoxnic regions
grep -Ff $work/retina-all-genes.tsv $here/refs/hg19_ncbiRefSeq.tsv | grep -v "NM_001375654.1" | awk -F"\t" '{n=split($10,start,","); split($11,end,","); for(i=1; i<n; i=i+1) {if($3 !~ /_/) {print $3 "\t" start[i] "\t" end[i] "\t" $13} } }' | sort | uniq | grep -v "_" | sort -V -k1,1 -k2,2 > $work/intervals.bed
grep -Ff $work/retina-AD-genes.tsv $here/refs/hg19_ncbiRefSeq.tsv | grep -v "NM_001375654.1" | awk -F"\t" '{n=split($10,start,","); split($11,end,","); for(i=1; i<n; i=i+1) {if($3 !~ /_/) {print $3 "\t" start[i] "\t" end[i] "\t" $13} } }' | sort | uniq | grep -v "_" | sort -V -k1,1 -k2,2 > $work/intervals.AD.bed

# OFF-PEAK outputs
cat $here/unsolved-2020/04_CNVs-results/CNVs-all.tsv $here/unsolved-2021/04_CNVs-results/CNVs-all.tsv $here/unsolved-2022/04_CNVs-results/CNVs-all.tsv | grep -v "Inf" > $work/02.OFF-PEAK.merge.tsv
cat $here/unsolved-2020/04_CNVs-results/CNVs-targets-only.tsv $here/unsolved-2021/04_CNVs-results/CNVs-targets-only.tsv $here/unsolved-2022/04_CNVs-results/CNVs-targets-only.tsv | grep -v "Inf" > $work/02.OFF-PEAK-ON.merge.tsv
cat $here/unsolved-2020/04_CNVs-results/CNVs-all.tsv $here/unsolved-2021/04_CNVs-results/CNVs-all.tsv $here/unsolved-2022/04_CNVs-results/CNVs-all.tsv $here/unsolved-2020/04_CNVs-results/CNVs-targets-only.tsv $here/unsolved-2021/04_CNVs-results/CNVs-targets-only.tsv $here/unsolved-2022/04_CNVs-results/CNVs-targets-only.tsv | grep -v "Inf" > $work/02.OFF-PEAK.all.tsv
cat $here/unsolved-2020/04_CNVs-results/CNVs-all.HQ.tsv $here/unsolved-2021/04_CNVs-results/CNVs-all.HQ.tsv $here/unsolved-2022/04_CNVs-results/CNVs-all.HQ.tsv $here/unsolved-2020/04_CNVs-results/CNVs-targets-only.HQ.tsv $here/unsolved-2021/04_CNVs-results/CNVs-targets-only.HQ.tsv $here/unsolved-2022/04_CNVs-results/CNVs-targets-only.HQ.tsv | grep -v "Inf" > $work/02.OFF-PEAK-HQ.all.tsv

# exclusion list for CNV more than 4 times
cat $work/02.OFF-PEAK.merge.tsv  | awk -F"\t" '{print $2 "\t" $3 "\t" $4}' | awk '{A[$0]++}END{for (i in A){if(A[i]>4){print i}}}' | sort -V -k1,1 -k2,2 | uniq > $work/02.OFF-PEAK.exclusion.tsv
cat $work/02.OFF-PEAK-ON.merge.tsv | awk -F"\t" '{print $2 "\t" $3 "\t" $4}' | awk '{A[$0]++}END{for (i in A){if(A[i]>4){print i}}}' | sort -V -k1,1 -k2,2 | uniq > $work/02.OFF-PEAK-ON.exclusion.tsv
cat $here/unsolved-202*/SavvyCNV/unsolved-202*.SavvyCNVON.out.tsv  | awk -F"\t" '{print $1 "\t" $2 "\t" $3}' | awk '{A[$0]++}END{for (i in A){if(A[i]>4){print i}}}' | sort -V -k1,1 -k2,2 | uniq > $work/02.SavvyCNVON.exclusion.tsv
cat $here/unsolved-202*/SavvyCNV/unsolved-202*.SavvyCNVOFF.out.tsv  | awk -F"\t" '{print $1 "\t" $2 "\t" $3}' | awk '{A[$0]++}END{for (i in A){if(A[i]>4){print i}}}' | sort -V -k1,1 -k2,2 | uniq > $work/02.SavvyCNVOFF.exclusion.tsv
for tool in gatk ExomeDepth conifer cnmops codex2 controlFREEC cnvkit
do
	cat $here/unsolved-202*/$tool/unsolved-202*.$tool.out.tsv  | awk -F"\t" '{print $1 "\t" $2 "\t" $3}' | awk '{A[$0]++}END{for (i in A){if(A[i]>4){print i}}}' | sort -V -k1,1 -k2,2 | uniq > $work/02.$tool.exclusion.tsv
done

# all CNVs without excluded ones
grep -h -v -Ff $work/02.OFF-PEAK.exclusion.tsv $work/02.OFF-PEAK.all.tsv > $work/02.OFF-PEAK.tsv
grep -h -v -Ff $work/02.OFF-PEAK-ON.exclusion.tsv $work/02.OFF-PEAK-HQ.all.tsv > $work/02.OFF-PEAK-HQ.tsv
grep -h -v -Ff $work/02.SavvyCNVON.exclusion.tsv $here/unsolved-202*/SavvyCNV/unsolved-202*.SavvyCNVON.out.tsv  > $work/02.SavvyCNVON.tsv
grep -h -v -Ff $work/02.SavvyCNVOFF.exclusion.tsv $here/unsolved-202*/SavvyCNV/unsolved-202*.SavvyCNVOFF.out.tsv  > $work/02.SavvyCNVOFF.tsv
for tool in gatk ExomeDepth conifer cnmops codex2 controlFREEC cnvkit
do
	grep -h -v -Ff $work/02.$tool.exclusion.tsv $here/unsolved-202*/$tool/unsolved-202*.$tool.out.tsv  > $work/02.$tool.tsv
done

# number after gene and common removal
awk -F"\t" '{print $2 "\t" $3 "\t" $4 "\t" $1}' $work/02.OFF-PEAK.tsv | sort -V -k1,1 -k2,2 | tail -n+7 > $work/02.OFF-PEAK.f.tsv
rm $work/02.OFF-PEAK.f2.tsv
touch $work/02.OFF-PEAK.f2.tsv
cut -f 4 $work/02.OFF-PEAK.f.tsv | sort | uniq | while read C
do
     awk -v C=${C} '($4==C)' $work/02.OFF-PEAK.f.tsv | sort -t $'\t' -k1,1 -k2,2n | bedtools merge >> $work/02.OFF-PEAK.f2.tsv
done
bedtools intersect -wa -a $work/02.OFF-PEAK.f2.tsv -b $work/intervals.bed | sort | uniq | wc -l

awk -F"\t" '{print $2 "\t" $3 "\t" $4 "\t" $1}' $work/02.OFF-PEAK-HQ.tsv | sort -V -k1,1 -k2,2 | tail -n+7 > $work/02.OFF-PEAK-HQ.f.tsv
rm $work/02.OFF-PEAK-HQ.f2.tsv
touch $work/02.OFF-PEAK-HQ.f2.tsv
cut -f 4 $work/02.OFF-PEAK-HQ.f.tsv | sort | uniq | while read C
do
     awk -v C=${C} '($4==C)' $work/02.OFF-PEAK-HQ.f.tsv | sort -t $'\t' -k1,1 -k2,2n | bedtools merge >> $work/02.OFF-PEAK-HQ.f2.tsv
done
bedtools intersect -wa -a $work/02.OFF-PEAK-HQ.f2.tsv -b $work/intervals.bed | sort | uniq | wc -l

cat $work/02.SavvyCNVON.tsv $work/02.SavvyCNVOFF.tsv | cut -f1-4 | sort -V -k1,1 -k2,2 > $work/02.Savvy.f.tsv
rm $work/02.Savvy.f2.tsv
touch $work/02.Savvy.f2.tsv
cut -f 4 $work/02.Savvy.f.tsv | sort | uniq | while read C
do
     awk -v C=${C} '($4==C)' $work/02.Savvy.f.tsv | sort -t $'\t' -k1,1 -k2,2n | bedtools merge >> $work/02.Savvy.f2.tsv
done
bedtools intersect -wa -a $work/02.Savvy.f2.tsv -b $work/intervals.bed | sort | uniq | wc -l

for tool in gatk ExomeDepth conifer cnmops codex2 controlFREEC cnvkit
do
	echo $tool
	bedtools intersect -wa -a $work/02.$tool.tsv -b $work/intervals.bed | sort | uniq | wc -l
done

rm $work/01.all.tsv
touch $work/01.all.tsv

for pat in "${patients[@]}"
do
	echo $pat

	pat2=$(echo $pat | sed 's/-/\./g')

	grep -h -P "$pat|$pat2" $work/02.OFF-PEAK.tsv | awk -F"\t" '{print $2 "\t" $3 "\t" $4 "\t" $1}' | grep -v "Inf" > $work/$pat.OFF-PEAK.tsv
	grep -h -P "$pat|$pat2" $work/02.OFF-PEAK-HQ.tsv | awk -F"\t" '{print $2 "\t" $3 "\t" $4 "\t" $1}' | grep -v "Inf" > $work/$pat.OFF-PEAK-HQ.tsv
	grep -h -P "$pat|$pat2" $work/02.SavvyCNVON.tsv > $work/$pat.SavvyON.tsv
	grep -h -P "$pat|$pat2" $work/02.SavvyCNVOFF.tsv > $work/$pat.SavvyOFF.tsv
	grep -h -P "$pat|$pat2" $work/02.gatk.tsv | awk -F"\t" '{if($1!="chrX") print $0; if($1=="chrX"){if($6==0) print $0}}' > $work/$pat.gatk.tsv
	grep -h -P "$pat|$pat2" $work/02.ExomeDepth.tsv > $work/$pat.ExomeDepth.tsv
	grep -h -P "$pat|$pat2" $work/02.conifer.tsv > $work/$pat.conifer.tsv
	grep -h -P "$pat|$pat2" $work/02.cnmops.tsv > $work/$pat.cnmops.tsv
	grep -h -P "$pat|$pat2" $work/02.codex2.tsv > $work/$pat.codex2.tsv
	grep -h -P "$pat|$pat2" $work/02.controlFREEC.tsv | awk -F"\t" '{if($1!="chrX") print $0; if($1=="chrX"){if($6==0) print $0}}' > $work/$pat.controlFREEC.tsv
	grep -h -P "$pat|$pat2" $work/02.cnvkit.tsv > $work/$pat.cnvkit.tsv
	cat $work/$pat.OFF-PEAK.tsv $work/$pat.OFF-PEAK-HQ.tsv $work/$pat.SavvyON.tsv $work/$pat.SavvyOFF.tsv $work/$pat.gatk.tsv $work/$pat.ExomeDepth.tsv | cut -f1-4 | sort -V -k1,1 -k2,2 > $work/$pat.all.tsv
	bedtools merge -i $work/$pat.all.tsv > $work/$pat.all.merge.tsv

	# AR - comphet CNV + variant
	bedtools intersect -wa -wb -a $work/$pat.all.merge.tsv -b $work/intervals.bed | cut -f1-3,7 | sort -V -k1,1 -k2,2 | uniq > $work/$pat.all.merge.Eye.tsv
	cut -f4 $work/$pat.all.merge.Eye.tsv | awk -F"\t" '{print "\t" $1 "\t"}' | uniq > $work/$pat.genes.tsv
	grep -Ff $work/$pat.genes.tsv $here2/CeGat_*/$pat/W*/c*/*.anno.0.1.HQ.0.01.22_7_3_5.SpAI.Pimpact-strong.tsv > $work/$pat.variants.tsv
	cut -f2 $work/$pat.variants.tsv | sort | uniq > $work/$pat.genes1.tsv

	# AR - comphet CNV + CNV
	bedtools intersect -wa -wb -a $work/$pat.all.merge.tsv -b $work/intervals.bed | cut -f1-3,7 | sort -V -k1,1 -k2,2 | uniq | cut -f4 | uniq -d > $work/$pat.genes2.tsv

	# AR - hom CNV
	grep -h -P "$pat|$pat2" $work/02.OFF-PEAK.tsv | awk -F"\t" '{if($21<0.125 || ($21>1.625 && $21<2.5) ) print $2 "\t" $3 "\t" $4 "\t" $1}' | grep -v "Inf" > $work/$pat.OFF-PEAK.hom.tsv
	grep -h -P "$pat|$pat2" $work/02.OFF-PEAK-HQ.tsv | awk -F"\t" '{if($21<0.125 || ($21>1.625 && $21<2.5) ) print $2 "\t" $3 "\t" $4 "\t" $1}' | grep -v "Inf" > $work/$pat.OFF-PEAK-HQ.hom.tsv
	grep -h -P "$pat|$pat2" $work/02.SavvyCNVON.tsv | awk -F"\t" '{if($6<0.25 || ($6>3.25 && $6<5) ) print $0}' > $work/$pat.SavvyON.hom.tsv
	grep -h -P "$pat|$pat2" $work/02.SavvyCNVOFF.tsv | awk -F"\t" '{if($6<0.25 || ($6>3.25 && $6<5) ) print $0}' > $work/$pat.SavvyOFF.hom.tsv
	grep -h -P "$pat|$pat2" $work/02.gatk.tsv | awk -F"\t" '{if($6<0.25 || ($6>3.25 && $6<5) ) print $0}' > $work/$pat.gatk.hom.tsv
	grep -h -P "$pat|$pat2" $work/$pat.ExomeDepth.tsv | awk -F"\t" '{if($6<0.25 || ($6>3.25 && $6<5) ) print $0}' > $work/$pat.ExomeDepth.hom.tsv
	cat $work/$pat.OFF-PEAK.hom.tsv $work/$pat.OFF-PEAK-HQ.hom.tsv $work/$pat.SavvyON.hom.tsv $work/$pat.SavvyOFF.hom.tsv $work/$pat.gatk.hom.tsv $work/$pat.ExomeDepth.hom.tsv | cut -f1-4 | sort -V -k1,1 -k2,2 > $work/$pat.all.hom.tsv
	bedtools merge -i $work/$pat.all.hom.tsv > $work/$pat.all.hom.merge.tsv
	bedtools intersect -wa -wb -a $work/$pat.all.hom.merge.tsv -b $work/intervals.bed | cut -f1-3,7 | sort -V -k1,1 -k2,2 | uniq > $work/$pat.all.hom.merge.Eye.tsv
	cut -f4 $work/$pat.all.hom.merge.Eye.tsv | sort | uniq > $work/$pat.genes3.tsv

	# AD - CNV
	bedtools intersect -wa -wb -a $work/$pat.all.merge.tsv -b $work/intervals.AD.bed | cut -f1-3,7 | sort -V -k1,1 -k2,2 | uniq > $work/$pat.all.merge.Eye.AD.tsv
	cut -f4 $work/$pat.all.merge.Eye.AD.tsv > $work/$pat.genes4.tsv

	# all
	cat $work/$pat.genes1.tsv $work/$pat.genes2.tsv $work/$pat.genes3.tsv $work/$pat.genes4.tsv | sort | uniq | awk -F"\t" '{print "\t" $1 "\t"}' > $work/$pat.genes5.tsv
	grep -Ff $work/$pat.genes5.tsv $here/refs/hg19_ncbiRefSeq.tsv | grep -v "NM_001375654.1" | awk -F"\t" '{n=split($10,start,","); split($11,end,","); for(i=1; i<n; i=i+1) {if($3 !~ /_/) {print $3 "\t" start[i] "\t" end[i] "\t" $13} } }' | sort | uniq | grep -v "_" | sort -V -k1,1 -k2,2 > $work/$pat.intervals.tsv

	echo $pat >> $work/01.all.tsv
	grep -h -P "$pat|$pat2" $work/02.OFF-PEAK.tsv | awk -F"\t" '{print $2 "\t" $3 "\t" $4 "\t" $1 "\tOFF-PEAK" "\t" $22 "\t" $21*2}' | grep -v "Inf" > $work/$pat.temp.tsv
	bedtools intersect -wa -wb -a $work/$pat.temp.tsv -b $work/$pat.intervals.tsv | cut -f1,2,3,4,5,6,7,11 | sort | uniq >> $work/01.all.tsv

	grep -h -P "$pat|$pat2" $work/02.OFF-PEAK-HQ.tsv | awk -F"\t" '{print $2 "\t" $3 "\t" $4 "\t" $1 "\tOFF-PEAK-HQ" "\t" $22 "\t" $21*2}' | grep -v "Inf" > $work/$pat.temp.tsv
	bedtools intersect -wa -wb -a $work/$pat.temp.tsv -b $work/$pat.intervals.tsv | cut -f1,2,3,4,5,6,7,11 | sort | uniq >> $work/01.all.tsv

	grep -h -P "$pat|$pat2" $work/02.SavvyCNVON.tsv  | awk -F"\t" '{print $1 "\t" $2 "\t" $3 "\t" $4 "\tSavvyON" "\t" $5 "\t" $6}' > $work/$pat.temp.tsv
	bedtools intersect -wa -wb -a $work/$pat.temp.tsv -b $work/$pat.intervals.tsv | cut -f1,2,3,4,5,6,7,11 | sort | uniq >> $work/01.all.tsv

	grep -h -P "$pat|$pat2" $work/02.SavvyCNVOFF.tsv | awk -F"\t" '{print $1 "\t" $2 "\t" $3 "\t" $4 "\tSavvyOFF" "\t" $5 "\t" $6}' > $work/$pat.temp.tsv
	bedtools intersect -wa -wb -a $work/$pat.temp.tsv -b $work/$pat.intervals.tsv | cut -f1,2,3,4,5,6,7,11 | sort | uniq >> $work/01.all.tsv

	grep -h -P "$pat|$pat2" $work/02.gatk.tsv | awk -F"\t" '{print $1 "\t" $2 "\t" $3 "\t" $4 "\tgatk" "\t" $5 "\t" $6}' > $work/$pat.temp.tsv
	bedtools intersect -wa -wb -a $work/$pat.temp.tsv -b $work/$pat.intervals.tsv | cut -f1,2,3,4,5,6,7,11 | sort | uniq >> $work/01.all.tsv

	grep -h -P "$pat|$pat2" $work/$pat.ExomeDepth.tsv | awk -F"\t" '{print $1 "\t" $2 "\t" $3 "\t" $4 "\tExomeDepth" "\t" $5 "\t" $6}' > $work/$pat.temp.tsv
	bedtools intersect -wa -wb -a $work/$pat.temp.tsv -b $work/$pat.intervals.tsv | cut -f1,2,3,4,5,6,7,11 | sort | uniq >> $work/01.all.tsv

	grep -h -P "$pat|$pat2" $work/$pat.conifer.tsv | awk -F"\t" '{print $1 "\t" $2 "\t" $3 "\t" $4 "\tconifer" "\t" $5 "\t" $6}' > $work/$pat.temp.tsv
	bedtools intersect -wa -wb -a $work/$pat.temp.tsv -b $work/$pat.intervals.tsv | cut -f1,2,3,4,5,6,7,11 | sort | uniq >> $work/01.all.tsv

	grep -h -P "$pat|$pat2" $work/$pat.cnmops.tsv | awk -F"\t" '{print $1 "\t" $2 "\t" $3 "\t" $4 "\tcnmops" "\t" $5 "\t" $6}' > $work/$pat.temp.tsv
	bedtools intersect -wa -wb -a $work/$pat.temp.tsv -b $work/$pat.intervals.tsv | cut -f1,2,3,4,5,6,7,11 | sort | uniq >> $work/01.all.tsv

	grep -h -P "$pat|$pat2" $work/$pat.codex2.tsv | awk -F"\t" '{print $1 "\t" $2 "\t" $3 "\t" $4 "\tcodex2" "\t" $5 "\t" $6}' > $work/$pat.temp.tsv
	bedtools intersect -wa -wb -a $work/$pat.temp.tsv -b $work/$pat.intervals.tsv | cut -f1,2,3,4,5,6,7,11 | sort | uniq >> $work/01.all.tsv

	grep -h -P "$pat|$pat2" $work/$pat.controlFREEC.tsv | awk -F"\t" '{print $1 "\t" $2 "\t" $3 "\t" $4 "\tcontrolFREEC" "\t" $5 "\t" $6}' > $work/$pat.temp.tsv
	bedtools intersect -wa -wb -a $work/$pat.temp.tsv -b $work/$pat.intervals.tsv | cut -f1,2,3,4,5,6,7,11 | sort | uniq >> $work/01.all.tsv

	grep -h -P "$pat|$pat2" $work/$pat.cnvkit.tsv | awk -F"\t" '{print $1 "\t" $2 "\t" $3 "\t" $4 "\tcnvkit" "\t" $5 "\t" $6}' > $work/$pat.temp.tsv
	bedtools intersect -wa -wb -a $work/$pat.temp.tsv -b $work/$pat.intervals.tsv | cut -f1,2,3,4,5,6,7,11 | sort | uniq >> $work/01.all.tsv

	grep -Ff $work/$pat.genes.tsv $here2/CeGat_*/$pat/W*/c*/*.anno.0.1.HQ.0.01.22_7_3_5.SpAI.Pimpact-strong.tsv >> $work/01.all.tsv
	echo "" >> $work/01.all.tsv
	echo "" >> $work/01.all.tsv

	rm $work/$pat.*

done



