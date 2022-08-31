
here=/home/mquinodo/SYNO/scripts_NGS_analysis/OFF-PEAK
batch=ICR96
work=$here/$batch/00_comparison-targets
mkdir -p $work

cat $here/$batch/*/$batch.SavvyCNVON.out.tsv > $here/$batch/SavvyCNV/$batch.SavvyCNV.out.tsv 
cp $here/$batch/*/$batch.*.out.tsv $work

# OFF-PEAK
tail -n+2 $work/CNVs-targets-only.tsv | awk -F"\t" '{if($5!="normal") print $2 "\t" $3 "\t" $4 "\t" $1 "\t" $21 "\t" $21*2}' | sed 's/X1/1/g' | sort -V -k1,1 -k2,2 > $work/$batch.OFF-PEAK.out.tsv
tail -n+2 $work/CNVs-targets-only.HQ.tsv | awk -F"\t" '{if($5!="normal") print $2 "\t" $3 "\t" $4 "\t" $1 "\t" $21 "\t" $21*2}' | sed 's/X1/1/g' | sort -V -k1,1 -k2,2 > $work/$batch.OFF-PEAK-HQ.out.tsv
folder=/home/mquinodo/$batch.temp
mkdir -p $folder
cp $work/* $folder

echo -e "batch\ttool\tTP\tFN\tFP\tsens\tspec" > $folder/results.tsv

# true set
grep -v Chromosome $folder/MLPA.tsv | grep -v "Normal" | awk -F"\t" '{print "chr" $8 "\t" $9 "\t" $10 "\t" $1 "\t" $6 "\t" $3}' | sort -V -k1,1 -k2,2 > $folder/ICR96-true.tsv

# tools
for tool in OFF-PEAK OFF-PEAK-HQ cnmops codex2 controlFREEC ExomeDepth SavvyCNV gatk cnvkit  # conifer
do
	#echo -e "ID\tTP\tFN\tFP\tTN\tsens\tspec\tpres\tnTP\tnFN\tnFP" > $folder/$tool.res.tsv
	echo -e "ID\tnTP\tnFN\tnFP" > $folder/$tool.res.tsv

	for pat in 17296 17297 17298 17299 17300 17301 17302 17303 17304 17305 17306 17307 17308 17309 17310 17311 17312 17313 17314 17315 17316 17317 17318 17319 17320 17321 17322 17323 17324 17325 17326 17327 17328 17329 17330 17331 17332 17333 17334 17335 17336 17337 17338 17339 17340 17341 17342 17343 17356 17357 17358 17359 17360 17361 17362 17363 17364 17365 17366 17367 17368 17369 17370 17371 17372 17373 17374 17375 17376 17377 17378 17379 17380 17381 17382 17383 17384 17385 17386 17387 17388 17389 17390 17391 17392 17393 17394 17395 17396 17397 17398 17399 17400 17401 17402 17403
	do
	echo $pat
	grep Deletion $folder/ICR96-true.tsv | grep -P "\t$pat\t" > $folder/test.$pat.del.tsv
	grep Duplication $folder/ICR96-true.tsv | grep -P "\t$pat\t" > $folder/test.$pat.dup.tsv
	awk -F"\t" '{if($6<2) print $0}' $folder/$batch.$tool.out.tsv | grep -P "\t$pat\t" > $folder/test.$pat.del-test.tsv
	awk -F"\t" '{if($6>2) print $0}' $folder/$batch.$tool.out.tsv | grep -P "\t$pat\t" > $folder/test.$pat.dup-test.tsv

	# take true CNV that are at least partially in tool
	bedtools intersect -wa -a $folder/test.$pat.del.tsv -b $folder/test.$pat.del-test.tsv | sort | uniq > $folder/test.$pat.del.inter.tsv
	bedtools intersect -wa -a $folder/test.$pat.dup.tsv -b $folder/test.$pat.dup-test.tsv | sort | uniq > $folder/test.$pat.dup.inter.tsv
	# take true CNV that are NOT in tool
	bedtools intersect -v -a $folder/test.$pat.del.tsv -b $folder/test.$pat.del-test.tsv > $folder/test.$pat.del.true-only.tsv
	bedtools intersect -v -a $folder/test.$pat.dup.tsv -b $folder/test.$pat.dup-test.tsv > $folder/test.$pat.dup.true-only.tsv
	# take tool CNV that are NOT TRUE true
	bedtools intersect -v -b $folder/test.$pat.del.tsv -a $folder/test.$pat.del-test.tsv > $folder/test.$pat.del.tool-only.tsv
	bedtools intersect -v -b $folder/test.$pat.dup.tsv -a $folder/test.$pat.dup-test.tsv > $folder/test.$pat.dup.tool-only.tsv
	
	nTP=$(cat $folder/test.$pat.del.inter.tsv $folder/test.$pat.dup.inter.tsv | wc -l)
	nFN=$(cat $folder/test.$pat.del.true-only.tsv $folder/test.$pat.dup.true-only.tsv | wc -l)
	nFP=$(cat $folder/test.$pat.del.tool-only.tsv $folder/test.$pat.dup.tool-only.tsv | wc -l)
	if [ -z "$nTP" ];then nTP=0;fi
	if [ -z "$nFN" ];then nFN=0;fi
	if [ -z "$nFP" ];then nFP=0;fi

	echo -e "$pat\t$nTP\t$nFN\t$nFP" >> $folder/$tool.res.tsv
	rm $folder/test.$pat.*
	done

	TP=$(cut -f2 $folder/$tool.res.tsv | tail -n+2 | awk -F'\t' '{sum+=$1;} END{print sum;}')
	FP=$(cut -f3 $folder/$tool.res.tsv | tail -n+2 | awk -F'\t' '{sum+=$1;} END{print sum;}')
	FN=$(cut -f4 $folder/$tool.res.tsv | tail -n+2 | awk -F'\t' '{sum+=$1;} END{print sum;}')
	T=$(($TP+$FN))
	P=$(($TP+$FP))
	sens=$(calc $TP/$T | sed 's/~//g' | sed -e 's/[ \t]*//')
	spec=$(calc $TP/$P | sed 's/~//g' | sed -e 's/[ \t]*//')
	echo -e "$batch\t$tool\t$TP\t$FN\t$FP\t$sens\t$spec" >> $folder/results.tsv
done

cp $folder/* $work
rm -rf $folder

