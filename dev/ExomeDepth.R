
args = commandArgs(trailingOnly=TRUE)

if (length(args)<4) {
  stop("At least 4 arguments must be supplied (patient's name, input file, output file).n", call.=FALSE)
} else {
	here=args[1]
	patients=strsplit(args[2],",")[[1]]
	type=args[3]
	out=args[4]
}

# here="/home/mquinodo/ICR96"
# patients=strsplit("17296,17297,17298,17299,17300,17301,17302,17303,17304,17305,17306,17307,17308,17309,17310,17311,17312,17313,17314,17315,17316,17317,17318,17319,17320,17321,17322,17323,17324,17325,17326,17327,17328,17329,17330,17331,17332,17333,17334,17335,17336,17337,17338,17339,17340,17341,17342,17343,17356,17357,17358,17359,17360,17361,17362,17363,17364,17365,17366,17367,17368,17369,17370,17371,17372,17373,17374,17375,17376,17377,17378,17379,17380,17381,17382,17383,17384,17385,17386,17387,17388,17389,17390,17391,17392,17393,17394,17395,17396,17397,17398,17399,17400,17401,17402,17403",",")[[1]]
# type=".dup.BQSR.bam"
# out="/home/mquinodo/SYNO/scripts_NGS_analysis/OFF-PEAK-train3/ICR96/ExomeDepth"

library(ExomeDepth)
data(exons.hg19)

# list of BAM files
bams<-vector(mode="character",length=length(patients))
for (i in 1:length(patients)){
	bams[i]=paste(here,"/",patients[i],type,sep="")
}

# parsing BAM files
ExomeCount <- getBamCounts(bed.frame = exons.hg19,
	bam.files = bams,
	include.chr = TRUE)

save(ExomeCount,file=paste(out,"/temp.RData",sep=""))
	
ExomeCount.dafr <- as(ExomeCount[, colnames(ExomeCount)], 'data.frame')
ExomeCount.dafr$chromosome <- gsub(as.character(ExomeCount.dafr$chromosome),
pattern = 'chr',
replacement = '') 

save(ExomeCount,ExomeCount.dafr,file=paste(out,"/autosomes.RData",sep=""))

######

data(exons.hg19.X)

# parsing BAM files
ExomeCount <- getBamCounts(bed.frame = exons.hg19.X,
	bam.files = bams,
	include.chr = TRUE)
ExomeCount.dafr <- as(ExomeCount[, colnames(ExomeCount)], 'data.frame')
ExomeCount.dafr$chromosome <- gsub(as.character(ExomeCount.dafr$chromosome),
pattern = 'chr',
replacement = '') 

save(ExomeCount,ExomeCount.dafr,file=paste(out,"/chrX.RData",sep=""))

######


load(file=paste(out,"/autosomes.RData",sep=""))

for (j in 1:length(patients)){
	print(patients[j])
	pat=gsub("-", ".", patients[j])
	file=paste(pat,type,sep="")
	# selecting reference set
	colnames(ExomeCount)=gsub("^X", "", colnames(ExomeCount))
	my.test <- ExomeCount[,file]
	my.ref.samples<-vector(mode="character",length=length(patients)-1)
	names=patients[-j]
	for (i in 1:length(names)){
		my.ref.samples[i]=paste(names[i],type,sep="")
	}
	my.ref.samples=gsub("-", ".", my.ref.samples)
	colnames(ExomeCount.dafr)=gsub("^X", "", colnames(ExomeCount.dafr))
	my.reference.set <- as.matrix(ExomeCount.dafr[, my.ref.samples])
	my.choice <- select.reference.set (test.counts = my.test,
		reference.counts = my.reference.set,
		bin.length = (ExomeCount.dafr$end - ExomeCount.dafr$start)/1000,
		n.bins.reduced = 10000)

	# making reference set
	my.matrix <- as.matrix( ExomeCount.dafr[, my.choice$reference.choice, drop = FALSE])

	my.reference.selected <- apply(X = my.matrix,
	MAR = 1,
	FUN = sum)

	# CNV calling
	all.exons <- new('ExomeDepth',
		test = my.test,
		reference = my.reference.selected,
		formula = 'cbind(test, reference) ~ 1')

	all.exons <- CallCNVs(x = all.exons,
		transition.probability = 10^-4,
		chromosome = ExomeCount.dafr$chromosome,
		start = ExomeCount.dafr$start,
		end = ExomeCount.dafr$end,
		name = exons.hg19$name)

	all.exons@CNV.calls$chromosome <- gsub(as.character(all.exons@CNV.calls$chromosome),
	pattern = 'chr',
	replacement = '') 

	# annotation
	data(Conrad.hg19)

	# all.exons <- AnnotateExtra(x = all.exons,
	# 	reference.annotation = Conrad.hg19.common.CNVs,
	# 	min.overlap = 0.1,
	# 	column.name = 'Conrad.hg19')

	exons.hg19.GRanges <- GenomicRanges::GRanges(seqnames = exons.hg19$chromosome,
		IRanges::IRanges(start=exons.hg19$start,end=exons.hg19$end),
		names = exons.hg19$name)
	
	# all.exons <- AnnotateExtra(x = all.exons,
	# 	reference.annotation = exons.hg19.GRanges,
	# 	min.overlap = 0.0001,
	# 	column.name = 'exons.hg19')

	# output
	output.file <- paste(out,"/exome_calls_",pat,".csv",sep="")
	write.csv(file = output.file,
		x = all.exons@CNV.calls,
		row.names = FALSE)
}

data(exons.hg19.X)
load(file=paste(out,"/chrX.RData",sep=""))

for (j in 1:length(patients)){
	pat=gsub("-", ".", patients[j])
	file=paste(pat,type,sep="")
	# selecting reference set
	colnames(ExomeCount)=gsub("^X", "", colnames(ExomeCount))
	my.test <- ExomeCount[,file]
	my.ref.samples<-vector(mode="character",length=length(patients)-1)
	names=patients[-j]
	for (i in 1:length(names)){
		my.ref.samples[i]=paste(names[i],type,sep="")
	}
	my.ref.samples=gsub("-", ".", my.ref.samples)
	colnames(ExomeCount.dafr)=gsub("^X", "", colnames(ExomeCount.dafr))
	
	my.reference.set <- as.matrix(ExomeCount.dafr[, my.ref.samples])
	my.choice <- select.reference.set (test.counts = my.test,
		reference.counts = my.reference.set,
		bin.length = (ExomeCount.dafr$end - ExomeCount.dafr$start)/1000,
		n.bins.reduced = 10000)

	# making reference set
	my.matrix <- as.matrix( ExomeCount.dafr[, my.choice$reference.choice, drop = FALSE])
		my.reference.selected <- apply(X = my.matrix,
		MAR = 1,
		FUN = sum)

	# CNV calling
	all.exons <- new('ExomeDepth',
		test = my.test,
		reference = my.reference.selected,
		formula = 'cbind(test, reference) ~ 1')

	all.exons <- CallCNVs(x = all.exons,
		transition.probability = 10^-4,
		chromosome = ExomeCount.dafr$chromosome,
		start = ExomeCount.dafr$start,
		end = ExomeCount.dafr$end,
		name = exons.hg19.X$name)

	all.exons@CNV.calls$chromosome <- gsub(as.character(all.exons@CNV.calls$chromosome),
	pattern = 'chr',
	replacement = '') 

	# annotation
	data(Conrad.hg19)

	# all.exons <- AnnotateExtra(x = all.exons,
	# 	reference.annotation = Conrad.hg19.common.CNVs,
	# 	min.overlap = 0.1,
	# 	column.name = 'Conrad.hg19')

	exons.hg19.GRanges <- GenomicRanges::GRanges(seqnames = exons.hg19.X$chromosome,
		IRanges::IRanges(start=exons.hg19.X$start,end=exons.hg19.X$end),
		names = exons.hg19.X$name)
	# all.exons <- AnnotateExtra(x = all.exons,
	# 	reference.annotation = exons.hg19.GRanges,
	# 	min.overlap = 0.0001,
	# 	column.name = 'exons.hg19')

	# output
	output.file <- paste(out,"/exome_calls_chrX_",pat,".csv",sep="")
	write.csv(file = output.file,
		x = all.exons@CNV.calls,
		row.names = FALSE)
}






