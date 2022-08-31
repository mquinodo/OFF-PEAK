
args = commandArgs(trailingOnly=TRUE)

if (length(args)<4) {
  stop("At least 4 arguments must be supplied (patient's name, input file, output file).n", call.=FALSE)
} else {
	here=args[1]
	patients=strsplit(args[2],",")[[1]]
	type=args[3]
	out=args[4]
  bedFile=args[5]
}

# here="/home/mquinodo/CeGat_2021-11"
# patients=strsplit("CHbasl0088,CHbasl0249,CHbasl0257,CHbasl0258,CHbasl0259,CHbasl0263,CHbasl0267,CHbasl0268,CHbasl0269,CHbasl0271,CHbasl0272,CHbasl0273,CHbasl0276,CHbasl0278,CHbasl0279,CHbasl0280,CHbasl0281,CHbasl0282,CHbasl0285,CHbasl0288,CHbasl0290,CHbasl0292,CHbasl0294,CHbasl0297,CHbasl0298,CHbasl0299,CHbasl0300,CHbasl0302,CHbasl0303,CHbasl0306,CHbasl0308,CHbasl0309,CHbasl0310,CHlaus0262,CHlaus0264,CHlaus0265,CHlaus0266,CHlaus0286,CHlaus0287,CHlaus0295,CHlaus0296,CHlaus0301,CHlaus0305",",")[[1]]
# type=".dup.BQSR.bam"
# out="/home/mquinodo/SYNO/scripts_NGS_analysis/OFF-PEAK-train3/CeGat_2021-11/codex2"

# here="/home/mquinodo/ICR96"
# patients=strsplit("17296,17297,17298,17299,17300,17301,17302,17303,17304,17305,17306,17307,17308,17309,17310,17311,17312,17313,17314,17315,17316,17317,17318,17319,17320,17321,17322,17323,17324,17325,17326,17327,17328,17329,17330,17331,17332,17333,17334,17335,17336,17337,17338,17339,17340,17341,17342,17343,17356,17357,17358,17359,17360,17361,17362,17363,17364,17365,17366,17367,17368,17369,17370,17371,17372,17373,17374,17375,17376,17377,17378,17379,17380,17381,17382,17383,17384,17385,17386,17387,17388,17389,17390,17391,17392,17393,17394,17395,17396,17397,17398,17399,17400,17401,17402,17403",",")[[1]]
# type=".dup.BQSR.bam"
# out="/home/mquinodo/SYNO/scripts_NGS_analysis/OFF-PEAK-train3/ICR96/codex2"
# bedFile="/home/mquinodo/SYNO/WES/ICR96/03_info/ICR96.noXY.bed"


library(CODEX2)
library(WES.1KG.WUGSC) # Load Toy data from the 1000 Genomes Project.

bamFile<-vector(mode="character",length=length(patients))
for (i in 1:length(patients)){
	bamFile[i]=paste(patients[i],type,sep="")
}
bamdir<-vector(mode="character",length=length(patients))
for (i in 1:length(patients)){
	bamdir[i]=paste(here,"/",patients[i],type,sep="")
}
sampname <- patients

bambedObj <- getbambed(bamdir = bamdir, bedFile = bedFile, 
                       sampname = sampname, projectname = "CODEX2_demo")

bamdir <- bambedObj$bamdir
sampname <- bambedObj$sampname
ref <- bambedObj$ref
projectname <- bambedObj$projectname

genome = BSgenome.Hsapiens.UCSC.hg19
gc <- getgc(ref, genome = genome)
mapp <- getmapp(ref, genome = genome)
values(ref) <- cbind(values(ref), DataFrame(gc, mapp)) 


coverageObj <- getcoverage(bambedObj, mapqthres = 20)
Y <- coverageObj$Y
write.csv(Y, file = paste(out,'/codex2_coverage.csv', sep=''), quote = FALSE)
head(Y[,1:5])


# QC
qcObj <- qc(Y, sampname, ref, cov_thresh = c(20, 4000),
            length_thresh = c(20, 2000), mapp_thresh = 0.9,
            gc_thresh = c(20, 80))
Y_qc <- qcObj$Y_qc; sampname_qc <- qcObj$sampname_qc
ref_qc <- qcObj$ref_qc; qcmat <- qcObj$qcmat; gc_qc <- ref_qc$gc
write.table(qcmat, file = paste(out, '/codex2_qcmat', '.txt', sep=''),
            sep = '\t', quote = FALSE, row.names = FALSE)


Y.nonzero <- Y_qc[apply(Y_qc, 1, function(x){!any(x==0)}),]
pseudo.sample <- apply(Y.nonzero,1,function(x){exp(1/length(x)*sum(log(x)))})
N <- apply(apply(Y.nonzero, 2, function(x){x/pseudo.sample}), 2, median)


normObj.null <- normalize_null(Y_qc = Y_qc,
                               gc_qc = gc_qc,
                               K = 1:5, N = N)
Yhat.null <- normObj.null$Yhat
AIC.null <- normObj.null$AIC
BIC.null <- normObj.null$BIC
RSS.null <- normObj.null$RSS

choiceofK(AIC.null, BIC.null, RSS.null, K = 1:5 , filename = paste(out,"/codex2_null_choiceofK.pdf",sep=""))

chr=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22")
finalcall.CBS <- segmentCBS(Y_qc,  
                            Yhat.null, optK = which.max(BIC.null),
                            K = 1:5,
                            sampname_qc = colnames(Y_qc),
                            ref_qc = ranges(ref_qc),
                            chr = chr, lmax = 400, mode = "integer")

write.table(as.data.frame(finalcall.CBS), file = paste(out,"/output.tsv",sep=""),sep="\t",quote=F)

filter1 <- finalcall.CBS$length_kb<=200
filter2 <- finalcall.CBS$length_kb/(finalcall.CBS$ed_exon-finalcall.CBS$st_exon+1)<50
finalcall.CBS.filter <- finalcall.CBS[filter1 & filter2, ]

filter3 <- finalcall.CBS.filter$lratio>40
filter4 <- (finalcall.CBS.filter$ed_exon-finalcall.CBS.filter$st_exon)>1
finalcall.CBS.filter=finalcall.CBS.filter[filter3|filter4,]

write.table(as.data.frame(finalcall.CBS.filter), file = paste(out,"/output.filter.tsv",sep=""),sep="\t",quote=F)

