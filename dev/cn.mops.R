
args = commandArgs(trailingOnly=TRUE)

if (length(args)<4) {
  stop("At least 4 arguments must be supplied (patient's name, input file, output file).n", call.=FALSE)
} else {
	here=args[1]
	patients=strsplit(args[2],",")[[1]]
	type=args[3]
	out=args[4]
}

# here="/home/mquinodo/SYNO/WES/EXOMES/CeGat_2020-03"
# patients=strsplit("CHlaus0001,CHlaus0002,CHlaus0003,CHlaus0004,CHlaus0005,CHlaus0006,CHlaus0007,CHlaus0008,CHlaus0009,CHbern0010",",")[[1]]
# type=".dup.BQSR.bam"
# out="/home/mquinodo/SYNO/scripts_NGS_analysis/OFF-PEAK-02/scripts/CeGat_2020-03"

library(cn.mops)

BAMFiles<-vector(mode="character",length=length(patients))
for (i in 1:length(patients)){
	#BAMFiles[i]=paste(here,"/",patients[i],"/bam/",patients[i],type,sep="")
	BAMFiles[i]=paste(here,"/",patients[i],type,sep="")
}
segments <- read.table("/home/mquinodo/SYNO/NGS/Twist/Twist_Exome_Target_hg19.cano.bed",sep="\t",as.is=TRUE)
gr <- GRanges(segments[,1],IRanges(segments[,2],segments[,3]))
X <- getSegmentReadCountsFromBAM(BAMFiles,GR=gr)
X2=as.data.frame(X)
s=apply(X2[,6:dim(X2)[2]],1,sum)
X3=X[-which(s==0),]
resCNMOPS <- exomecn.mops(X3)
resCNMOPS <- calcIntegerCopyNumbers(resCNMOPS)


segm <- as.data.frame(segmentation(resCNMOPS))
CNVs <- as.data.frame(cnvs(resCNMOPS))
CNVRegions <- as.data.frame(cnvr(resCNMOPS))

write.csv(segm,file=paste(out,"/cnmops.segmentation.csv",sep=""))
write.csv(CNVs,file=paste(out,"/cnmops.cnvs.csv",sep=""))
write.csv(CNVRegions,file=paste(out,"/cnmops.cnvr.csv",sep=""))



