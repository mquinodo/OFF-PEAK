
#!/usr/bin/env Rscript

library(optparse)

option_list = list(
  make_option(c("--output"), type="character", default="NA", 
              help="output folder", metavar="character"),

  make_option(c("--data"), type="character", default="NA", 
              help="file with number of reads per target", metavar="character"),

  make_option(c("--mincor"), type="numeric", default="0.9", 
              help="minimum correlation between samples for CNV detection", metavar="numeric"),

  make_option(c("--minsignal"), type="numeric", default="2500", 
              help="minimum signal per target", metavar="numeric"),

  make_option(c("--maxvar"), type="numeric", default="-0.2", 
              help="maximum variance per target", metavar="numeric"),

  make_option(c("--leaveoneout"), type="numeric", default="1", 
              help="leave-one-out PCA method, yes=1, no=0", metavar="numeric"),

  make_option(c("--downsample"), type="numeric", default="20000", 
              help="number of targets to downsample for noise reduction", metavar="numeric"),

  make_option(c("--nbFake"), type="numeric", default="500", 
              help="number of fake CNVs for noise reduction", metavar="numeric"),

  make_option(c("--stopPC"), type="numeric", default="0.0001", 
              help="stopping criteria for PC removal", metavar="numeric"),

  make_option(c("--minZ"), type="numeric", default="4", 
              help="minimum abs(Z-score) for single-exon CNVs", metavar="numeric"),

  make_option(c("--minOfftarget"), type="numeric", default="1000", 
              help="minimum size for off-targets with no exons", metavar="numeric"),

  make_option(c("--databasefile"), type="character", default="NA", 
              help="RData file for CNV annotation", metavar="character"),

  make_option(c("--chromosome-plots"), action="store_true", default=FALSE, 
              help="Will do coverage plots per chromosome"),

  make_option(c("--genome-plots"), action="store_true", default=FALSE, 
              help="Will do genome-wide coverage plots"),

  make_option(c("--nb-plots"), type="numeric", default="10", 
              help="Number of plots for CNVs per individual")
); 
 
opt_parser = OptionParser(option_list=option_list);
args = parse_args(opt_parser);

folder=as.character(args[1])
data=as.character(args[2])
mincor=as.numeric(args[3])
minsignal=as.numeric(args[4])
maxvar=as.numeric(args[5])
leaveoneout=as.numeric(args[6])
downsample=as.numeric(args[7])
nbFake=as.numeric(args[8])
stopPC=as.numeric(args[9])
minZ=as.numeric(args[10])
minOfftarget=as.numeric(args[11])
databasefile=as.character(args[12])
chromoPlots=as.logical(args[13])
genomePlots=as.logical(args[14])
nbPlots=as.numeric(args[15])

if(data=="NA"){
  stop("You need to include the files containing the readcount per target with the --data option. Exit.")
}
if(folder=="NA"){
  stop("You need to provide an output directory with the --output option. Exit.")
}
if(databasefile=="NA"){
  stop("You need to include the database file with the --databasefile option. Exit.")
}


# folder="/home/mquinodo/SYNO/WES/EXOMES/Carmen_Igor/12_OFF-PEAK"
# data="/home/mquinodo/SYNO/WES/EXOMES/Carmen_Igor/12_OFF-PEAK/ALL.target.tsv"
# mincor=0.9
# minsignal=2500
# maxvar=-0.2
# leaveoneout=1
# downsample=20000
# nbFake=500
# stopPC=0.0001
# minZ=4
# minOfftarget=1000
# databasefile="/home/mquinodo/SYNO/scripts_NGS_analysis/OFF-PEAK-train5/refs/data-hg19.RData"
# genomePlots=TRUE
# chromoPlots=TRUE
# nbPlots=10



# loading libraries
library(gplots) # for heatmap.2 function
library(ExomeDepth) # for C_hmm function
library(pROC) # for roc function
library(caTools) # for runmean function

# load files for CNV annotation
load(databasefile)

# reading input files with readcounts
dataALL=read.table(file=data,header=T)

# removing off-targets smaller than minOfftarget and without exons inside
taken=which(grepl("NM",dataALL[,4])==T | as.numeric(dataALL[,3])-as.numeric(dataALL[,2])>minOfftarget)
dataALL=dataALL[taken,]

# separating GC content from other data
GC=dataALL[,5]
dataALL=dataALL[,-5]

# creating output directories
dir.create(folder, showWarnings = FALSE)
dir.create(paste(folder,"/01_general-stats",sep=""), showWarnings = FALSE)
dir.create(paste(folder,"/02_BED-files",sep=""), showWarnings = FALSE)
dir.create(paste(folder,"/03_Samples-info",sep=""), showWarnings = FALSE)
dir.create(paste(folder,"/04_CNVs-results",sep=""), showWarnings = FALSE)
dir.create(paste(folder,"/05_RData-files",sep=""), showWarnings = FALSE)
dir.create(paste(folder,"/06_PC-plots",sep=""), showWarnings = FALSE)

# saving parameters used
n=c("folder","data","mincor","minsignal","maxvar","leaveoneout","downsample","nbFake","stopPC","minZ","minOfftarget","databasefile")
out=cbind(n,args[1:12])
write.table(out,file=paste(folder,"/01_general-stats/log_parameters.tsv",sep=""),quote=F,sep="\t",row.names=F,col.names=F)

# important functions

# graphical representation
plotcnv <- function(chr,begin,end,ID,data,pdf,side,offtar,UseCano) {

  load(data)

  all=allPLOT

  if(offtar==F){all=all[which(grepl("Off-target",all[,4])==F),]}

  pdf(file=pdf,width=10,height=12)
  par(mfrow=c(1,1))

  all=all[which(all[,15]=="Yes"),]

  {

    # taking regions of interest
    sel=which(all[,1]==chr & as.numeric(all[,2])>=begin & as.numeric(all[,3])<=end)
    sel=max(sel[1]-side,0):min(dim(all)[1],(sel[length(sel)]+side))

    d=all[sel,]

    if(length(which(d[,12]==0))>0) {d=d[-which(d[,12]==0),]}

    pos=unique(d[,c(2,3,4)])
    pos=cbind(pos,pos,pos)
    pos[,4]=0
    pos[,5]=0
    pos[,6]=0
    pos[,7]=0
    pos[,8]=0
    pos[,9]=0
    pos=as.data.frame(pos)
    for (i in c(1,2,4:8)){
      pos[,i]=as.numeric(pos[,i])
    }
    for(i in 1:dim(pos)[1]){
      # corrected counts
      pos[i,4]=as.numeric(d[which(as.numeric(d[,2])==pos[i,1]),11]) # ID
      pos[i,5]=as.numeric(d[which(as.numeric(d[,2])==pos[i,1]),12]) # mean
      pos[i,6]=as.numeric(d[which(as.numeric(d[,2])==pos[i,1]),13]) # SD
      pos[i,7]=min(pos[i,5]-pos[i,6],pos[i,4])/pos[i,5]
      pos[i,8]=max(pos[i,5]+pos[i,6],pos[i,4])/pos[i,5]
      pos[i,9]="Yes"
    }

    ma=max(1.5*1.05,(pos[,4]/pos[,5])*1.5)
    mi=(-0.05)
    plot((pos[,4]/pos[,5])+1000,xlim=c(0-length(pos[,4])*0.24,length(pos[,4])*1.02),xaxs="i",yaxs="i",ylim=c(mi-(ma-mi)*0.5,ma),main=paste("Sample: ",ID,"  /  Position: ",chr,":",begin,"-",end,"\nHigh quality targets",sep=""),xlab="",ylab="",yaxt='n',xaxt='n',cex.main=1.6)
    
    # red box
    d1=which(pos[,1]>=begin & pos[,2]<=end)
    #polygon(c(min(d1)-0.45,max(d1)+0.45,min(d1)-0.45,max(d1)+0.45),c(-10,-10,1000,1000),col=rgb(1,0.9,1),border=NA)
    rect(min(d1)-0.45,-10,max(d1)+0.45,5,col=rgb(1,0.9,1),border=NA)
    abline(h=mi-(ma-mi)*0.5)

    # # boxes for dicarded
    # for (i in 1:dim(pos)[1]){
    #   if(pos[i,9]=="No"){
    #     polygon(c(i-0.5,i+0.5,i-0.5,i+0.5),c(-10,-10,1000,1000),col="lightgreen",border=NA)
    #   }
    # }
    
    # Y axis
    ya=seq(0,10,0.5)
    ya=ya[which(ya>mi & ya<ma)]
    axis(2,at=ya)

    # sd grey shape
    l=1:dim(pos)[1]
    low=1-2*(pos[,6]/pos[,5])
    high=1+2*(pos[,6]/pos[,5])
    polygon(c(l,rev(l)),c(low,rev(high)),col="gray85",border=NA)
    #polygon(c(-1000,1000,-1000,1000),c(-10,-10,mi,mi),col="white",border=NA)
    rect(-1000,-10,1000,mi,col="white",border=NA)
    d1=which(pos[,1]>=begin & pos[,2]<=end)
    #polygon(c(min(d1)-0.45,max(d1)+0.45,min(d1)-0.45,max(d1)+0.45),c(-1000,-1000,mi,mi),col=rgb(1,0.9,1),border=NA)
    rect(min(d1)-0.45,-5,max(d1)+0.45,mi,col=rgb(1,0.9,1),border=NA)
    
    abline(v=0-length(pos[,4])*0.24)
    abline(v=length(pos[,4])*1.02)
    abline(h=mi-(ma-mi)*0.5)
    abline(h=1)
    
    for (i in 1:length(ya)){
      if(ya[i]!=1){abline(h=ya[i],lty=2,col=1)}
    }
    
    lines(pos[,4]/pos[,5], lty = 2)
    abline(h=mi)
    tar=pos[which(grepl("Off-target",pos[,3])==F & pos[,9]=="Yes"),]
    anti=pos[which(grepl("Off-target",pos[,3])==T & pos[,9]=="Yes"),]
    tarD=pos[which(grepl("Off-target",pos[,3])==F & pos[,9]=="No"),]
    antiD=pos[which(grepl("Off-target",pos[,3])==T & pos[,9]=="No"),]

    points(which(grepl("Off-target",pos[,3])==F & pos[,9]=="Yes"),tar[,4]/tar[,5],pch=15,col="green")
    points(which(grepl("Off-target",pos[,3])==T & pos[,9]=="Yes"),anti[,4]/anti[,5],pch=17,col="blue")
    points(which(grepl("Off-target",pos[,3])==F & pos[,9]=="No"),tarD[,4]/tarD[,5],pch=15,col="red")
    points(which(grepl("Off-target",pos[,3])==T & pos[,9]=="No"),antiD[,4]/antiD[,5],pch=17,col="orange")

    legend(0-length(pos[,4])*0.23,ma/2,yjust=0.5,cex=0.75,legend=c("On-targets","Off-targets","Selected sample","2*SD controls","CNV region"),col=c("green","blue","black","darkgrey",rgb(1,0.9,1),"white"),pch=c(15,17,NA,15,15,NA),pt.cex=c(1,1,NA,2,2,NA),lty=c(NA,NA,2,NA,NA,NA),bg="white",y.intersp=1.5)
    
    pos2=pos[which(grepl("NM",pos[,3])==T),]
    
    if(dim(pos2)[1]>0){
      
      mtext('Observed / Expected read ratio', side=2, line=2.5,at=mi+(ma-mi)/2,cex=1.3)
      mtext('Transcripts', side=2, line=2.5,at=mi-(ma-mi)*0.15,cex=1.3)
      mtext('and exons', side=2, line=1,at=mi-(ma-mi)*0.15,cex=1.3)
      mtext("* indicates canonical isoforms", side=2, line=0.25,at=mi-(ma-mi)*0.15,cex=0.5)
      mtext('ClinVar and', side=2, line=2.5,at=mi-(ma-mi)*0.4,cex=1.3)
      mtext('gnomAD CNVs', side=2, line=1,at=mi-(ma-mi)*0.4,cex=1.3)

      pos2[,8]="Yes"
      temp=pos2
      temp=temp[-(1:dim(pos2)[1]),]
      for(i in 1:dim(pos2)[1]){
        n=strsplit(pos2[i,3],",")[[1]]
        if(grepl("Off-target",pos2[i,3])){pos2[i,8]="No"}
        for (j in 1:length(n)){
          new=pos2[i,]
          new[3]=n[j]
          new[6]=strsplit(n[j],"_")[[1]][1]
          temp=rbind(temp,new)
        }
      }
      pos2=temp
      pos2=pos2[which(grepl("NM",pos2[,3])==T),]

      pos2[,7]="No"
      for (i in 1:dim(pos2)[1]){
        if(grepl("part",pos2[i,3])==T){pos2[i,7]="Yes"}
      }
      for (i in 1:dim(pos2)[1]){
        if(pos2[i,7]=="Yes"){
          pos2[which(pos2[,1]==pos2[i,1] & pos2[,2]==pos2[i,2]),7]="Yes"
        }
      }

      for (i in 1:dim(pos2)[1]){
        pos2[i,6]=strsplit(strsplit(pos2[i,3],"exon")[[1]][2],"-")[[1]][1]
      }

      genes=pos2

      if(dim(pos2)[1]>0){
        pos2[,5]=pos2[,8]
        for(i in 1:dim(pos2)[1]){
          genes[i,4]=paste(strsplit(pos2[,3],"_")[[i]][1:3],collapse="_")
          iso=paste(strsplit(pos2[,3],"_")[[i]][2:3],collapse="_")
          if(is.element(iso,as.vector(cano)[[1]])){
            genes[i,4]=paste(genes[i,4],"*",collapse="",sep="")
          }
        }

        listgenes=unique(genes[which(genes[,4]!=""),4])
        if(UseCano==TRUE){listgenes=listgenes[which(grepl("\\*",listgenes)==T)]}
        
        low=mi-(ma-mi)*0.3
        high=mi
        
        posg=seq(low,high,length.out=length(listgenes)+2)
        posg=posg[2:(length(posg)-1)]



        niso=length(listgenes)
        nreg=length(sel)

        if(niso<4){
          policeexon=0.6
          policeiso=0.6
          hexon1=20
          hexon2=25
        }

        if(niso>=4 & niso<7){
          policeexon=0.5
          policeiso=0.6
          hexon1=30
          hexon2=45
        }

        if(niso>=7 & niso<20){
          policeexon=0.4
          policeiso=0.6
          hexon1=60
          hexon2=75
        }

        if(niso>=20 & niso<30){
          policeexon=0.2
          policeiso=0.4
          hexon1=80
          hexon2=95
        }

        if(niso>=30 & niso<40){
          policeexon=0.2
          policeiso=0.3
          hexon1=100
          hexon2=115
        }

        if(niso>=40){
          policeexon=0.2
          policeiso=0.2
          hexon1=120
          hexon2=135
        }

        if(nreg>150){
          policeexon=0.1
        }

        if(nreg>70){
          policeexon=0.2
        }

        
        for(j in 1:length(listgenes)){
          pos3=pos2[which(genes[,4]==listgenes[j]),]
          pos3=pos3[!duplicated(pos3[,c('Begin','End')]),]

          for(i in 1:dim(pos3)[1]){
            posi=which(pos[,1]==pos3[i,1])
            center=posg[j]
            if(i>1){
              lines(c(posi,posiold),c(center,center),col=j)
            }
            posiold=posi
          }
          
          for(i in 1:dim(pos3)[1]){
            posi=which(pos[,1]==pos3[i,1])
            center=posg[j]
            b3=0.3
            t=(high-low)/15
            if(pos3[i,5]=="Yes"){
              h=(high-low)/hexon1
              polygon(c(posi-b3,posi+b3,posi+b3,posi-b3,posi-b3),c(center-h,center-h,center+h,center+h,center-h),col=j,border=NA)
              if(pos3[i,7]!="Yes"){
                polygon(c(posi-b3,posi+b3,posi+b3,posi-b3,posi-b3),c(center-h,center-h,center+h,center+h,center-h),col=j,border=NA)
              }
              if(i>1){if(pos3[i,7]=="Yes" & pos3[i-1,7]=="Yes" & pos3[i,6]==pos3[i-1,6]){
                polygon(c(posi-b3-0.5,posi+b3,posi+b3,posi-b3-0.5,posi-b3-0.5),c(center-h,center-h,center+h,center+h,center-h),col=j,border=NA)
              }}
              if(pos3[i,7]!="Yes"){
                text(posi,center,pos3[i,6],cex=policeexon,adj=c(0.5,0.5),col="lightgray")
              } else {
                n=which(pos3[,6]==pos3[i,6])
                n1=which(pos[,1]==pos3[which(pos3[,6]==pos3[i,6])[1],1])
                n2=which(pos[,1]==pos3[which(pos3[,6]==pos3[i,6])[length(which(pos3[,6]==pos3[i,6]))],1])
                text((n2+n1)/2,center,pos3[i,6],cex=policeexon,adj=c(0.5,0.5),col="lightgray")
              }
            }
            if(pos3[i,5]=="No"){
              h=(high-low)/hexon2
              if(pos3[i,7]!="Yes"){
                polygon(c(posi-b3,posi+b3,posi+b3,posi-b3,posi-b3),c(center-h,center-h,center+h,center+h,center-h),col=j,border=NA)
              }
              if(i>1){if(pos3[i,7]=="Yes" & pos3[i-1,7]=="Yes" & pos3[i,6]==pos3[i-1,6]){
                polygon(c(posi-b3-0.5,posi+b3,posi+b3,posi-b3-0.5,posi-b3-0.5),c(center-h,center-h,center+h,center+h,center-h),col=j,border=NA)
              }}
              if(pos3[i,7]!="Yes"){
                text(posi,center,pos3[i,6],cex=policeexon,adj=c(0.5,0.5),col="lightgray")
              } else {
                n=which(pos3[,6]==pos3[i,6])
                n1=which(pos[,1]==pos3[which(pos3[,6]==pos3[i,6])[1],1])
                n2=which(pos[,1]==pos3[which(pos3[,6]==pos3[i,6])[length(which(pos3[,6]==pos3[i,6]))],1])
                text((n2+n1)/2,center,pos3[i,6],cex=policeexon,adj=c(0.5,0.5),col="lightgray")
              }
            }
            
            posiold=posi
          }
          
          posi=which(is.element(pos[,1],pos3[,1]))
          if(grepl("\\*",listgenes[j])==T){
              text(0-(length(pos[,4])*0.24)/2,posg[j],listgenes[j],cex=policeiso,adj=c(0.5,0.5),col=j,font=2)
            } else {
              text(0-(length(pos[,4])*0.24)/2,posg[j],listgenes[j],cex=policeiso,adj=c(0.5,0.5),col=j)
            }
          abline(h=ma)
          abline(h=mi-(ma-mi)*0.3)
        }
      }
    }

    s1=which(ClinVar[,1]==chr & as.numeric(ClinVar[,2])>=min(as.numeric(d[,2])) & as.numeric(ClinVar[,3])<=max(as.numeric(d[,3])) & grepl("deletion",ClinVar[,4])==T)
    s2=which(ClinVar[,1]==chr & as.numeric(ClinVar[,2])>=min(as.numeric(d[,2])) & as.numeric(ClinVar[,3])<=max(as.numeric(d[,3])) & grepl("duplication",ClinVar[,4])==T)
    s3=which(gnomAD[,1]==chr & as.numeric(gnomAD[,2])>=min(as.numeric(d[,2])) & as.numeric(gnomAD[,3])<=max(as.numeric(d[,3])) & grepl("DEL",gnomAD[,4])==T)
    s3=c(s3,which(gnomAD[,1]==chr & as.numeric(gnomAD[,3])>=min(as.numeric(d[,2])) & as.numeric(gnomAD[,3])<=max(as.numeric(d[,3])) & grepl("DEL",gnomAD[,4])==T))
    s3=c(s3,which(gnomAD[,1]==chr & as.numeric(gnomAD[,2])>=min(as.numeric(d[,2])) & as.numeric(gnomAD[,2])<=max(as.numeric(d[,3])) & grepl("DEL",gnomAD[,4])==T))
    s3=c(s3,which(gnomAD[,1]==chr & as.numeric(gnomAD[,2])<=min(as.numeric(d[,2])) & as.numeric(gnomAD[,3])>=max(as.numeric(d[,3])) & grepl("DEL",gnomAD[,4])==T))
    s4=which(gnomAD[,1]==chr & as.numeric(gnomAD[,2])>=min(as.numeric(d[,2])) & as.numeric(gnomAD[,3])<=max(as.numeric(d[,3])) & grepl("DUP",gnomAD[,4])==T)
    s4=c(s4,which(gnomAD[,1]==chr & as.numeric(gnomAD[,3])>=min(as.numeric(d[,2])) & as.numeric(gnomAD[,3])<=max(as.numeric(d[,3])) & grepl("DUP",gnomAD[,4])==T))
    s4=c(s4,which(gnomAD[,1]==chr & as.numeric(gnomAD[,2])>=min(as.numeric(d[,2])) & as.numeric(gnomAD[,2])<=max(as.numeric(d[,3])) & grepl("DUP",gnomAD[,4])==T))
    s4=c(s4,which(gnomAD[,1]==chr & as.numeric(gnomAD[,2])<=min(as.numeric(d[,2])) & as.numeric(gnomAD[,3])>=max(as.numeric(d[,3])) & grepl("DUP",gnomAD[,4])==T))

    s1=unique(s1)
    s2=unique(s2)
    s3=unique(s3)
    s4=unique(s4)

    s=c(s1,s2,s3,s4)
    s=cbind(s,s,s,s,s,s)
    s[,1]=c(ClinVar[s1,2],ClinVar[s2,2],gnomAD[s3,2],gnomAD[s4,2])
    s[,2]=c(ClinVar[s1,3],ClinVar[s2,3],gnomAD[s3,3],gnomAD[s4,3])
    s[,3]=c(rep(1,length(s1)),rep(2,length(s2)),rep(3,length(s3)),rep(4,length(s4)))


    if(length(s)>0){
      for (i in 1:dim(s)[1]){
        s[i,4]=suppressWarnings(min(which(as.numeric(d[,2])>as.numeric(s[i,1]) & as.numeric(d[,3])<as.numeric(s[i,2]))))
        s[i,5]=suppressWarnings(max(which(as.numeric(d[,2])>as.numeric(s[i,1]) & as.numeric(d[,3])<as.numeric(s[i,2]))))
      }
      if(length(which(s[,4]=="Inf"))>0){s=s[-which(s[,4]=="Inf"),]}
    }

    l0=2

    if(length(s)>6){
      if(dim(s)[1]>10){l0=1}
      yc=seq(from=mi-(ma-mi)*0.5,to=mi-(ma-mi)*0.3,by=(mi-(ma-mi)*0.3-(mi-(ma-mi)*0.5))/(dim(s)[1]+3))
      for (i in 1:dim(s)[1]){
        lines(c(s[i,4]-0.5,s[i,5]+0.5),c(yc[i+2],yc[i+2]),col=s[i,3],lwd=l0)
      }
    }
    if(length(s)==6){
      yc=seq(from=mi-(ma-mi)*0.5,to=mi-(ma-mi)*0.3,by=(mi-(ma-mi)*0.3-(mi-(ma-mi)*0.5))/(2))
      i=1
      lines(c(s[4]-0.5,s[5]+0.5),c(yc[2],yc[2]),col=s[3],lwd=2)
    }

    legend(0-length(pos[,4])*0.23,mi-(ma-mi)*0.4,yjust=0.5,cex=0.8,legend=c("ClinVar patho. del.","ClinVar patho. dup.","gnomAD >1% del.","gnomAD >1% dup."),col=c(1,2,3,4),pch=c(NA,NA,NA,NA),lwd=c(2,2,2,2),bg="white",y.intersp=1.5)

    
  }

  invisible(dev.off())
}

# graphical representation
plotcnvall <- function(chr,begin,end,ID,data,pdf,side,offtar,UseCano) {

  load(data)

  all=allPLOT

  if(offtar==F){all=all[which(grepl("Off-target",all[,4])==F),]}

  pdf(file=pdf,width=10,height=12)
  par(mfrow=c(1,1))

  {
    # taking regions of interest
    sel=which(all[,1]==chr & as.numeric(all[,2])>=begin & as.numeric(all[,3])<=end)
    sel=max(sel[1]-side,0):min(dim(all)[1],(sel[length(sel)]+side))

    d=all[sel,]

    if(length(which(d[,12]==0))>0) {d=d[-which(d[,12]==0),]}

    pos=unique(d[,c(2,3,4)])
    pos=cbind(pos,pos,pos)
    pos[,4]=0
    pos[,5]=0
    pos[,6]=0
    pos[,7]=0
    pos[,8]=0
    pos[,9]=0
    pos=as.data.frame(pos)
    for (i in c(1,2,4:8)){
      pos[,i]=as.numeric(pos[,i])
    }
    for(i in 1:dim(pos)[1]){
      # corrected counts
      pos[i,4]=as.numeric(d[which(as.numeric(d[,2])==pos[i,1]),11]) # ID
      pos[i,5]=as.numeric(d[which(as.numeric(d[,2])==pos[i,1]),12]) # mean
      pos[i,6]=as.numeric(d[which(as.numeric(d[,2])==pos[i,1]),13]) # SD
      pos[i,7]=min(pos[i,5]-pos[i,6],pos[i,4])/pos[i,5]
      pos[i,8]=max(pos[i,5]+pos[i,6],pos[i,4])/pos[i,5]
      pos[i,9]=d[which(as.numeric(d[,2])==pos[i,1]),15]
    }

    ma=max(1.5*1.05,(pos[,4]/pos[,5])*1.5)
    ma=min(ma,5)
    mi=(-0.05)
    plot((pos[,4]/pos[,5])+1000,xlim=c(0-length(pos[,4])*0.24,length(pos[,4])*1.02),xaxs="i",yaxs="i",ylim=c(mi-(ma-mi)*0.5,ma),main=paste("Sample: ",ID,"  /  Position: ",chr,":",begin,"-",end,"\nAll targets",sep=""),xlab="",ylab="",yaxt='n',xaxt='n',cex.main=1.6)

    # red box
    d1=which(pos[,1]>=begin & pos[,2]<=end)
    #polygon(c(min(d1)-0.45,max(d1)+0.45,min(d1)-0.45,max(d1)+0.45),c(-10,-10,1000,1000),col=rgb(1,0.9,1),border=NA)
    rect(min(d1)-0.45,-10,max(d1)+0.45,5,col=rgb(1,0.9,1),border=NA)
    abline(h=mi-(ma-mi)*0.5)

    # boxes for dicarded
    for (i in 1:dim(pos)[1]){
      if(pos[i,9]=="No"){
        rect(i-0.5,-1000,i+0.5,1000,col="lemonchiffon",border=NA)
      }
    }
    
    # Y axis
    ya=seq(0,10,0.5)
    ya=ya[which(ya>mi & ya<ma)]
    axis(2,at=ya)

    # sd grey shape
    l=1:dim(pos)[1]
    low=1-2*(pos[,6]/pos[,5])
    high=1+2*(pos[,6]/pos[,5])
    polygon(c(l,rev(l)),c(low,rev(high)),col="gray85",border=NA)
    #polygon(c(-1000,1000,-1000,1000),c(-10,-10,mi,mi),col="white",border=NA)
    rect(-1000,-10,1000,mi,col="white",border=NA)
    d1=which(pos[,1]>=begin & pos[,2]<=end)
    #polygon(c(min(d1)-0.45,max(d1)+0.45,min(d1)-0.45,max(d1)+0.45),c(-1000,-1000,mi,mi),col=rgb(1,0.9,1),border=NA)
    rect(min(d1)-0.45,-5,max(d1)+0.45,mi,col=rgb(1,0.9,1),border=NA)
    for (i in 1:dim(pos)[1]){
      if(pos[i,9]=="No"){
        rect(i-0.5,-1000,i+0.5,mi,col="lemonchiffon",border=NA)
      }
    }

    abline(v=0-length(pos[,4])*0.24)
    abline(v=length(pos[,4])*1.02)
    abline(h=mi-(ma-mi)*0.5)
    abline(h=1)
    
    for (i in 1:length(ya)){
      if(ya[i]!=1){abline(h=ya[i],lty=2,col=1)}
    }
    
    lines(pos[,4]/pos[,5], lty = 2)
    abline(h=mi)
    tar=pos[which(grepl("Off-target",pos[,3])==F & pos[,9]=="Yes"),]
    anti=pos[which(grepl("Off-target",pos[,3])==T & pos[,9]=="Yes"),]
    tarD=pos[which(grepl("Off-target",pos[,3])==F & pos[,9]=="No"),]
    antiD=pos[which(grepl("Off-target",pos[,3])==T & pos[,9]=="No"),]

    points(which(grepl("Off-target",pos[,3])==F & pos[,9]=="Yes"),tar[,4]/tar[,5],pch=15,col="green")
    points(which(grepl("Off-target",pos[,3])==T & pos[,9]=="Yes"),anti[,4]/anti[,5],pch=17,col="blue")
    points(which(grepl("Off-target",pos[,3])==F & pos[,9]=="No"),tarD[,4]/tarD[,5],pch=15,col="red")
    points(which(grepl("Off-target",pos[,3])==T & pos[,9]=="No"),antiD[,4]/antiD[,5],pch=17,col="orange")

    legend(0-length(pos[,4])*0.23,ma/2,yjust=0.5,cex=0.75,legend=c("On-targets","Off-targets","Low qual. on-targets","Low qual. off-targets","Selected sample","2*SD controls","CNV region","Low quality"),col=c("green","blue","red","orange","black","darkgrey",rgb(1,0.9,1),"lemonchiffon"),pch=c(15,17,15,17,NA,15,15,15),pt.cex=c(1,1,1,1,NA,2,2,2),lty=c(NA,NA,NA,NA,2,NA,NA,NA),bg="white",y.intersp=1.5)
    
    pos2=pos[which(grepl("NM",pos[,3])==T),]
    
    if(dim(pos2)[1]>0){
      
      mtext('Observed / Expected read ratio', side=2, line=2.5,at=mi+(ma-mi)/2,cex=1.3)
      mtext('Transcripts', side=2, line=2.5,at=mi-(ma-mi)*0.15,cex=1.3)
      mtext('and exons', side=2, line=1,at=mi-(ma-mi)*0.15,cex=1.3)
      mtext("* indicates canonical isoforms", side=2, line=0.25,at=mi-(ma-mi)*0.15,cex=0.5)
      mtext('ClinVar and', side=2, line=2.5,at=mi-(ma-mi)*0.4,cex=1.3)
      mtext('gnomAD CNVs', side=2, line=1,at=mi-(ma-mi)*0.4,cex=1.3)

      pos2[,8]="Yes"
      temp=pos2
      temp=temp[-(1:dim(pos2)[1]),]
      for(i in 1:dim(pos2)[1]){
        n=strsplit(pos2[i,3],",")[[1]]
        if(grepl("Off-target",pos2[i,3])){pos2[i,8]="No"}
        for (j in 1:length(n)){
          new=pos2[i,]
          new[3]=n[j]
          new[6]=strsplit(n[j],"_")[[1]][1]
          temp=rbind(temp,new)
        }
      }
      pos2=temp
      pos2=pos2[which(grepl("NM",pos2[,3])==T),]

      pos2[,7]="No"
      for (i in 1:dim(pos2)[1]){
        if(grepl("part",pos2[i,3])==T){pos2[i,7]="Yes"}
      }
      for (i in 1:dim(pos2)[1]){
        if(pos2[i,7]=="Yes"){
          pos2[which(pos2[,1]==pos2[i,1] & pos2[,2]==pos2[i,2]),7]="Yes"
        }
      }

      for (i in 1:dim(pos2)[1]){
        pos2[i,6]=strsplit(strsplit(pos2[i,3],"exon")[[1]][2],"-")[[1]][1]
      }

      genes=pos2

      if(dim(pos2)[1]>0){
        pos2[,5]=pos2[,8]
        for(i in 1:dim(pos2)[1]){
          genes[i,4]=paste(strsplit(pos2[,3],"_")[[i]][1:3],collapse="_")
          iso=paste(strsplit(pos2[,3],"_")[[i]][2:3],collapse="_")
          if(is.element(iso,as.vector(cano)[[1]])){
            genes[i,4]=paste(genes[i,4],"*",collapse="",sep="")
          }
        }

        listgenes=unique(genes[which(genes[,4]!=""),4])
        if(UseCano==TRUE){listgenes=listgenes[which(grepl("\\*",listgenes)==T)]}
        
        low=mi-(ma-mi)*0.3
        high=mi
        
        posg=seq(low,high,length.out=length(listgenes)+2)
        posg=posg[2:(length(posg)-1)]



        niso=length(listgenes)
        nreg=length(sel)

        if(niso<4){
          policeexon=0.6
          policeiso=0.6
          hexon1=20
          hexon2=25
        }

        if(niso>=4 & niso<7){
          policeexon=0.5
          policeiso=0.6
          hexon1=30
          hexon2=45
        }

        if(niso>=7 & niso<20){
          policeexon=0.4
          policeiso=0.6
          hexon1=60
          hexon2=75
        }

        if(niso>=20 & niso<30){
          policeexon=0.2
          policeiso=0.4
          hexon1=80
          hexon2=95
        }

        if(niso>=30 & niso<40){
          policeexon=0.2
          policeiso=0.3
          hexon1=100
          hexon2=115
        }

        if(niso>=40){
          policeexon=0.2
          policeiso=0.2
          hexon1=120
          hexon2=135
        }

        if(nreg>150){
          policeexon=0.1
        }

        if(nreg>70){
          policeexon=0.2
        }

        
        for(j in 1:length(listgenes)){
          pos3=pos2[which(genes[,4]==listgenes[j]),]
          pos3=pos3[!duplicated(pos3[,c('Begin','End')]),]

          for(i in 1:dim(pos3)[1]){
            posi=which(pos[,1]==pos3[i,1])
            center=posg[j]
            if(i>1){
              lines(c(posi,posiold),c(center,center),col=j)
            }
            posiold=posi
          }
          
          for(i in 1:dim(pos3)[1]){
            posi=which(pos[,1]==pos3[i,1])
            center=posg[j]
            b3=0.3
            t=(high-low)/15
            if(pos3[i,5]=="Yes"){
              h=(high-low)/hexon1
              polygon(c(posi-b3,posi+b3,posi+b3,posi-b3,posi-b3),c(center-h,center-h,center+h,center+h,center-h),col=j,border=NA)
              if(pos3[i,7]!="Yes"){
                polygon(c(posi-b3,posi+b3,posi+b3,posi-b3,posi-b3),c(center-h,center-h,center+h,center+h,center-h),col=j,border=NA)
              }
              if(i>1){if(pos3[i,7]=="Yes" & pos3[i-1,7]=="Yes" & pos3[i,6]==pos3[i-1,6]){
                polygon(c(posi-b3-0.5,posi+b3,posi+b3,posi-b3-0.5,posi-b3-0.5),c(center-h,center-h,center+h,center+h,center-h),col=j,border=NA)
              }}
              if(pos3[i,7]!="Yes"){
                text(posi,center,pos3[i,6],cex=policeexon,adj=c(0.5,0.5),col="lightgray")
              } else {
                n=which(pos3[,6]==pos3[i,6])
                n1=which(pos[,1]==pos3[which(pos3[,6]==pos3[i,6])[1],1])
                n2=which(pos[,1]==pos3[which(pos3[,6]==pos3[i,6])[length(which(pos3[,6]==pos3[i,6]))],1])
                text((n2+n1)/2,center,pos3[i,6],cex=policeexon,adj=c(0.5,0.5),col="lightgray")
              }
            }
            if(pos3[i,5]=="No"){
              h=(high-low)/hexon2
              if(pos3[i,7]!="Yes"){
                polygon(c(posi-b3,posi+b3,posi+b3,posi-b3,posi-b3),c(center-h,center-h,center+h,center+h,center-h),col=j,border=NA)
              }
              if(i>1){if(pos3[i,7]=="Yes" & pos3[i-1,7]=="Yes" & pos3[i,6]==pos3[i-1,6]){
                polygon(c(posi-b3-0.5,posi+b3,posi+b3,posi-b3-0.5,posi-b3-0.5),c(center-h,center-h,center+h,center+h,center-h),col=j,border=NA)
              }}
              if(pos3[i,7]!="Yes"){
                text(posi,center,pos3[i,6],cex=policeexon,adj=c(0.5,0.5),col="lightgray")
              } else {
                n=which(pos3[,6]==pos3[i,6])
                n1=which(pos[,1]==pos3[which(pos3[,6]==pos3[i,6])[1],1])
                n2=which(pos[,1]==pos3[which(pos3[,6]==pos3[i,6])[length(which(pos3[,6]==pos3[i,6]))],1])
                text((n2+n1)/2,center,pos3[i,6],cex=policeexon,adj=c(0.5,0.5),col="lightgray")
              }
            }
            
            posiold=posi
          }
          
          posi=which(is.element(pos[,1],pos3[,1]))
          if(grepl("\\*",listgenes[j])==T){
              text(0-(length(pos[,4])*0.24)/2,posg[j],listgenes[j],cex=policeiso,adj=c(0.5,0.5),col=j,font=2)
            } else {
              text(0-(length(pos[,4])*0.24)/2,posg[j],listgenes[j],cex=policeiso,adj=c(0.5,0.5),col=j)
            }
          abline(h=ma)
          abline(h=mi-(ma-mi)*0.3)
        }
      }
    }

    s1=which(ClinVar[,1]==chr & as.numeric(ClinVar[,2])>=min(as.numeric(d[,2])) & as.numeric(ClinVar[,3])<=max(as.numeric(d[,3])) & grepl("deletion",ClinVar[,4])==T)
    s2=which(ClinVar[,1]==chr & as.numeric(ClinVar[,2])>=min(as.numeric(d[,2])) & as.numeric(ClinVar[,3])<=max(as.numeric(d[,3])) & grepl("duplication",ClinVar[,4])==T)
    s3=which(gnomAD[,1]==chr & as.numeric(gnomAD[,2])>=min(as.numeric(d[,2])) & as.numeric(gnomAD[,3])<=max(as.numeric(d[,3])) & grepl("DEL",gnomAD[,4])==T)
    s3=c(s3,which(gnomAD[,1]==chr & as.numeric(gnomAD[,3])>=min(as.numeric(d[,2])) & as.numeric(gnomAD[,3])<=max(as.numeric(d[,3])) & grepl("DEL",gnomAD[,4])==T))
    s3=c(s3,which(gnomAD[,1]==chr & as.numeric(gnomAD[,2])>=min(as.numeric(d[,2])) & as.numeric(gnomAD[,2])<=max(as.numeric(d[,3])) & grepl("DEL",gnomAD[,4])==T))
    s3=c(s3,which(gnomAD[,1]==chr & as.numeric(gnomAD[,2])<=min(as.numeric(d[,2])) & as.numeric(gnomAD[,3])>=max(as.numeric(d[,3])) & grepl("DEL",gnomAD[,4])==T))
    s4=which(gnomAD[,1]==chr & as.numeric(gnomAD[,2])>=min(as.numeric(d[,2])) & as.numeric(gnomAD[,3])<=max(as.numeric(d[,3])) & grepl("DUP",gnomAD[,4])==T)
    s4=c(s4,which(gnomAD[,1]==chr & as.numeric(gnomAD[,3])>=min(as.numeric(d[,2])) & as.numeric(gnomAD[,3])<=max(as.numeric(d[,3])) & grepl("DUP",gnomAD[,4])==T))
    s4=c(s4,which(gnomAD[,1]==chr & as.numeric(gnomAD[,2])>=min(as.numeric(d[,2])) & as.numeric(gnomAD[,2])<=max(as.numeric(d[,3])) & grepl("DUP",gnomAD[,4])==T))
    s4=c(s4,which(gnomAD[,1]==chr & as.numeric(gnomAD[,2])<=min(as.numeric(d[,2])) & as.numeric(gnomAD[,3])>=max(as.numeric(d[,3])) & grepl("DUP",gnomAD[,4])==T))

    s1=unique(s1)
    s2=unique(s2)
    s3=unique(s3)
    s4=unique(s4)

    s=c(s1,s2,s3,s4)
    s=cbind(s,s,s,s,s,s)
    s[,1]=c(ClinVar[s1,2],ClinVar[s2,2],gnomAD[s3,2],gnomAD[s4,2])
    s[,2]=c(ClinVar[s1,3],ClinVar[s2,3],gnomAD[s3,3],gnomAD[s4,3])
    s[,3]=c(rep(1,length(s1)),rep(2,length(s2)),rep(3,length(s3)),rep(4,length(s4)))


    if(length(s)>0){
      for (i in 1:dim(s)[1]){
        s[i,4]=suppressWarnings(min(which(as.numeric(d[,2])>as.numeric(s[i,1]) & as.numeric(d[,3])<as.numeric(s[i,2]))))
        s[i,5]=suppressWarnings(max(which(as.numeric(d[,2])>as.numeric(s[i,1]) & as.numeric(d[,3])<as.numeric(s[i,2]))))
      }
      if(length(which(s[,4]=="Inf"))>0){s=s[-which(s[,4]=="Inf"),]}
    }

    l0=2

    if(length(s)>6){
      if(dim(s)[1]>10){l0=1}
      yc=seq(from=mi-(ma-mi)*0.5,to=mi-(ma-mi)*0.3,by=(mi-(ma-mi)*0.3-(mi-(ma-mi)*0.5))/(dim(s)[1]+3))
      for (i in 1:dim(s)[1]){
        lines(c(s[i,4]-0.5,s[i,5]+0.5),c(yc[i+2],yc[i+2]),col=s[i,3],lwd=l0)
      }
    }
    if(length(s)==6){
      yc=seq(from=mi-(ma-mi)*0.5,to=mi-(ma-mi)*0.3,by=(mi-(ma-mi)*0.3-(mi-(ma-mi)*0.5))/(2))
      i=1
      lines(c(s[4]-0.5,s[5]+0.5),c(yc[2],yc[2]),col=s[3],lwd=2)
    }

    legend(0-length(pos[,4])*0.23,mi-(ma-mi)*0.4,yjust=0.5,cex=0.8,legend=c("ClinVar patho. del.","ClinVar patho. dup.","gnomAD >1% del.","gnomAD >1% dup."),col=c(1,2,3,4),pch=c(NA,NA,NA,NA),lwd=c(2,2,2,2),bg="white",y.intersp=1.5)

    
  }

  invisible(dev.off())
}


# find real CNVs
PCAcorrection <- function(remPC,PCAmodel,dataCOR2tCOR,sdCOR2,meanCOR2,pat,dataCOR,numCOR){

  projection = predict(PCAmodel,dataCOR2tCOR)
  if(length(remPC)>0){
    dataPC4 <- projection[,remPC] %*% t(PCAmodel$rotation[,remPC])
    dataCOR3tCOR=dataCOR2tCOR-dataPC4
  } else {
    dataCOR3tCOR=dataCOR2tCOR
  }

  # rescaling
  dataPC5=dataCOR3tCOR
  for(i in 1:(numCOR-4)){
    dataPC5[i,]=dataCOR3tCOR[i,]*sdCOR2+meanCOR2
  }
  dataPC5t=t(dataPC5)
  numPC5=dim(dataPC5t)[2]
  
  # computing z-deviation from other samples and ratio
  meanCP5=apply(dataPC5t,1,mean)
  sdPC5=apply(dataPC5t,1,sd)
  z=dataPC5t
  ratio=dataPC5t
  for(i in 1:(numPC5)){
    # computing z-score and ratio
    z[,i]=(dataPC5t[,i]-meanCP5)/sdPC5
    ratio[,i]=dataPC5t[,i]
  }
  
  # removing z-scores above or below 2 (95%)
  dataPC6=dataPC5t
  for(i in 1:(numPC5)){
    dataPC6[which(z[,i]>2 | z[,i]<(-2)),i]=NA
  }
  
  # recomputing z-scores without outliers
  meanPC6=apply(dataPC6,1,mean,na.rm=T)
  sdPC6=apply(dataPC6,1,sd,na.rm=T)
  z2=dataPC5t
  ratio2=dataPC5t
  for(i in 1:(numPC5)){
    z2[,i]=(dataPC5t[,i]-meanPC6)/sdPC6
    ratio2[,i]=dataPC5t[,i]/meanPC6
  }

  sel=which(colnames(dataCOR)==pat)-4

  # computation of p-value of binomial distribution
  binom=1:dim(z2)[1]
  binom[1:dim(z2)[1]]=1
  z3=abs(z2[,sel])
  binom=2*pnorm(-abs(z3))

  # putting output in nice format
  all <- array(0, dim=c(dim(dataPC5t)[1],20))
  colnames(all)=c("Chromosome","Begin","End","Targets","ID","Counts_unc","Average_unc","SD_unc","Ratio_unc","Z-score_cor","Counts","Average_counts","SD_counts","Ratio","P-val_cor","Nb_targets","PCs_removed","Nb_CNV_HQ","Nb_samples_used","Rank_sample")
  meanCOR=apply(dataCOR[,5:numCOR],1,mean)
  meanPC6=apply(dataPC6,1,mean,na.rm=T)

  sel=which(colnames(dataCOR)==pat)
  
  all[,1]=dataCOR[,1] # chr
  all[,2]=dataCOR[,2] # begin
  all[,3]=dataCOR[,3] # end
  all[,4]=dataCOR[,4] # name target
  all[,5]=colnames(dataCOR)[which(colnames(dataCOR)==pat)] # name sample
  all[,6]=round(dataCOR[,sel],digits=0) # raw counts
  all[,7]=round(meanCOR,digits=0) # average raw counts
  all[,8]=round(apply(dataCOR[,5:numCOR],1,sd),digits=2) # SD raw counts
  all[,9]=round(as.numeric(all[,6])/as.numeric(all[,7]),digits=3) # raw ratio
  all[,10]=round(z2[,sel-4],digits=2) # corrected z-score
  all[,11]=round(dataPC5t[,sel-4],digits=0) # corrected counts
  all[,12]=round(meanPC6,digits=0) # corrected average counts
  all[,13]=round(apply(dataPC6,1,sd,na.rm=T),digits=2) # corrected SD counts
  all[,14]=round(as.numeric(all[,11])/as.numeric(all[,12]),digits=3) # corrected ratio
  all[,15]=signif(binom,digits=2) # corrected binomial p-value
  all[,16]=1 # then putting number of targets for later
  all[,17]=max(remPC) # number of removed PCs
  all[,18]=length(which(abs(as.numeric(all[,10]))>6)) # number of CNVs with abs(Z-score) > 6
  all[,19]=dim(dataPC5t)[2] # number samples taken as reference
  all[,20]=1 # rank in the sample for later

  all=as.data.frame(all)
    
  return(all)
}

# HMM function from ExomeDepth
viterbi.hmm <- function(transitions, loglikelihood, positions, expected.CNV.length) {
  if ( nrow(transitions) != ncol(transitions) ) stop("Transition matrix is not square")
  if ( length(positions) != nrow(loglikelihood) ) {
    stop("The number of positions are not matching the number of rows of the likelihood matrix", length(positions), " and ", nrow(loglikelihood))
  }
  
  nstates <- nrow(transitions)
  nobs <- nrow(loglikelihood)
  
  res <- .Call("C_hmm", nstates, nobs, transitions, loglikelihood, positions, as.double(expected.CNV.length), PACKAGE = 'ExomeDepth')
  dimnames(res[[2]])[[2]] <- c('start.p', 'end.p', 'type', 'nexons')
  res[[2]] <- as.data.frame(res[[2]])
  names(res) <- c('Viterbi.path', 'calls')
  
  return(res)
}

# return AUC for fake CNVs
returnAUC <- function(all,fakeCNV) {
  targetDEL=unique(as.numeric(fakeCNV[which(fakeCNV[,3]=="hetdel"),6]))
  targetDUP=unique(as.numeric(fakeCNV[which(fakeCNV[,3]=="hetdup"),6]))

  posDEL=all[targetDEL,3]
  posDUP=all[targetDUP,3]

  p1=unique(as.numeric(fakeCNV[,6]))
  p1=p1[which(p1>0)]
  p2=unique(as.numeric(fakeCNV[,6]))
  p2=p2[which(p2>0)]

  neg=all[-c(p1,p2),3]

  obs=c(rep("pos",length(posDEL)),rep("neg",length(neg)))
  obs=as.factor(obs)
  z=c(as.numeric(posDEL),as.numeric(neg))
  rf.roc<-roc(obs,z,quiet=T)
  a1=auc(rf.roc)
  
  obs=c(rep("pos",length(posDUP)),rep("neg",length(neg)))
  obs=as.factor(obs)
  z=c(as.numeric(posDUP),as.numeric(neg))
  rf.roc<-roc(obs,z,quiet=T)
  a2=auc(rf.roc)
  return(c(a1,a2))
}

# Correlation between samples (heatmap and text output)
{
  # res is matrix for correlations
  num=dim(dataALL)[2]
  if(num-4<=5){
    stop(paste("You need to analyze at least 6 samples to run OFF-PEAK. You provided ",num-4," samples. Exit.",sep=""))
  }
  correlation=matrix(nrow=(num-4),ncol=(num-4))
  colnames(correlation)=colnames(dataALL)[5:dim(dataALL)[2]]
  rownames(correlation)=colnames(dataALL)[5:dim(dataALL)[2]]

  # selecting targets with signal and noise limits
  meanALL=apply(dataALL[,5:num],1,mean)
  sdALL=apply(dataALL[,5:num],1,sd)
  ok=which(meanALL>minsignal & dataALL[,1]!="chrX" & log10(sdALL/meanALL)<maxvar)

  # downsampling targets to 10'000 for quick computation
  if(length(ok)>10000){
    ok=ok[sample(1:length(ok),10000,replace=F)]
  }
  selection=dataALL[ok,5:num]

  # computing correlations
  for (i in 1:(num-4)){
    for(j in 1:(num-4)){
      correlation[i,j]=cor(selection[,i],selection[,j])
    }
  }

  # plotting heatmap of correlation and writing output to text file
  pdf(file=paste(folder,"/01_general-stats/Heatmap-correlations-all.pdf",sep=""),width=7,height=7)
  par(mfrow=c(1,1))
  heatmap.2(correlation,cexRow=0.1,cexCol=0.1,breaks=seq(min(correlation), 1, length.out=71),trace="none", col = colorRampPalette(c("red","orange", "yellow", "green"))(n = 70))
  dev.off()
  correlation2=cbind(rownames(correlation),correlation)
  correlation2=rbind(colnames(correlation2),correlation2)
  write.table(correlation2,file=paste(folder,"/01_general-stats/Pairwise-correlations-all.tsv",sep=""),quote=F,sep="\t",row.names=F,col.names=F)

  # taking maximal correlation per sample
  maxCOR=1:dim(correlation)[1]
  for (i in 1:dim(correlation)[1]){
    maxCOR[i]=max(correlation[i,-i])
  }
  maxCOR=as.matrix(maxCOR)
  rownames(maxCOR)=colnames(dataALL[,5:num])

  # plotting maximal correlations and writing output to text file
  pdf(file=paste(folder,"/01_general-stats/Maximal-correlation-per-sample.pdf",sep=""),width=14,height=7)
  par(mfrow=c(1,2))
  plot(sort(maxCOR),xlab="Samples",ylab="Maximum correlation",main="Maximal correlation per sample (sorted)")
  ns=min(20,num-4)
  xx<-barplot(maxCOR[sort(maxCOR,index.return=T)$ix[1:ns],],las=2,cex.names=0.6,xpd=FALSE,ylim=c(min(maxCOR[sort(maxCOR,index.return=T)$ix[1:ns],])/1.05,1),main="Samples with lowest maximal correlation")
  text(x = xx, y = min(maxCOR[sort(maxCOR,index.return=T)$ix[1:ns],])/1.02, label = round(maxCOR[sort(maxCOR,index.return=T)$ix[1:ns],],digits=4), pos = 1, cex = 0.8, col = "red", srt = 90)
  dev.off()
  maxCOR2=cbind(rownames(maxCOR),maxCOR)
  colnames(maxCOR2)=c("ID","maxCOR")
  write.table(maxCOR2,file=paste(folder,"/01_general-stats/Maximal-correlation-per-sample.tsv",sep=""),quote=F,sep="\t",row.names=F)

}

######################################################

# declaration of variables for output
output <- array(0, dim=c(0,20))
colnames(output)=c("Chr","begin","end","name","ind","raw-counts","raw-average","raw-sd","raw-ratio","cor-z-score","cor-counts","cor-average","cor-sd","cor-ratio","cor-p-val","rank","mean-fake","nb-Z5","nb-samples","nb-SV")
BOTHsave=output
TARsave=output
num=dim(dataALL)[2]
sampleinfo=matrix(ncol=14,nrow=0)
colnames(sampleinfo)=c("ID","PCs_removed","Nb_CNV_HQ","Nb_samples_used","Ratio-mean","Ratio-SD","Z-ratio-SD","Ratio-1%","Ratio-50%","Ratio-99%","Nb<0.25","Nb<0.75","Nb>1.25","Nb>1.75")

for(pat in colnames(dataALL)[5:num]){

  print(paste("Analyzing sample: ",pat,sep=""))

  # creating output directories
  dir.create(paste(folder,"/",pat,sep=""), showWarnings = FALSE)
  dir.create(paste(folder,"/",pat,"/plots_on-targets_off-targets",sep=""), showWarnings = FALSE)
  dir.create(paste(folder,"/",pat,"/plots_on-targets-only",sep=""), showWarnings = FALSE)
  if(genomePlots){dir.create(paste(folder,"/",pat,"/plots_genome",sep=""), showWarnings = FALSE)}
  if(chromoPlots){dir.create(paste(folder,"/",pat,"/plots_chromosomes",sep=""), showWarnings = FALSE)}
  
  print(paste("  Step 1: selecting samples with correlation > ",mincor,sep=""))

  {

    selpat=which(colnames(dataALL)==pat)
    rem=which(correlation[selpat-4,]<mincor)
    
    if(num-length(rem)-4<=15){
      mini=min(num-4-1,15)
      print(paste("WARNING: less than 15 other samples with R2>",mincor,", taking ",mini," samples with highest correlation. This can lead to decreased performances.",sep=""))
      rem=which(num-4-rank(correlation[selpat-4,])>15)
    }
    if(length(rem)>0){
      dataSEL=dataALL[,-(rem+4)]
      numSEL=dim(dataSEL)[2]
    } else {
      dataSEL=dataALL
      numSEL=dim(dataALL)[2]
    }
    print(paste("    Samples selected excluding the one analyzed: ",dim(dataSEL)[2]-5," out of ",dim(dataALL)[2]-5,sep=""))

  }

  #######

  print("  Step 2: normalize by sample and filtering of targets/antitargets")

  {

    # remove targets without signal and without variance
    meanSEL=apply(dataSEL[,5:numSEL],1,mean)
    sdSEL=apply(dataSEL[,5:numSEL],1,sd)
    taken=which(meanSEL>0 & sdSEL!=0)

    dataCOR=dataSEL[taken,]
    numCOR=dim(dataCOR)[2]

    # normalize to GC content
    meanCOR=apply(dataCOR[,5:numCOR],1,mean)
    # loop on samples
    for (i in 5:numCOR){
      # for autosomes
      sel=which(dataCOR[,1]!="chrX")
      ratio=dataCOR[sel,i]/meanCOR[sel]
      fit=lm(ratio~GC[taken[sel]])
      ratioCOR=ratio-GC[taken[sel]]*fit$coefficients[2]
      ratioCOR=ratioCOR+median(ratio)-median(ratioCOR)
      correction=ratio/ratioCOR
      correction[which(correction<0.1)]=0.1
      correction[which(correction>1.9)]=1.9
      correction[which(is.na(correction))]=1
      dataCOR[sel,i]=dataCOR[sel,i]/correction

      # for chrX
      if (length(which(dataCOR[,1]=="chrX"))>0){
        sel=which(dataCOR[,1]=="chrX")
        ratio=dataCOR[sel,i]/meanCOR[sel]
        fit=lm(ratio~GC[taken[sel]])
        ratioCOR=ratio-GC[taken[sel]]*fit$coefficients[2]
        ratioCOR=ratioCOR+median(ratio)-median(ratioCOR)
        correction=ratio/ratioCOR
        correction[which(correction<0.1)]=0.1
        correction[which(correction>1.9)]=1.9
        correction[which(is.na(correction))]=1
        dataCOR[sel,i]=dataCOR[sel,i]/correction
      }
    }

    # normalize by sample for autosomes
    auto=which(dataCOR[,1]!="chrX")
    sumCOR=apply(dataCOR[auto,5:numCOR],2,sum)
    for (i in 5:numCOR){
      dataCOR[auto,i]=dataCOR[auto,i]*max(sumCOR)/sumCOR[i-4]
    }

    # normalize by sample for chrX
    if (length(which(dataCOR[,1]=="chrX"))>0){
      chrX=which(dataCOR[,1]=="chrX")
      sumCOR=apply(dataCOR[chrX,5:numCOR],2,sum)
      for (i in 5:numCOR){
        dataCOR[chrX,i]=dataCOR[chrX,i]*max(sumCOR)/sumCOR[i-4]
      }
    }

    # filtering for signal and noise limits
    meanCOR=apply(dataCOR[,5:numCOR],1,mean)
    sdCOR=apply(dataCOR[,5:numCOR],1,sd)
    dataPLOT=dataCOR
    taken2=which(meanCOR>minsignal & log10(sdCOR/meanCOR)<maxvar)
    takenPLOT=dataCOR[taken2,4]
    dataCOR=dataCOR[taken2,]

    print(paste("    Intervals selected: ",dim(dataCOR)[1]," out of ",dim(dataSEL)[1],sep=""))

  }

  ######

  print("  Step 3: optimization of PC removal with fake CNVs")

  {

    dataPC=dataCOR
    numPC=dim(dataPC)[2]

    # downsample targets and off-targets to find optimal number of PCs to remove
    if(downsample<dim(dataPC)[1]){
      dataPC=dataPC[sample(dim(dataPC)[1],downsample,replace=F),]
    }

    # creating fake CNVs (heterozygous deletions and duplications)
    neach=nbFake
    samples=sample(dim(dataPC)[1], neach*2,replace=F)
    targets=rep(which(colnames(dataPC)==pat), neach*2)
    for (i in 1:neach){
      dataPC[(samples[i]),targets[i]]=dataPC[(samples[i]),targets[i]]*0.5
    }
    for (i in (neach+1):(neach*2)){
      dataPC[(samples[i]),targets[i]]=dataPC[(samples[i]),targets[i]]*1.5
    }
    hetdel=dataPC[samples[1:neach],4]
    hetdel=cbind(hetdel,colnames(dataPC)[targets[1:neach]])
    hetdel=cbind(hetdel,"hetdel","target",1)
    hetdup=dataPC[samples[(neach+1):(neach*2)],4]
    hetdup=cbind(hetdup,colnames(dataPC)[targets[(neach+1):(neach*2)]])
    hetdup=cbind(hetdup,"hetdup","target",1)
    fakeCNV=rbind(hetdel,hetdup)

    dataPC2=dataPC[,5:numPC]
    dataPC2t=t(dataPC2)
    dataPC2tCOR=dataPC2t

    # normalize data
    sdPC <- apply(dataPC2t, 2, sd) # calculate standard deviation
    meanPC <- apply(dataPC2t, 2, mean) # calculate mean
    for (i in 5:numPC){
      dataPC2tCOR[i-4,]=(dataPC2t[i-4,]-meanPC)/sdPC
    }

    # doing the PCA and removing analyzed sample if leaveoneout==1
    if(leaveoneout==1){
      dataPC3=dataPC2tCOR[-which(rownames(dataPC2tCOR)==pat),]
      sdPC3=apply(dataPC3,2,sd)
      dataPC3=dataPC3[,which(sdPC3>0)]
      dataPC2tCOR=dataPC2tCOR[,which(sdPC3>0)]
      dataPC2t=dataPC2t[,which(sdPC3>0)]
      dataPC=dataPC[which(sdPC3>0),]
      sdPC <- apply(dataPC2t, 2, sd) # calculate standard deviation
      meanPC <- apply(dataPC2t, 2, mean) # calculate mean
    } else {
      dataPC3=dataPC2tCOR
    }
    PCAmodel <- prcomp(dataPC3, center = TRUE,scale. = TRUE)

    # matrix for PC removal result
    PCresults<-matrix(nrow=dim(dataPC3)[1],ncol=3)
    stop=FALSE

    # remove increasing number of PCs and record AUCs of the z-score of fake CNVs
    for (pc in 0:(dim(dataPC3)[1]-1)){
      if(stop==FALSE){
        trial=1:pc
        if(pc==0){
          remPC=c()
          trial=c()
        }

        projection = predict(PCAmodel,dataPC2tCOR)
        if(length(trial)>0){
          dataPC4 <- projection[,trial] %*% t(PCAmodel$rotation[,trial])
          dataPC3tCOR=dataPC2tCOR-dataPC4
        } else {
          dataPC3tCOR=dataPC2tCOR
        }

        # rescaling for all
        dataPC5=dataPC3tCOR
        for(i in 1:(numPC-4)){
          dataPC5[i,]=dataPC3tCOR[i,]*sdPC+meanPC
        }
        dataPC5t=t(dataPC5)
        
        # computing z-deviation
        meanPC5t=apply(dataPC5t,1,mean)
        sdPC5t=apply(dataPC5t,1,sd)
        z=dataPC5t
        ratio=dataPC5t
        for(i in 1:(numPC-4)){
          # computing z-score and ratio
          z[,i]=(dataPC5t[,i]-meanPC5t)/sdPC5t
          ratio[,i]=dataPC5t[,i]
        }
        
        # removing z-scores above or below 2 (95%)
        dataPC6=dataPC5t
        for(i in 1:(numPC-4)){
          dataPC6[which(z[,i]>2 | z[,i]<(-2)),i]=NA
        }
        
        # recomputing z-scores without outliers
        meanPC6=apply(dataPC6,1,mean,na.rm=T)
        sdPC6=apply(dataPC6,1,sd,na.rm=T)
        z6=dataPC5t
        ratio6=dataPC5t
        for(i in 1:(numPC-4)){
          z6[,i]=(dataPC5t[,i]-meanPC6)/sdPC6
          ratio6[,i]=dataPC5t[,i]/meanPC6
        }

        # putting output in nice format
        all <- matrix(0, ncol=3, nrow=dim(dataPC5t)[1])
        colnames(all)=c("name","ind","z-score")
        i=which(colnames(dataPC)==pat)-4
        all[,1]=dataPC[,4]
        all[,3]=z6[,i]
        all=as.data.frame(all)

        fakeCNV=cbind(fakeCNV,fakeCNV[,5])
        fakeCNV[,6]=samples
        
        PCresults[pc+1,]=c(pc,returnAUC(all,fakeCNV))
        score=(PCresults[,2]+PCresults[,3])/2

        # stopping if score with 10 less PCs removed is higher
        if(pc>9){
          if(score[pc+1]<(score[pc-9])){
            stop=T
          }
        }
      }
    }

    score=(PCresults[,2]+PCresults[,3])/2
    maxPC=which(score>max(score,na.rm=T)-stopPC)[1]-1
    remPC=1:maxPC
    if(maxPC==0){remPC=1}
    if(stopPC==(-1)){remPC=c()}

    # plotting AUCs for different number of PC removed
    pdf(file=paste(folder,"/06_PC-plots/PC-removal-",pat,".pdf",sep=""),width=7,height=7)
    plot(PCresults[,1],PCresults[,2],ylim=c(min(PCresults[,2:3],na.rm=T),max(PCresults[,2:3],na.rm=T)),type='l',xlab="Number of removed PCs",ylab="AUC",main=paste("PC removal for ",pat,sep=""))
    lines(PCresults[,1],PCresults[,3],col=2)
    lines(PCresults[,1],(PCresults[,2]+PCresults[,3])/2,col=3)
    abline(v=max(remPC))
    legend("bottomright",cex=0.8,legend=c("Heterozygous deletions","Heterozygous duplications","Both","Number of PCs to remove"),col=c("black","red","green","blue"),pch=NA,lty=c(1,1,1,1),lwd=c(1,1,1,2))
    dev.off()

  }

  ######

  print("  Step 4: computation of CNVs")

  {

    dataCOR2=dataCOR[,5:numCOR]
    dataCOR2t=t(dataCOR2)
    dataCOR2tCOR=dataCOR2t

    # normalize data and helping detection of homozygous deletion
    sdCOR2 <- apply(dataCOR2t, 2, sd) # calculate standard deviation
    meanCOR2 <- apply(dataCOR2t, 2, mean) # calculate mean
    homdelrat=0.05
    homdelz=-3
    homdelout=-10
    sel1=which(rownames(dataCOR2t)==pat)
    sel2=which(dataCOR2t[sel1,]/meanCOR2<homdelrat & (dataCOR2t[sel1,]-meanCOR2)/sdCOR2<homdelz)
    for (i in 5:numCOR){
      dataCOR2tCOR[i-4,]=(dataCOR2t[i-4,]-meanCOR2)/sdCOR2 
    }
    dataCOR2tCOR[sel1,sel2]=homdelout

    # doing the PCA
    if(leaveoneout==1){
      dataCOR3=dataCOR2tCOR[-which(rownames(dataCOR2tCOR)==pat),]
      sdCOR3=apply(dataCOR3,2,sd)
      dataCOR3=dataCOR3[,which(sdCOR3>0)]
      dataCOR2tCOR=dataCOR2tCOR[,which(sdCOR3>0)]
      dataCOR2t=dataCOR2t[,which(sdCOR3>0)]
      dataCOR=dataCOR[which(sdCOR3>0),]
      sdCOR2 <- apply(dataCOR2t, 2, sd) # calculate standard deviation
      meanCOR2 <- apply(dataCOR2t, 2, mean) # calculate mean
    } else {
      dataCOR3=dataCOR2tCOR
      sdCOR2 <- apply(dataCOR2t, 2, sd) 
      meanCOR2 <- apply(dataCOR2t, 2, mean)
    }
    PCAmodel <- prcomp(dataCOR3, center = TRUE,scale. = TRUE)

    var = PCAmodel$sdev^2 / sum(PCAmodel$sdev^2)
    pdf(file=paste(folder,"/06_PC-plots","/Variance-explained-",pat,".pdf",sep=""),width=14,height=7)
    par(mfrow=c(1,2))
    plot(var,cex=0.5,xlab="PC number",ylab="Percent variance explained",main=paste("Variance explained for ",pat,sep=""))
    if(length(remPC)>0){abline(v=max(remPC),col=4)} else {abline(v=0,col=4)}
    legend("topright",cex=0.8,legend=c("Number of PCs to remove"),col=c("blue"),pch=NA,lty=c(1),lwd=c(2))
    if(length(remPC)>0){
      plot(var,cex=0.5,xlab="PC number",ylab="Percent variance explained",main=paste("Zoom for ",pat,sep=""),xlim=c(0,min(max(remPC)+10,dim(dataCOR)[2]-4)))
    } else {
      plot(var,cex=0.5,xlab="PC number",ylab="Percent variance explained",main=paste("Zoom for ",pat,sep=""),xlim=c(0,min(0+10,dim(dataCOR)[2]-4)))
    }
    if(length(remPC)>0){abline(v=max(remPC),col=4)} else {abline(v=0,col=4)}
    legend("topright",cex=0.8,legend=c("Number of PCs to remove"),col=c("blue"),pch=NA,lty=c(1),lwd=c(2))
    dev.off()

    # rem=PCs to remove
    all <- PCAcorrection(remPC,PCAmodel,dataCOR2tCOR,sdCOR2,meanCOR2,pat,dataCOR,numCOR)
    all=as.data.frame(all)

    # putting at zero negative corrected counts and ratios
    if(length(which(as.numeric(all[,11])<0))>0){
      all[which(as.numeric(all[,11])<0),11]=0
    }
    if(length(which(as.numeric(all[,14])<0))>0){
      all[which(as.numeric(all[,14])<0),14]=0
    }
    select=which(all[,6]==0)
    if(length(select)>0){
      all[select,11]=0
      all[select,10]=round((as.numeric(all[select,11])-as.numeric(all[select,12]))/as.numeric(all[select,13]),digits=2)
    }

    # removing targets with change of ratio more than 0.75
    F1=mean(as.numeric(all[,10]))
    F2=sd(as.numeric(all[,10]))
    m=as.numeric(all[,14])
    m[which(m<0)]=0
    all=all[which(abs(as.numeric(all[,9])-m)<0.75),]
    
    save(all,file=paste(folder,"/05_RData-files/data-",pat,".RData",sep=""))

    # computating information about sample for sample output file
    sampmean=round(mean(as.numeric(all[,14])),digits=3)
    sampsd=round(sd(as.numeric(all[,14])),digits=3)
    sampquantile=quantile(as.numeric(all[,14]),probs=c(0.01,0.5,0.99))
    sampnb=c(length(which(as.numeric(all[,14])<0.25)),length(which(as.numeric(all[,14])<0.75)),length(which(as.numeric(all[,14])>1.25)),length(which(as.numeric(all[,14])>1.75)))
    sampleinfo=rbind(sampleinfo,c(all[1,c(5,17,18,19)],sampmean,sampsd,0,sampquantile[1:3],sampnb))
    varav=mean(as.numeric(sampleinfo[,6]))
    varsd=sd(as.numeric(sampleinfo[,6]))
    sampleinfo[,7]=(as.numeric(sampleinfo[,6])-varav)/varsd
    lowqual=which(abs((as.numeric(sampleinfo[,6])-varav)/varsd)>2)
    if(length(lowqual)>0){
        HQ=as.character(sampleinfo[-lowqual,1])
    } else {
        HQ=as.character(sampleinfo[,1])
    }

    dataPLOT2=dataPLOT[,5:numCOR]
    dataPLOT2t=t(dataPLOT2)
    dataPLOT2tCOR=dataPLOT2t

    # doing the PCA
    if(leaveoneout==1){
      dataPLOT3=dataPLOT2tCOR[-which(rownames(dataPLOT2tCOR)==pat),]
      sdPLOT3=apply(dataPLOT3,2,sd)
      dataPLOT3=dataPLOT3[,which(sdPLOT3>0)]
      dataPLOT2tCOR=dataPLOT2tCOR[,which(sdPLOT3>0)]
      dataPLOT2t=dataPLOT2t[,which(sdPLOT3>0)]
      dataPLOT=dataPLOT[which(sdPLOT3>0),]
      sdPLOT2 <- apply(dataPLOT2t, 2, sd) # calculate standard deviation
      meanPLOT2 <- apply(dataPLOT2t, 2, mean) # calculate mean
    } else {
      dataPLOT3=dataPLOT2tCOR
      sdPLOT2 <- apply(dataPLOT2t, 2, sd)
      meanPLOT2 <- apply(dataPLOT2t, 2, mean)
    }
    PCAmodel <- prcomp(dataPLOT3, center = TRUE,scale. = TRUE)

    # rem=PCs to remove
    allPLOT <- PCAcorrection(remPC,PCAmodel,dataPLOT2tCOR,sdPLOT2,meanPLOT2,pat,dataPLOT,numCOR)
    allPLOT=as.data.frame(allPLOT)

    # putting at zero negative corrected counts and ratios
    if(length(which(as.numeric(allPLOT[,11])<0))>0){
      allPLOT[which(as.numeric(allPLOT[,11])<0),11]=0
    }
    if(length(which(as.numeric(allPLOT[,14])<0))>0){
      allPLOT[which(as.numeric(allPLOT[,14])<0),14]=0
    }
    select=which(allPLOT[,6]==0)
    if(length(select)>0){
      allPLOT[select,11]=0
      allPLOT[select,10]=round((as.numeric(allPLOT[select,11])-as.numeric(allPLOT[select,12]))/as.numeric(allPLOT[select,13]),digits=2)
    }

    allPLOT=allPLOT[which(is.element(allPLOT[,4],all[,4])==F),]
    allPLOT[,15]="No"

    allPLOT2=all
    allPLOT2[,15]="Yes"

    allPLOT=rbind(allPLOT,allPLOT2)

    allPLOT=allPLOT[order(as.numeric(allPLOT[,2]),decreasing=FALSE),]
    allPLOT=allPLOT[order(allPLOT[,1],decreasing=FALSE),]

    save(allPLOT,file=paste(folder,"/05_RData-files/data-plot-",pat,".RData",sep=""))

  }

  ######

  print("  Step 5: genome and chromosome plots")

  {

    # genome plot by target
    raw=dataALL[,selpat]/apply(dataALL[,5:num],1,mean)
    if(genomePlots){
      n1=20
      n2=200
      n3=5000
      pdf(file=paste(folder,"/",pat,"/plots_genome/",pat,"-genome_by-target.pdf",sep=""),width=16,height=16)
      par(mfrow=c(2,1))
      plot(-10,-10,type='l',ylim=c(0,2),xlim=c(0,dim(dataALL)[1]),col=4,xaxt='n',ylab="Observed / Expected read ratio",xlab="Chromosomes",main=paste("Genome plot by target (raw) for",pat,sep=" "),xaxs = "i",yaxs = "i")
      loc=0
      for(chr in unique(dataALL[,1])){
        temp=loc
        loc=max(which(dataALL[,1]==chr))
        text((loc+temp)/2,-0.1,chr,cex=1,xpd=NA,srt=90)
        lines((temp+1):loc,runmean(raw[which(dataALL[,1]==chr)],n1),col=4)
        lines((temp+1):loc,runmean(raw[which(dataALL[,1]==chr)],n2),col=2)
        lines((temp+1):loc,runmean(raw[which(dataALL[,1]==chr)],n3),col=3)
        abline(v=loc)
      }
      abline(v=0)
      abline(h=0.5,lty=2,col="darkgrey")
      abline(h=1,lty=2,col=1)
      abline(h=1.5,lty=2,col="darkgrey")
      legend("topright", legend = c("Average over 10","Average over 200","Average over 5000") ,pch = c(NA,NA,NA), lty = c(1,1,1), lwd = c(2,2,1),col = c(4,2,3),bg='white',y.intersp=1.2,cex=1.2,ncol=1)

      plot(-10,-10,type='l',ylim=c(0,2),xlim=c(0,dim(all)[1]),col=4,xaxt='n',ylab="Observed / Expected read ratio",xlab="Chromosomes",main=paste("Genome plot by target (filtered and corrected) for",pat,sep=" "),xaxs = "i",yaxs = "i")
      loc=0
      for(chr in unique(all[,1])){
        temp=loc
        loc=max(which(all[,1]==chr))
        lines((temp+1):loc,runmean(all[which(all[,1]==chr),14],n1),col=4)
        lines((temp+1):loc,runmean(all[which(all[,1]==chr),14],n2),col=2)
        lines((temp+1):loc,runmean(all[which(all[,1]==chr),14],n3),col=3)
        abline(v=loc)
        text((loc+temp)/2,-0.1,chr,cex=1,xpd=NA,srt=90)
      }
      abline(v=0)
      abline(h=0.5,lty=2,col="darkgrey")
      abline(h=1,lty=2,col=1)
      abline(h=1.5,lty=2,col="darkgrey")
      legend("topright", legend = c("Average over 10","Average over 200","Average over 5000") ,pch = c(NA,NA,NA), lty = c(1,1,1), lwd = c(2,2,1),col = c(4,2,3),bg='white',y.intersp=1.2,cex=1.2,ncol=1)
      dev.off()

      # genome plot by coordinate
      n1=20
      n2=200
      n3=5000
      size=0
      for (chr in unique(dataALL[,1])){
        size=size+max(as.numeric(dataALL[which(dataALL[,1]==chr),3]))
      }
      pdf(file=paste(folder,"/",pat,"/plots_genome/",pat,"-genome_by-coord.pdf",sep=""),width=16,height=16)
      par(mfrow=c(2,1))
      plot(-10,-10,type='l',ylim=c(0,2),xlim=c(0,size),col=4,xaxt='n',ylab="Observed / Expected read ratio",xlab="Chromosomes",main=paste("Genome plot by coordinate (raw) for",pat,sep=" "),xaxs = "i",yaxs = "i")
      loc=0
      for(chr in unique(dataALL[,1])){
        temp=loc
        loc=loc+max(as.numeric(dataALL[which(dataALL[,1]==chr),3]))
        text((loc+temp)/2,-0.1,chr,cex=1,xpd=NA,srt=90)
        lines(temp+as.numeric(dataALL[which(dataALL[,1]==chr),2]),runmean(raw[which(dataALL[,1]==chr)],n1),col=4)
        lines(temp+as.numeric(dataALL[which(dataALL[,1]==chr),2]),runmean(raw[which(dataALL[,1]==chr)],n2),col=2)
        lines(temp+as.numeric(dataALL[which(dataALL[,1]==chr),2]),runmean(raw[which(dataALL[,1]==chr)],n3),col=3)
        abline(v=loc)
      }
      abline(v=0)
      abline(h=0.5,lty=2,col="darkgrey")
      abline(h=1,lty=2,col=1)
      abline(h=1.5,lty=2,col="darkgrey")
      legend("topright", legend = c("Average over 10","Average over 200","Average over 5000") ,pch = c(NA,NA,NA), lty = c(1,1,1), lwd = c(2,2,1),col = c(4,2,3),bg='white',y.intersp=1.2,cex=1.2,ncol=1)

      size=0
      for (chr in unique(all[,1])){
        size=size+max(as.numeric(all[which(all[,1]==chr),3]))
      }
      plot(-10,-10,type='l',ylim=c(0,2),xlim=c(0,size),col=4,xaxt='n',ylab="Observed / Expected read ratio",xlab="Chromosomes",main=paste("Genome plot by coordinate (filtered and corrected) for",pat,sep=" "),xaxs = "i",yaxs = "i")
      loc=0
      for(chr in unique(all[,1])){
        temp=loc
        loc=loc+max(as.numeric(all[which(all[,1]==chr),3]))
        lines(temp+as.numeric(all[which(all[,1]==chr),2]),runmean(all[which(all[,1]==chr),14],n1),col=4)
        lines(temp+as.numeric(all[which(all[,1]==chr),2]),runmean(all[which(all[,1]==chr),14],n2),col=2)
        lines(temp+as.numeric(all[which(all[,1]==chr),2]),runmean(all[which(all[,1]==chr),14],n3),col=3)
        abline(v=loc)
        text((loc+temp)/2,-0.1,chr,cex=1,xpd=NA,srt=90)
      }
      abline(v=0)
      abline(h=0.5,lty=2,col="darkgrey")
      abline(h=1,lty=2,col=1)
      abline(h=1.5,lty=2,col="darkgrey")
      legend("topright", legend = c("Average over 10","Average over 200","Average over 5000") ,pch = c(NA,NA,NA), lty = c(1,1,1), lwd = c(2,2,1),col = c(4,2,3),bg='white',y.intersp=1.2,cex=1.2,ncol=1)
      dev.off()
    }

    if(chromoPlots){
      # CHROMOSOME plots by coord
      n1=20
      n2=200
      n3=5000
      for(chr in unique(dataALL[,1])){
        pdf(file=paste(folder,"/",pat,"/plots_chromosomes/",pat,"-",chr,"_by-coord.pdf",sep=""),width=16,height=16)
        par(mfrow=c(2,1))
        plot(-10,-10,type='l',ylim=c(0,2),xlim=c(0,max(as.numeric(dataALL[which(dataALL[,1]==chr),3]))),col=4,ylab="Observed / Expected read ratio",xlab="Position on chromosome",main=paste(chr,"plot by coordinate (raw) for",pat,sep=" "),xaxs = "i",yaxs = "i")
        loc=max(as.numeric(dataALL[which(dataALL[,1]==chr),3]))
        lines(as.numeric(dataALL[which(dataALL[,1]==chr),2]),runmean(raw[which(dataALL[,1]==chr)],n1),col=4)
        lines(as.numeric(dataALL[which(dataALL[,1]==chr),2]),runmean(raw[which(dataALL[,1]==chr)],n2),col=2)
        lines(as.numeric(dataALL[which(dataALL[,1]==chr),2]),runmean(raw[which(dataALL[,1]==chr)],n3),col=3)
        abline(v=0)
        abline(h=0.5,lty=2,col="darkgrey")
        abline(h=1,lty=2,col=1)
        abline(h=1.5,lty=2,col="darkgrey")
        legend("topright", legend = c("Average over 10","Average over 200","Average over 5000") ,pch = c(NA,NA,NA), lty = c(1,1,1), lwd = c(2,2,1),col = c(4,2,3),bg='white',y.intersp=1.2,cex=1.2,ncol=1)

        plot(-10,-10,type='l',ylim=c(0,2),xlim=c(0,max(as.numeric(all[which(all[,1]==chr),3]))),col=4,ylab="Observed / Expected read ratio",xlab="Position on chromosome",main=paste(chr,"plot by coordinate (filtered and corrected) for",pat,sep=" "),xaxs = "i",yaxs = "i")
        loc=max(as.numeric(all[which(all[,1]==chr),3]))
        lines(as.numeric(all[which(all[,1]==chr),2]),runmean(all[which(all[,1]==chr),14],n1),col=4)
        lines(as.numeric(all[which(all[,1]==chr),2]),runmean(all[which(all[,1]==chr),14],n2),col=2)
        lines(as.numeric(all[which(all[,1]==chr),2]),runmean(all[which(all[,1]==chr),14],n3),col=3)
        abline(v=0)
        abline(h=0.5,lty=2,col="darkgrey")
        abline(h=1,lty=2,col=1)
        abline(h=1.5,lty=2,col="darkgrey")
        legend("topright", legend = c("Average over 10","Average over 200","Average over 5000") ,pch = c(NA,NA,NA), lty = c(1,1,1), lwd = c(2,2,1),col = c(4,2,3),bg='white',y.intersp=1.2,cex=1.2,ncol=1)
        dev.off()
      }

      # CHROMOSOME plots by targets
      n1=20
      n2=200
      n3=5000
      for(chr in unique(dataALL[,1])){
        pdf(file=paste(folder,"/",pat,"/plots_chromosomes/",pat,"-",chr,"_by-target.pdf",sep=""),width=16,height=16)
        par(mfrow=c(2,1))
        loc=length(which(dataALL[,1]==chr))
        plot(-10,-10,type='l',ylim=c(0,2),xlim=c(0,loc),col=4,xaxt='n',ylab="Observed / Expected read ratio",xlab=chr,main=paste(chr,"plot by target (raw) for",pat,sep=" "),xaxs = "i",yaxs = "i")
        lines(1:loc,runmean(raw[which(dataALL[,1]==chr)],n1),col=4)
        lines(1:loc,runmean(raw[which(dataALL[,1]==chr)],n2),col=2)
        lines(1:loc,runmean(raw[which(dataALL[,1]==chr)],n3),col=3)
        abline(v=0)
        abline(h=0.5,lty=2,col="darkgrey")
        abline(h=1,lty=2,col=1)
        abline(h=1.5,lty=2,col="darkgrey")
        legend("topright", legend = c("Average over 10","Average over 200","Average over 5000") ,pch = c(NA,NA,NA), lty = c(1,1,1), lwd = c(2,2,1),col = c(4,2,3),bg='white',y.intersp=1.2,cex=1.2,ncol=1)
        
        loc=length(which(all[,1]==chr))
        plot(-10,-10,type='l',ylim=c(0,2),xlim=c(0,loc),col=4,xaxt='n',ylab="Observed / Expected read ratio",xlab=chr,main=paste(chr,"plot by target (filtered and corrected) for",pat,sep=" "),xaxs = "i",yaxs = "i")
        lines(1:loc,runmean(all[which(all[,1]==chr),14],n1),col=4)
        lines(1:loc,runmean(all[which(all[,1]==chr),14],n2),col=2)
        lines(1:loc,runmean(all[which(all[,1]==chr),14],n3),col=3)
        abline(v=0)
        abline(h=0.5,lty=2,col="darkgrey")
        abline(h=1,lty=2,col=1)
        abline(h=1.5,lty=2,col="darkgrey")
        legend("topright", legend = c("Average over 10","Average over 200","Average over 5000") ,pch = c(NA,NA,NA), lty = c(1,1,1), lwd = c(2,2,1),col = c(4,2,3),bg='white',y.intersp=1.2,cex=1.2,ncol=1)
        dev.off()
      }
    }

  }

  ######

  print("  Step 6: merging of consecutive CNVs for all targets")

  {

    # column for SD quality
    all=cbind(all,all[,20])
    all[,21]=0

    # cor-counts cor-average  cor-sd cor-ratio cor-p-val rank
    data=all[,c(11,12,13,14,15,16)]
    data[,4]=((as.numeric(data[,1])-as.numeric(data[,2]))/as.numeric(data[,3])) # z-score from normal
    data[,5]=((as.numeric(data[,1])-as.numeric(data[,2])/2)/as.numeric(data[,3])) # z-score from het. deletion
    data[,6]=((as.numeric(data[,1])-as.numeric(data[,2])*1.5)/as.numeric(data[,3])) # z-score from het. duplication
    data[,4]=log10(2*pnorm(-abs(data[,4]))) # log10 p-value from z-score
    data[,5]=log10(2*pnorm(-abs(data[,5]))) # log10 p-value from z-score
    data[,6]=log10(2*pnorm(-abs(data[,6]))) # log10 p-value from z-score
    data=as.data.frame(cbind(data[,4],data[,5],data[,6]))
    colnames(data)=c("cn2","cn1","cn3")

    # transition probabilities for the merging
    p1=0.00005
    p2=0.5
    trans=matrix(c(1-2*p1,p1,p1, 1-p2,p2,0 ,1-p2,0,p2),nrow=3,byrow=T)

    events=matrix(ncol=5,nrow=0)
    events=as.data.frame(events)

    # loop over chromosomes and detect consecutive CNVs with viterbi algorithm
    for(chr in unique(all[,1])){
      data2=data[which(all[,1]==chr),]
      posi=as.numeric(all[which(all[,1]==chr),2])

      len=50000

      data2[which(data2[,1]==-Inf),1]=(-10)
      data2[which(data2[,2]==-Inf),2]=(-10)
      data2[which(data2[,3]==-Inf),3]=(-10)

      my.calls <- viterbi.hmm(as.matrix(trans),as.matrix(data2),as.integer(posi),len)

      if(dim(my.calls$calls)[1]>0){
        temp=matrix(ncol=5,nrow=dim(my.calls$calls)[1])
        temp=as.data.frame(temp)
        temp[,1]=chr # chromosome
        temp[,2]=posi[my.calls$calls[,1]] # begin
        temp[,3]=posi[my.calls$calls[,2]] # end
        temp[,4]=my.calls$calls[,3] # type 1=del, 2=dup
        temp[,5]=my.calls$calls[,4] # number of targets involved
        events=rbind(events,temp)
      }
    }
    colnames(events)=c("chr","begin","end","state","nbtargets")

    # removing single target events
    events=events[which(events[,5]>1),]
    nCNV=dim(events)[1]
    allMULTI=all[1:nCNV,]

    # merging consecutive CNVs
    if(nCNV>0){
      chr="chr1"
      suball=all[which(all[,1]==chr),]
      for (i in 1:nCNV){
        if(events[i,1]!=chr){
          chr=events[i,1]
          suball=all[which(all[,1]==chr),]
        }
        # targets inside the interval
        sel=which(suball[,1]==events[i,1] & as.numeric(suball[,2])>=as.numeric(events[i,2]) & as.numeric(suball[,2])<=as.numeric(events[i,3]))
        allMULTI[i,1]=suball[sel[1],1] # chromosome
        allMULTI[i,2]=min(as.numeric(suball[sel,2])) # begin
        allMULTI[i,3]=max(as.numeric(suball[sel,3])) # end
        order=sort(as.numeric(suball[sel,2]),index.return=T)$ix
        allMULTI[i,4]=paste(suball[sel[order],4],sep=",",collapse=";") # name
        allMULTI[i,5]=suball[sel[1],5] # sample ID
        allMULTI[i,6]=sum(as.numeric(suball[sel,6])) # raw-count
        allMULTI[i,7]=sum(as.numeric(suball[sel,7])) # raw-average
        allMULTI[i,8]=round(sqrt(mean((as.numeric(suball[sel,8])^2))),digits=2) # raw-sd
        allMULTI[i,9]=round(as.numeric(allMULTI[i,6])/as.numeric(allMULTI[i,7]),digits=3) # raw-ratio
        allMULTI[i,11]=sum(as.numeric(suball[sel,11])) # cor-counts
        allMULTI[i,12]=sum(as.numeric(suball[sel,12])) # cor-average
        allMULTI[i,13]=round(sqrt(mean((as.numeric(suball[sel,13])^2))),digits=2) # cor-sd
        allMULTI[i,10]=round((as.numeric(allMULTI[i,11])-as.numeric(allMULTI[i,12]))/as.numeric(allMULTI[i,13]),digits=2) # z-score
        allMULTI[i,14]=round(as.numeric(allMULTI[i,11])/as.numeric(allMULTI[i,12]),digits=3) # cor-ratio
        allMULTI[i,15]=signif(2*pnorm(-abs(as.numeric(allMULTI[i,10]))),digits=2) # p-value (z-score)
        allMULTI[i,16]=length(sel) # number target
        allMULTI[i,17]=suball[sel[1],17] # number of PCs removed
        allMULTI[i,18]=suball[sel[1],18] # nb-z6
        allMULTI[i,19]=suball[sel[1],19] # nb-samples

        # computation of CQ quality score
        nsel=length(sel)
        nsub=dim(suball)[1]
        if(nsel>=3){
          l1=max(sel[1]-11,1)
          l2=max(sel[1]-2,1)
          l3=sel[2]
          l4=sel[nsel]-1
          l5=min(sel[nsel]+2,nsub)
          l6=min(sel[nsel]+11,nsub)
          p1=sum(abs(as.numeric(suball[l1:l2,14])-1))
          m=mean(pmax(as.numeric(suball[l3:l4,14]),0))
          p2=sum(abs(pmax(as.numeric(suball[l3:l4,14]),0)-m))
          p3=sum(abs(as.numeric(suball[l5:l6,14])-1))
          allMULTI[i,21]=p1/length(l1:l2)/2+p2/length(l3:l4)+p3/length(l5:l6)/2
        }
        if(nsel<=3){
          l1=max(sel[1]-11,1)
          l2=max(sel[1]-2,1)
          l3=sel[1]
          l4=sel[nsel]
          l5=min(sel[nsel]+2,nsub)
          l6=min(sel[nsel]+11,nsub)
          p1=sum(abs(as.numeric(suball[l1:l2,14])-1))
          m=mean(pmax(as.numeric(suball[l3:l4,14]),0))
          p2=sum(abs(pmax(as.numeric(suball[l3:l4,14]),0)-m))
          p3=sum(abs(as.numeric(suball[l5:l6,14])-1))
          allMULTI[i,21]=p1/length(l1:l2)/2+p2/length(l3:l4)+p3/length(l5:l6)/2
        }
      }
    }

    # taking only corrected Z-score above 4
    allSINGLE=all[which(abs(as.numeric(all[,10]))>((minZ*F2)+F1) | (as.numeric(all[,14])<0.1 & as.numeric(all[,9])<0.1)),]
    if(dim(allSINGLE)[1]>0){allSINGLE[,16]=1} # number of targets = 1 for isolated

    # CQ score for unitargets
    if(dim(allSINGLE)[1]>0){
      chr="chr1"
      suball=all[which(all[,1]==chr),]
      for (i in 1:dim(allSINGLE)[1]){
        if(allSINGLE[i,1]!=chr){
          chr=allSINGLE[i,1]
          suball=all[which(all[,1]==chr),]
        }
        # targets inside the interval
        sel=which(suball[,1]==allSINGLE[i,1] & as.numeric(suball[,2])>=as.numeric(allSINGLE[i,2]) & as.numeric(suball[,3])<=as.numeric(allSINGLE[i,3]))
        nsel=length(sel)
        nsub=dim(suball)[1]
        if(nsel==1){
          l1=max(sel[1]-11,1)
          l2=max(sel[1]-2,1)
          l3=sel[1]
          l4=sel[nsel]
          l5=min(sel[nsel]+2,nsub)
          l6=min(sel[nsel]+11,nsub)
          p1=sum(abs(as.numeric(suball[l1:l2,14])-1))
          m=mean(pmax(as.numeric(suball[l3:l4,14]),0))
          p2=sum(abs(pmax(as.numeric(suball[l3:l4,14]),0)-m))
          p3=sum(abs(as.numeric(suball[l5:l6,14])-1))
          allSINGLE[i,21]=p1/length(l1:l2)/2+p2/length(l3:l4)+p3/length(l5:l6)/2
        }
      }
    }

    # removing single targets that are included in larger CNVs
    take=allSINGLE[,1]
    take[1:length(take)]="Yes"
    chr="chr1"
    suballSINGLE=allSINGLE[which(allSINGLE[,1]==chr),]
    for (i in 1:dim(allMULTI)[1]){
      if(allMULTI[i,1]!=chr){
        chr=allMULTI[i,1]
        suballSINGLE=allSINGLE[which(allSINGLE[,1]==chr),]
      }
      sel=which(suballSINGLE[,1]==allMULTI[i,1] & ( (as.numeric(suballSINGLE[,2])>=as.numeric(allMULTI[i,2]) & as.numeric(suballSINGLE[,2])<=as.numeric(allMULTI[i,3])) | (as.numeric(suballSINGLE[,3])>=as.numeric(allMULTI[i,2]) & as.numeric(suballSINGLE[,3])<=as.numeric(allMULTI[i,3])) ) )
      take[which(allSINGLE[,1]==chr)[sel]]="No"
    }
    allSINGLE=allSINGLE[which(take=="Yes"),]

    # selecting multi target events with corrected Z-score > 2 or ratio and corrected ratio below 0.1
    allMULTI=allMULTI[which(abs((as.numeric(allMULTI[,10])-F1)/F2/as.numeric(allMULTI[,16]))>2 | (as.numeric(allMULTI[,14])<0.1 & as.numeric(allMULTI[,9])<0.1)),]

    # merging single and multiple target events
    allBOTH=rbind(allMULTI,allSINGLE)

    # correct Z-score and sorting by abs(z-score)
    allBOTH[,10]=(as.numeric(allBOTH[,10])-F1)/F2
    allBOTH=allBOTH[sort(abs(as.numeric(allBOTH[,10])),index.return=T,decreasing=T)$ix,]

    # annotating rank
    if(dim(allBOTH)[1]>0){
      allBOTH[,20]=1:length(allBOTH[,20])# rank 
    }

    if(dim(allBOTH)[1]>1000){
      allBOTH=allBOTH[1:1000,]
    }

    BOTHsave=rbind(BOTHsave,allBOTH)

  }

  ######

  print("  Step 7: merging of consecutive CNVs for targets only")

  {

    all=all[which(grepl("Off-target",all[,4])==F),]
    save(all,file=paste(folder,"/05_RData-files/data-targets-",pat,".RData",sep=""))

    # cor-counts cor-average   cor-sd cor-ratio cor-p-val rank
    data=all[,c(11,12,13,14,15,16)]
    data[,4]=((as.numeric(data[,1])-as.numeric(data[,2]))/as.numeric(data[,3])) # z-score from normal
    data[,5]=((as.numeric(data[,1])-as.numeric(data[,2])/2)/as.numeric(data[,3])) # z-score from het. deletion
    data[,6]=((as.numeric(data[,1])-as.numeric(data[,2])*1.5)/as.numeric(data[,3])) # z-score from het. duplication
    data[,4]=log10(2*pnorm(-abs(data[,4]))) # log10 p-value from z-score
    data[,5]=log10(2*pnorm(-abs(data[,5]))) # log10 p-value from z-score
    data[,6]=log10(2*pnorm(-abs(data[,6]))) # log10 p-value from z-score

    data=as.data.frame(cbind(data[,4],data[,5],data[,6]))
    colnames(data)=c("cn2","cn1","cn3")

    # transition probabilities for the merging
    p1=0.00005
    p2=0.5
    trans=matrix(c(1-2*p1,p1,p1, 1-p2,p2,0 ,1-p2,0,p2),nrow=3,byrow=T)

    events=matrix(ncol=5,nrow=0)
    events=as.data.frame(events)

    # loop over chromosomes and detect consecutive CNVs with viterbi algorithm
    for(chr in unique(all[,1])){
      data2=data[which(all[,1]==chr),]
      posi=as.numeric(all[which(all[,1]==chr),2])

      len=50000

      data2[which(data2[,1]==-Inf),1]=(-10)
      data2[which(data2[,2]==-Inf),2]=(-10)
      data2[which(data2[,3]==-Inf),3]=(-10)

      my.calls <- viterbi.hmm(as.matrix(trans),as.matrix(data2),as.integer(posi),len)

      if(dim(my.calls$calls)[1]>0){
        temp=matrix(ncol=5,nrow=dim(my.calls$calls)[1])
        temp=as.data.frame(temp)
        temp[,1]=chr # chromosome
        temp[,2]=posi[my.calls$calls[,1]] # begin
        temp[,3]=posi[my.calls$calls[,2]] # end
        temp[,4]=my.calls$calls[,3] # type 1=del, 2=dup
        temp[,5]=my.calls$calls[,4] # number of targets involved
        events=rbind(events,temp)
      }
    }
    colnames(events)=c("chr","begin","end","state","nbtargets")

    # removing single target events
    events=events[which(events[,5]>1),]
    nCNV=dim(events)[1]
    allMULTI=all[1:nCNV,]

    # merging consecutive CNVs
    if(nCNV>0){
      chr="chr1"
      suball=all[which(all[,1]==chr),]
      for (i in 1:nCNV){
        if(events[i,1]!=chr){
          chr=events[i,1]
          suball=all[which(all[,1]==chr),]
        }
        # targets inside the interval
        sel=which(suball[,1]==events[i,1] & as.numeric(suball[,2])>=as.numeric(events[i,2]) & as.numeric(suball[,2])<=as.numeric(events[i,3]))
        allMULTI[i,1]=suball[sel[1],1] # chromosome
        allMULTI[i,2]=min(as.numeric(suball[sel,2])) # begin
        allMULTI[i,3]=max(as.numeric(suball[sel,3])) # end
        order=sort(as.numeric(suball[sel,2]),index.return=T)$ix
        allMULTI[i,4]=paste(suball[sel[order],4],sep=",",collapse=";") # name
        allMULTI[i,5]=suball[sel[1],5] # sample ID
        allMULTI[i,6]=sum(as.numeric(suball[sel,6])) # raw-count
        allMULTI[i,7]=sum(as.numeric(suball[sel,7])) # raw-average
        allMULTI[i,8]=round(sqrt(mean((as.numeric(suball[sel,8])^2))),digits=2) # raw-sd
        allMULTI[i,9]=round(as.numeric(allMULTI[i,6])/as.numeric(allMULTI[i,7]),digits=3) # raw-ratio
        allMULTI[i,11]=sum(as.numeric(suball[sel,11])) # cor-counts
        allMULTI[i,12]=sum(as.numeric(suball[sel,12])) # cor-average
        allMULTI[i,13]=round(sqrt(mean((as.numeric(suball[sel,13])^2))),digits=2) # cor-sd
        allMULTI[i,10]=round((as.numeric(allMULTI[i,11])-as.numeric(allMULTI[i,12]))/as.numeric(allMULTI[i,13]),digits=2) # z-score
        allMULTI[i,14]=round(as.numeric(allMULTI[i,11])/as.numeric(allMULTI[i,12]),digits=3) # cor-ratio
        allMULTI[i,15]=signif(2*pnorm(-abs(as.numeric(allMULTI[i,10]))),digits=2) # p-value (z-score)
        allMULTI[i,16]=length(sel) # number target
        allMULTI[i,17]=suball[sel[1],17] # number of PCs removed
        allMULTI[i,18]=suball[sel[1],18] # nb-z6
        allMULTI[i,19]=suball[sel[1],19] # nb-samples

        # CQ score computation
        nsel=length(sel)
        nsub=dim(suball)[1]
        if(nsel>=3){
          l1=max(sel[1]-11,1)
          l2=max(sel[1]-2,1)
          l3=sel[2]
          l4=sel[nsel]-1
          l5=min(sel[nsel]+2,nsub)
          l6=min(sel[nsel]+11,nsub)
          p1=sum(abs(as.numeric(suball[l1:l2,14])-1))
          m=mean(pmax(as.numeric(suball[l3:l4,14]),0))
          p2=sum(abs(pmax(as.numeric(suball[l3:l4,14]),0)-m))
          p3=sum(abs(as.numeric(suball[l5:l6,14])-1))
          allMULTI[i,21]=p1/length(l1:l2)/2+p2/length(l3:l4)+p3/length(l5:l6)/2
        }
        if(nsel<=3){
          l1=max(sel[1]-11,1)
          l2=max(sel[1]-2,1)
          l3=sel[1]
          l4=sel[nsel]
          l5=min(sel[nsel]+2,nsub)
          l6=min(sel[nsel]+11,nsub)
          p1=sum(abs(as.numeric(suball[l1:l2,14])-1))
          m=mean(pmax(as.numeric(suball[l3:l4,14]),0))
          p2=sum(abs(pmax(as.numeric(suball[l3:l4,14]),0)-m))
          p3=sum(abs(as.numeric(suball[l5:l6,14])-1))
          allMULTI[i,21]=p1/length(l1:l2)/2+p2/length(l3:l4)+p3/length(l5:l6)/2
        }
      }
    }

    # taking only corrected Z-score above 4
    allSINGLE=all[which(abs(as.numeric(all[,10]))>((minZ*F2)+F1) | (as.numeric(all[,14])<0.1 & as.numeric(all[,9])<0.1)),]
    if(dim(allSINGLE)[1]>0){allSINGLE[,16]=1} # number of targets = 1 for isolated

    # CQ score for unitargets
    if(dim(allSINGLE)[1]>0){
      chr="chr1"
      suball=all[which(all[,1]==chr),]
      for (i in 1:dim(allSINGLE)[1]){
        if(allSINGLE[i,1]!=chr){
          chr=allSINGLE[i,1]
          suball=all[which(all[,1]==chr),]
        }
        # targets inside the interval
        sel=which(suball[,1]==allSINGLE[i,1] & as.numeric(suball[,2])>=as.numeric(allSINGLE[i,2]) & as.numeric(suball[,3])<=as.numeric(allSINGLE[i,3]))
        nsel=length(sel)
        nsub=dim(suball)[1]
        if(nsel==1){
          l1=max(sel[1]-11,1)
          l2=max(sel[1]-2,1)
          l3=sel[1]
          l4=sel[nsel]
          l5=min(sel[nsel]+2,nsub)
          l6=min(sel[nsel]+11,nsub)
          p1=sum(abs(as.numeric(suball[l1:l2,14])-1))
          m=mean(pmax(as.numeric(suball[l3:l4,14]),0))
          p2=sum(abs(pmax(as.numeric(suball[l3:l4,14]),0)-m))
          p3=sum(abs(as.numeric(suball[l5:l6,14])-1))
          allSINGLE[i,21]=p1/length(l1:l2)/2+p2/length(l3:l4)+p3/length(l5:l6)/2
        }
      }
    }

    # removing single targets that are included in larger CNVs
    take=allSINGLE[,1]
    take[1:length(take)]="Yes"
    chr="chr1"
    suballSINGLE=allSINGLE[which(allSINGLE[,1]==chr),]
    for (i in 1:dim(allMULTI)[1]){
      if(allMULTI[i,1]!=chr){
        chr=allMULTI[i,1]
        suballSINGLE=allSINGLE[which(allSINGLE[,1]==chr),]
      }
      sel=which(suballSINGLE[,1]==allMULTI[i,1] & ( (as.numeric(suballSINGLE[,2])>=as.numeric(allMULTI[i,2]) & as.numeric(suballSINGLE[,2])<=as.numeric(allMULTI[i,3])) | (as.numeric(suballSINGLE[,3])>=as.numeric(allMULTI[i,2]) & as.numeric(suballSINGLE[,3])<=as.numeric(allMULTI[i,3])) ) )
      take[which(allSINGLE[,1]==chr)[sel]]="No"
    }
    allSINGLE=allSINGLE[which(take=="Yes"),]

    # keeping only corrected Z-score above 3 for single targets
    allMULTI=allMULTI[which(abs((as.numeric(allMULTI[,10])-F1)/F2/as.numeric(allMULTI[,16]))>2 | (as.numeric(allMULTI[,14])<0.1 & as.numeric(allMULTI[,9])<0.1)),]
    
    # merging single and multiple target event
    allBOTH=rbind(allMULTI,allSINGLE)

    # correct Z-score and sorting by abs(z-score)
    allBOTH[,10]=(as.numeric(allBOTH[,10])-F1)/F2
    allBOTH=allBOTH[sort(abs(as.numeric(allBOTH[,10])),index.return=T,decreasing=T)$ix,]

    # annotating rank
    if(dim(allBOTH)[1]>0){
      allBOTH[,20]=1:length(allBOTH[,20])# rank 
    }

    # taking only first 1000 if more (to avoid long computation for poor quality samples)
    if(dim(allBOTH)[1]>1000){
      allBOTH=allBOTH[1:1000,]
    }

    TARsave=rbind(TARsave,allBOTH)

  }

}


print("Writing outputs and plots for all targets")
save(BOTHsave,file=paste(folder,"/05_RData-files/data-BOTHsave.RData",sep=""))

if(dim(BOTHsave)[1]>0){

  BOTHsave=cbind(BOTHsave,BOTHsave[,1:10],BOTHsave[,21],BOTHsave[,17:21],BOTHsave[,17:19])
  colnames(BOTHsave)[21:40]=c("Nb_overlapping_samples","Genes","Type","Ploidy","PQ","gnomAD-CNV_1%","ClinVar-patho_included","ncRNAs","gnomAD-CNV_ALL","ClinVar-patho_overlap","Exons","CQ","QUAL","Begin-min","End-max","Genes-possible","Exons-possible","Nb-exon-before-chr","Nb-exon-after-chr","RefSeq_Functional_Element")

  # annotating nb of samples with overlapping (80%) CNVs with HQ-sample
  for (i in 1:dim(BOTHsave)[1]){
    med=median(as.numeric(unique(BOTHsave[,c(5,18)])[,2]))
    n=which(is.element(BOTHsave[,5],HQ)==TRUE & BOTHsave[,5]!=BOTHsave[i,5] & BOTHsave[,1]==BOTHsave[i,1] & sign(as.numeric(BOTHsave[,10]))==sign(as.numeric(BOTHsave[i,10])) & ((as.numeric(BOTHsave[,2])>=as.numeric(BOTHsave[i,2]) & as.numeric(BOTHsave[,2])<=as.numeric(BOTHsave[i,3])) | (as.numeric(BOTHsave[,3])>=as.numeric(BOTHsave[i,2]) & as.numeric(BOTHsave[,3])<=as.numeric(BOTHsave[i,3])) | (as.numeric(BOTHsave[,2])<=as.numeric(BOTHsave[i,2]) & as.numeric(BOTHsave[,3])>=as.numeric(BOTHsave[i,3])) | (as.numeric(BOTHsave[,2])>=as.numeric(BOTHsave[i,2]) & as.numeric(BOTHsave[,3])<=as.numeric(BOTHsave[i,3])) ))
    n2=c()
    if(length(n)>0){
      for(j in 1:length(n)){
        overlap=(min(as.numeric(BOTHsave[i,3]),as.numeric(BOTHsave[n[j],3]))-max(as.numeric(BOTHsave[i,2]),as.numeric(BOTHsave[n[j],2])))/(as.numeric(BOTHsave[i,3])-as.numeric(BOTHsave[i,2]))
        if(overlap>0.8){
          n2=c(n2,n[j])
        }
      }
    }
    n1=length(unique(BOTHsave[n2,5]))
    BOTHsave[i,21]=max(n1,0)
  }

  # annotating genes and exons
  for (pat in colnames(dataALL)[5:num]){
    load(file=paste(folder,"/05_RData-files/data-",pat,".RData",sep=""))
    for (chr in unique(BOTHsave[which(BOTHsave[,5]==pat),1])){
      sel=which(BOTHsave[,5]==pat & BOTHsave[,1]==chr)
      allchr=all[which(all[,1]==chr),]
      exonschr=exons[which(exons[,1]==chr),]
      if(length(sel)>0){
        temp=BOTHsave
        BOTHsave=BOTHsave[sel,]
        for (i in 1:length(sel)){
          m1=which(allchr[,1]==BOTHsave[i,1] & as.numeric(allchr[,3])<=as.numeric(BOTHsave[i,2]))
          m2=which(allchr[,1]==BOTHsave[i,1] & as.numeric(allchr[,2])>=as.numeric(BOTHsave[i,3]))
          if(length(m1)>0){l1=max(as.numeric(allchr[m1,3]))} else {l1=as.numeric(BOTHsave[i,2])}
          if(length(m2)>0){l2=min(as.numeric(allchr[m2,2]))} else {l2=as.numeric(BOTHsave[i,3])}
          k1=exonschr[which(exonschr[,1]==BOTHsave[i,1] & (  (as.numeric(exonschr[,2])>=l1 & as.numeric(exonschr[,2])<=l2) | (as.numeric(exonschr[,3])>=l1 & as.numeric(exonschr[,3])<=l2) | (as.numeric(exonschr[,2])<=l1 & as.numeric(exonschr[,3])>=l2)  )),4]
          k2=strsplit(k1,"_")
          k3=unique(unlist(k2))
          k4=k3[-which(grepl("exon",k3)==T | grepl("chr",k3)==T | grepl("\\.",k3)==T | grepl("Off-target",k3)==T | k3=="NM")]
          k5=paste(k4,collapse=",")
          k5[which(k5=="")]="-"
          BOTHsave[i,36]=k5
          k6=paste(k1,collapse=";")
          BOTHsave[i,37]=k6
          BOTHsave[i,34]=l1
          BOTHsave[i,35]=l2

          k01=which(exonschr[,1]==BOTHsave[i,1] & as.numeric(exonschr[,2])>=as.numeric(BOTHsave[i,2]) & as.numeric(exonschr[,2])<=as.numeric(BOTHsave[i,3]))
          k02=which(exonschr[,1]==BOTHsave[i,1] & as.numeric(exonschr[,3])>=as.numeric(BOTHsave[i,2]) & as.numeric(exonschr[,3])<=as.numeric(BOTHsave[i,3]))
          k03=which(exonschr[,1]==BOTHsave[i,1] & as.numeric(exonschr[,2])<=as.numeric(BOTHsave[i,2]) & as.numeric(exonschr[,3])>=as.numeric(BOTHsave[i,3]))
          k1=exonschr[unique(c(k01,k02,k03)),4]
          k2=strsplit(k1,"_")
          k3=unique(unlist(k2))
          k4=k3[-which(grepl("exon",k3)==T | grepl("chr",k3)==T | grepl("\\.",k3)==T | grepl("Off-target",k3)==T | k3=="NM")]
          k5=paste(k4,collapse=",")
          k5[which(k5=="")]="-"
          BOTHsave[i,22]=k5
          k6=paste(k1,collapse=";")
          if(k6==""){k6="-"}
          BOTHsave[i,31]=k6
        }
        temp[sel,]=BOTHsave
        BOTHsave=temp
      }
    }
  }

  # annotation of ploidy and quality
  for (i in 1:dim(BOTHsave)[1]){

    c0=(as.numeric(BOTHsave[i,11])-0*as.numeric(BOTHsave[i,12]))/as.numeric(BOTHsave[i,13])
    c1=(as.numeric(BOTHsave[i,11])-0.5*as.numeric(BOTHsave[i,12]))/as.numeric(BOTHsave[i,13])
    c2=(as.numeric(BOTHsave[i,11])-1*as.numeric(BOTHsave[i,12]))/as.numeric(BOTHsave[i,13])
    c3=(as.numeric(BOTHsave[i,11])-1.5*as.numeric(BOTHsave[i,12]))/as.numeric(BOTHsave[i,13])
    c4=(as.numeric(BOTHsave[i,11])-2*as.numeric(BOTHsave[i,12]))/as.numeric(BOTHsave[i,13])
    c5=(as.numeric(BOTHsave[i,11])-2.5*as.numeric(BOTHsave[i,12]))/as.numeric(BOTHsave[i,13])
    c6=abs(c(c0,c1,c2,c3,c4,c5))

    type=c("deletion","deletion","normal","duplication","duplication","duplication")
    ploidy=c(0,1,2,3,4,">4")
    BOTHsave[i,23]=type[which.min(c6)]
    BOTHsave[i,24]=ploidy[which.min(c6)]
    BOTHsave[i,25]=round(abs(c2)-min(c6),digits=1)/as.numeric(BOTHsave[i,16])

  }
  BOTHsave=BOTHsave[which(BOTHsave[,24]!=2),]

  # annotation of nb-CNV-HQ
  for (sample in unique(BOTHsave[,5])){
    n=length(which(BOTHsave[,5]==sample & abs(as.numeric(BOTHsave[,10]))>6 & as.numeric(BOTHsave[,21])==0 & BOTHsave[,24]!="2" & as.numeric(BOTHsave[,25])>2))
    BOTHsave[which(BOTHsave[,5]==sample),18]=n
    sampleinfo[which(sampleinfo[,1]==sample),3]=n
  }
  write.table(sampleinfo,file=paste(folder,"/03_Samples-info/Samples_quality_info.tsv",sep=""),quote=F,sep="\t",row.names=F)

  # annotation of frequent gnomAD-SV-1%
  for (i in 1:dim(BOTHsave)[1]){
    j1=which(gnomAD[,1]==BOTHsave[i,1] & as.numeric(gnomAD[,2])<as.numeric(BOTHsave[i,2]) & as.numeric(gnomAD[,3])>as.numeric(BOTHsave[i,3]))
    if(length(j1)>0){
        BOTHsave[i,26]=paste(gnomAD[j1,4],collapse=",")
      } else {
        BOTHsave[i,26]="-"
      }
  }

  # annotation of ClinVar Pathogenic
  del=grepl("deletion",ClinVar[,4])
  dup=grepl("duplication",ClinVar[,4])
  for (i in 1:dim(BOTHsave)[1]){
    type=BOTHsave[i,23]
    if(BOTHsave[i,23]=="duplication"){
      j1=which(ClinVar[,1]==BOTHsave[i,1] & as.numeric(ClinVar[,2])>as.numeric(BOTHsave[i,2]) & as.numeric(ClinVar[,3])<as.numeric(BOTHsave[i,3]) & dup==T)
      j2=which(ClinVar[,1]==BOTHsave[i,1] & ((as.numeric(ClinVar[,2])>as.numeric(BOTHsave[i,2]) & as.numeric(ClinVar[,3])<as.numeric(BOTHsave[i,3])) | (as.numeric(ClinVar[,2])<as.numeric(BOTHsave[i,2]) & as.numeric(ClinVar[,3])>as.numeric(BOTHsave[i,2]))  | (as.numeric(ClinVar[,2])<as.numeric(BOTHsave[i,3]) & as.numeric(ClinVar[,3])>as.numeric(BOTHsave[i,3])) | (as.numeric(ClinVar[,2])<as.numeric(BOTHsave[i,2]) & as.numeric(ClinVar[,3])>as.numeric(BOTHsave[i,3]))) & dup==T)
    }
    if(BOTHsave[i,23]=="deletion"){
      j1=which(ClinVar[,1]==BOTHsave[i,1] & as.numeric(ClinVar[,2])>as.numeric(BOTHsave[i,2]) & as.numeric(ClinVar[,3])<as.numeric(BOTHsave[i,3]) & del==T)
      j2=which(ClinVar[,1]==BOTHsave[i,1] & ((as.numeric(ClinVar[,2])>as.numeric(BOTHsave[i,2]) & as.numeric(ClinVar[,3])<as.numeric(BOTHsave[i,3])) | (as.numeric(ClinVar[,2])<as.numeric(BOTHsave[i,2]) & as.numeric(ClinVar[,3])>as.numeric(BOTHsave[i,2]))  | (as.numeric(ClinVar[,2])<as.numeric(BOTHsave[i,3]) & as.numeric(ClinVar[,3])>as.numeric(BOTHsave[i,3])) | (as.numeric(ClinVar[,2])<as.numeric(BOTHsave[i,2]) & as.numeric(ClinVar[,3])>as.numeric(BOTHsave[i,3]))) & del==T)
    }
    if(length(j1)>0){
        BOTHsave[i,27]=paste(ClinVar[j1,4],collapse=",")
    } else {
        BOTHsave[i,27]="-"
    }
    if(length(j2)>0){
        BOTHsave[i,30]=paste(ClinVar[j2,4],collapse=",")
    } else {
        BOTHsave[i,30]="-"
    }
  }

  # annotation of ncRNA
  for (i in 1:dim(BOTHsave)[1]){
    j1=which(ncRNA[,1]==BOTHsave[i,1] & as.numeric(ncRNA[,2])>as.numeric(BOTHsave[i,2]) & as.numeric(ncRNA[,3])<as.numeric(BOTHsave[i,3]))
    j2=which(ncRNA[,1]==BOTHsave[i,1] & as.numeric(ncRNA[,2])<as.numeric(BOTHsave[i,2]) & as.numeric(ncRNA[,3])>as.numeric(BOTHsave[i,2]))
    j3=which(ncRNA[,1]==BOTHsave[i,1] & as.numeric(ncRNA[,2])<as.numeric(BOTHsave[i,3]) & as.numeric(ncRNA[,3])>as.numeric(BOTHsave[i,3]))
    j4=which(ncRNA[,1]==BOTHsave[i,1] & as.numeric(ncRNA[,2])<as.numeric(BOTHsave[i,2]) & as.numeric(ncRNA[,3])>as.numeric(BOTHsave[i,3]))
    j5=c(j1,j2,j3,j4)
    if(length(j5)>0){
        BOTHsave[i,28]=paste(ncRNA[j5,4],collapse=",")
      } else {
        BOTHsave[i,28]="-"
      }
  }

  # annotation of all gnomAD-CNV
  for (i in 1:dim(BOTHsave)[1]){
    j1=which(gnomAD2[,1]==BOTHsave[i,1] & as.numeric(gnomAD2[,2])<as.numeric(BOTHsave[i,2]) & as.numeric(gnomAD2[,3])>as.numeric(BOTHsave[i,3]))
    if(length(j1)>0){
        BOTHsave[i,29]=paste(gnomAD2[j1,4],collapse=",")
      } else {
        BOTHsave[i,29]="-"
      }
  }

  # annotation of number of targets from chromosome ends
  ex=unique(exons[,1:3])
  for (i in 1:dim(BOTHsave)[1]){
    BOTHsave[i,38]=length(which(ex[,1]==BOTHsave[i,1] & as.numeric(ex[,2])<as.numeric(BOTHsave[i,2])))
    BOTHsave[i,39]=length(which(ex[,1]==BOTHsave[i,1] & as.numeric(ex[,3])>as.numeric(BOTHsave[i,3])))
  }

  # RefSeq Functional Element Features
  for (i in 1:dim(BOTHsave)[1]){
    j1=which(FEfeats[,1]==BOTHsave[i,1] & as.numeric(FEfeats[,2])>as.numeric(BOTHsave[i,2]) & as.numeric(FEfeats[,3])<as.numeric(BOTHsave[i,3]))
    j2=which(FEfeats[,1]==BOTHsave[i,1] & as.numeric(FEfeats[,2])<as.numeric(BOTHsave[i,2]) & as.numeric(FEfeats[,3])>as.numeric(BOTHsave[i,2]))
    j3=which(FEfeats[,1]==BOTHsave[i,1] & as.numeric(FEfeats[,2])<as.numeric(BOTHsave[i,3]) & as.numeric(FEfeats[,3])>as.numeric(BOTHsave[i,3]))
    j4=which(FEfeats[,1]==BOTHsave[i,1] & as.numeric(FEfeats[,2])<as.numeric(BOTHsave[i,2]) & as.numeric(FEfeats[,3])>as.numeric(BOTHsave[i,3]))
    j5=c(j1,j2,j3,j4)
    if(length(j5)>0){
        BOTHsave[i,40]=paste(FEfeats[j5,4],collapse=",")
      } else {
        BOTHsave[i,40]="-"
      }
  }

  # QUAL
  BOTHsave[,32]=round((1/as.numeric(BOTHsave[,32])),digits=2)
  BOTHsave[which(BOTHsave[,32]=="Inf"),32]=10000
  BOTHsave[,33]=round(abs(as.numeric(BOTHsave[,10]))/as.numeric(BOTHsave[,16]),digits=2)
  BOTHsave[,25]=round(as.numeric(BOTHsave[,25]),digits=2)

  # sorting by QUAL
  BOTHsave=BOTHsave[sort(abs(as.numeric(BOTHsave[,10])),index.return=T,decreasing=T)$ix,]

  # writing outputs
  cols=c(5,1:3,34:35,23,24,16,4,22,31,36,28,40,21,26,29,11:14,10,25,32,33,20,18,27,30,38,39)
  sel=which(BOTHsave[,23]!="normal")
  write.table(BOTHsave[sel,cols],file=paste(folder,"/04_CNVs-results/CNVs-all.tsv",sep=""),quote=F,sep="\t",row.names=F)
  sel=which(abs(as.numeric(BOTHsave[,10]))>5 & BOTHsave[,24]!="2" & as.numeric(BOTHsave[,25])>2 & as.numeric(BOTHsave[,32])>2.5 & as.numeric(BOTHsave[,33])>3)
  write.table(BOTHsave[sel,cols],file=paste(folder,"/04_CNVs-results/CNVs-all.HQ.tsv",sep=""),quote=F,sep="\t",row.names=F)
  sel=which(abs(as.numeric(BOTHsave[,10]))>5 & BOTHsave[,24]!="2" & as.numeric(BOTHsave[,25])>2 & as.numeric(BOTHsave[,32])>2.5 & as.numeric(BOTHsave[,33])>3 & as.numeric(BOTHsave[,21])==0)
  write.table(BOTHsave[sel,cols],file=paste(folder,"/04_CNVs-results/CNVs-all.HQ.unique.tsv",sep=""),quote=F,sep="\t",row.names=F)

  varav=mean(as.numeric(sampleinfo[,6]))
  varsd=sd(as.numeric(sampleinfo[,6]))
  lowqual=which(abs((as.numeric(sampleinfo[,6])-varav)/varsd)>2)
  if(length(lowqual)>0){
      good=as.character(sampleinfo[-lowqual,1])
  } else {
      good=as.character(sampleinfo[,1])
  }
  sel=which(abs(as.numeric(BOTHsave[,10]))>5 & BOTHsave[,24]!="2" & as.numeric(BOTHsave[,25])>2 & as.numeric(BOTHsave[,32])>2.5 & as.numeric(BOTHsave[,33])>3 & as.numeric(BOTHsave[,21])==0 & is.element(BOTHsave[,5],good))
  write.table(BOTHsave[sel,cols],file=paste(folder,"/04_CNVs-results/CNVs-all.HQ.unique.HQ-sample.tsv",sep=""),quote=F,sep="\t",row.names=F)
  sel=which(abs(as.numeric(BOTHsave[,10]))>5 & BOTHsave[,24]!="2" & as.numeric(BOTHsave[,25])>2 & as.numeric(BOTHsave[,32])>2.5 & as.numeric(BOTHsave[,33])>3 & is.element(BOTHsave[,5],good))
  write.table(BOTHsave[sel,cols],file=paste(folder,"/04_CNVs-results/CNVs-all.HQ.HQ-sample.tsv",sep=""),quote=F,sep="\t",row.names=F)
  sel=which(BOTHsave[,23]!="normal" & is.element(BOTHsave[,5],good))
  write.table(BOTHsave[sel,cols],file=paste(folder,"/04_CNVs-results/CNVs-all.HQ-sample.tsv",sep=""),quote=F,sep="\t",row.names=F)

  write.table(as.character(sampleinfo[lowqual,1]),file=paste(folder,"/03_Samples-info/Samples_low-quality.tsv",sep=""),quote=F,sep="\t",row.names=F,col.names=F)

  # BED file for IGV
  bed=BOTHsave[,c(1:3,5,23)]
  for(i in 1:dim(bed)[1]){
    bed[i,4]=paste(bed[i,4],"_",bed[i,5],sep="")
  }
  write.table(bed[,1:4],file=paste(folder,"/02_BED-files/CNVs-all-IGV.bed",sep=""),quote=F,sep="\t",row.names=F,col.names=F)

  # BED file for AnnotSV
  bed=BOTHsave[,c(1:3,5,23)]
  bed[which(bed[,5]=="duplication"),5]="DUP"
  bed[which(bed[,5]=="deletion"),5]="DEL"
  write.table(bed[which(bed[,5]=="DUP" | bed[,5]=="DEL"),c(1,2,3,5,4)],file=paste(folder,"/02_BED-files/CNVs-all-AnnotSV.bed",sep=""),quote=F,sep="\t",row.names=F,col.names=F)

  # writing individual outputs and doing plots for best CNVs
  for (pat in colnames(dataALL)[5:num]){
    write.table(BOTHsave[which(BOTHsave[,5]==pat),cols],file=paste(folder,"/",pat,"/",pat,"-merged.tsv",sep=""),quote=F,sep="\t",row.names=F)
    allBOTH=BOTHsave[which(BOTHsave[,5]==pat),]
    allBOTH=allBOTH[which(abs(as.numeric(allBOTH[,10]))>5),]
    allBOTH=allBOTH[sort(abs(as.numeric(allBOTH[,10])),index.return=T,decreasing=T)$ix,]

    #plotting best CNVs
    if (dim(allBOTH)[1]>0){
      nb=length(which(abs(as.numeric(allBOTH[,10]))>5))
      if (nb>0){
        for (i in 1:min(dim(allBOTH)[1],nb,nbPlots)){
          if(as.numeric(allBOTH[i,16])<1000){
            chr=allBOTH[i,1]
            begin=as.numeric(allBOTH[i,2])
            end=as.numeric(allBOTH[i,3])
            data=paste(folder,"/05_RData-files/data-plot-",pat,".RData",sep="")
            side=20
            z=allBOTH[i,10]

            pdf=paste(folder,"/",pat,"/plots_on-targets_off-targets/",pat,"-",chr,"-",begin,"-",end,".HQ.pdf",sep="")
            plotcnv(chr,begin,end,pat,data,pdf,side,TRUE,FALSE)

            pdf=paste(folder,"/",pat,"/plots_on-targets_off-targets/",pat,"-",chr,"-",begin,"-",end,".all.pdf",sep="")
            plotcnvall(chr,begin,end,pat,data,pdf,side,TRUE,FALSE)
          }
        }
      }
    }

  }

}



print("Writing outputs and plots for on-targets only")
save(TARsave,file=paste(folder,"/05_RData-files/data-TARsave.RData",sep=""))


if(dim(TARsave)[1]>0){

  TARsave=cbind(TARsave,TARsave[,1:10],TARsave[,21],TARsave[,17:21],TARsave[,17:19])
  colnames(TARsave)[21:40]=c("Nb_overlapping_samples","Genes","Type","Ploidy","Ploidy-qual","gnomAD-CNV_1%","ClinVar-patho_included","ncRNAs","gnomAD-CNV_ALL","ClinVar-patho_overlap","Exons","CNV-qual","QUAL","Begin-min","End-max","Genes-possible","Exons-possible","Nb-exon-before-chr","Nb-exon-after-chr","RefSeq_Functional_Element")

  # annotating nb of samples with overlapping (80%) CNVs with HQ-sample
  for (i in 1:dim(TARsave)[1]){
    med=median(as.numeric(unique(TARsave[,c(5,18)])[,2]))
    n=which(is.element(TARsave[,5],HQ)==TRUE & TARsave[,5]!=TARsave[i,5] & TARsave[,1]==TARsave[i,1] & sign(as.numeric(TARsave[,10]))==sign(as.numeric(TARsave[i,10])) & ((as.numeric(TARsave[,2])>=as.numeric(TARsave[i,2]) & as.numeric(TARsave[,2])<=as.numeric(TARsave[i,3])) | (as.numeric(TARsave[,3])>=as.numeric(TARsave[i,2]) & as.numeric(TARsave[,3])<=as.numeric(TARsave[i,3])) | (as.numeric(TARsave[,2])<=as.numeric(TARsave[i,2]) & as.numeric(TARsave[,3])>=as.numeric(TARsave[i,3])) | (as.numeric(TARsave[,2])>=as.numeric(TARsave[i,2]) & as.numeric(TARsave[,3])<=as.numeric(TARsave[i,3])) ))
    n2=c()
    if(length(n)>0){
      for(j in 1:length(n)){
        overlap=(min(as.numeric(TARsave[i,3]),as.numeric(TARsave[n[j],3]))-max(as.numeric(TARsave[i,2]),as.numeric(TARsave[n[j],2])))/(as.numeric(TARsave[i,3])-as.numeric(TARsave[i,2]))
        #overlap2=(min(as.numeric(TARsave[i,3]),as.numeric(TARsave[n[j],3]))-max(as.numeric(TARsave[i,2]),as.numeric(TARsave[n[j],2])))/(as.numeric(TARsave[n[j],3])-as.numeric(TARsave[n[j],2]))
        if(overlap>0.8){
          n2=c(n2,n[j])
        }
      }
    }
    n1=length(unique(TARsave[n2,5]))
    TARsave[i,21]=max(n1,0)
  }

  # annotating genes and exons
  for (pat in colnames(dataALL)[5:num]){
    load(file=paste(folder,"/05_RData-files/data-targets-",pat,".RData",sep=""))
    for (chr in unique(TARsave[which(TARsave[,5]==pat),1])){
      sel=which(TARsave[,5]==pat & TARsave[,1]==chr)
      allchr=all[which(all[,1]==chr),]
      exonschr=exons[which(exons[,1]==chr),]
      if(length(sel)>0){
        temp=TARsave
        TARsave=TARsave[sel,]
        for (i in 1:length(sel)){
          m1=which(allchr[,1]==TARsave[i,1] & as.numeric(allchr[,3])<=as.numeric(TARsave[i,2]))
          m2=which(allchr[,1]==TARsave[i,1] & as.numeric(allchr[,2])>=as.numeric(TARsave[i,3]))
          if(length(m1)>0){l1=max(as.numeric(allchr[m1,3]))} else {l1=as.numeric(TARsave[i,2])}
          if(length(m2)>0){l2=min(as.numeric(allchr[m2,2]))} else {l2=as.numeric(TARsave[i,3])}
          k1=exonschr[which(exonschr[,1]==TARsave[i,1] & (  (as.numeric(exonschr[,2])>=l1 & as.numeric(exonschr[,2])<=l2) | (as.numeric(exonschr[,3])>=l1 & as.numeric(exonschr[,3])<=l2) | (as.numeric(exonschr[,2])<=l1 & as.numeric(exonschr[,3])>=l2)  )),4]
          k2=strsplit(k1,"_")
          k3=unique(unlist(k2))
          k4=k3[-which(grepl("exon",k3)==T | grepl("chr",k3)==T | grepl("\\.",k3)==T | grepl("Off-target",k3)==T | k3=="NM")]
          k5=paste(k4,collapse=",")
          k5[which(k5=="")]="-"
          TARsave[i,36]=k5
          k6=paste(k1,collapse=";")
          TARsave[i,37]=k6
          TARsave[i,34]=l1
          TARsave[i,35]=l2

          k01=which(exonschr[,1]==TARsave[i,1] & as.numeric(exonschr[,2])>=as.numeric(TARsave[i,2]) & as.numeric(exonschr[,2])<=as.numeric(TARsave[i,3]))
          k02=which(exonschr[,1]==TARsave[i,1] & as.numeric(exonschr[,3])>=as.numeric(TARsave[i,2]) & as.numeric(exonschr[,3])<=as.numeric(TARsave[i,3]))
          k03=which(exonschr[,1]==TARsave[i,1] & as.numeric(exonschr[,2])<=as.numeric(TARsave[i,2]) & as.numeric(exonschr[,3])>=as.numeric(TARsave[i,3]))
          k1=exonschr[unique(c(k01,k02,k03)),4]
          k2=strsplit(k1,"_")
          k3=unique(unlist(k2))
          k4=k3[-which(grepl("exon",k3)==T | grepl("chr",k3)==T | grepl("\\.",k3)==T | grepl("Off-target",k3)==T | k3=="NM")]
          k5=paste(k4,collapse=",")
          k5[which(k5=="")]="-"
          TARsave[i,22]=k5
          k6=paste(k1,collapse=";")
          if(k6==""){k6="-"}
          TARsave[i,31]=k6
        }
        temp[sel,]=TARsave
        TARsave=temp
      }
    }
  }

  # annotation of ploidy and quality
  for (i in 1:dim(TARsave)[1]){

    c0=(as.numeric(TARsave[i,11])-0*as.numeric(TARsave[i,12]))/as.numeric(TARsave[i,13])
    c1=(as.numeric(TARsave[i,11])-0.5*as.numeric(TARsave[i,12]))/as.numeric(TARsave[i,13])
    c2=(as.numeric(TARsave[i,11])-1*as.numeric(TARsave[i,12]))/as.numeric(TARsave[i,13])
    c3=(as.numeric(TARsave[i,11])-1.5*as.numeric(TARsave[i,12]))/as.numeric(TARsave[i,13])
    c4=(as.numeric(TARsave[i,11])-2*as.numeric(TARsave[i,12]))/as.numeric(TARsave[i,13])
    c5=(as.numeric(TARsave[i,11])-2.5*as.numeric(TARsave[i,12]))/as.numeric(TARsave[i,13])
    c6=abs(c(c0,c1,c2,c3,c4,c5))

    type=c("deletion","deletion","normal","duplication","duplication","duplication")
    ploidy=c(0,1,2,3,4,">4")
    TARsave[i,23]=type[which.min(c6)]
    TARsave[i,24]=ploidy[which.min(c6)]
    TARsave[i,25]=round(abs(c2)-min(c6),digits=1)/as.numeric(TARsave[i,16])

  }
  TARsave=TARsave[which(TARsave[,24]!=2),]

  # annotation of nb-CNV-HQ
  for (sample in unique(TARsave[,5])){
    n=length(which(TARsave[,5]==sample & abs(as.numeric(TARsave[,10]))>6 & as.numeric(TARsave[,21])==0 & as.numeric(TARsave[,24])!=2 & as.numeric(TARsave[,25])>2))
    TARsave[which(TARsave[,5]==sample),18]=n
  }

  # annotation of frequent gnomAD-SV-1%
  for (i in 1:dim(TARsave)[1]){
    j1=which(gnomAD[,1]==TARsave[i,1] & as.numeric(gnomAD[,2])<as.numeric(TARsave[i,2]) & as.numeric(gnomAD[,3])>as.numeric(TARsave[i,3]))
    if(length(j1)>0){
        TARsave[i,26]=paste(gnomAD[j1,4],collapse=",")
      } else {
        TARsave[i,26]="-"
      }
  }

  # annotation of ClinVar Pathogenic
  del=grepl("deletion",ClinVar[,4])
  dup=grepl("duplication",ClinVar[,4])
  for (i in 1:dim(TARsave)[1]){
    type=TARsave[i,23]
    if(TARsave[i,23]=="duplication"){
      j1=which(ClinVar[,1]==TARsave[i,1] & as.numeric(ClinVar[,2])>as.numeric(TARsave[i,2]) & as.numeric(ClinVar[,3])<as.numeric(TARsave[i,3]) & dup==T)
      j2=which(ClinVar[,1]==TARsave[i,1] & ((as.numeric(ClinVar[,2])>as.numeric(TARsave[i,2]) & as.numeric(ClinVar[,3])<as.numeric(TARsave[i,3])) | (as.numeric(ClinVar[,2])<as.numeric(TARsave[i,2]) & as.numeric(ClinVar[,3])>as.numeric(TARsave[i,2]))  | (as.numeric(ClinVar[,2])<as.numeric(TARsave[i,3]) & as.numeric(ClinVar[,3])>as.numeric(TARsave[i,3])) | (as.numeric(ClinVar[,2])<as.numeric(TARsave[i,2]) & as.numeric(ClinVar[,3])>as.numeric(TARsave[i,3]))) & dup==T)
    }
    if(TARsave[i,23]=="deletion"){
      j1=which(ClinVar[,1]==TARsave[i,1] & as.numeric(ClinVar[,2])>as.numeric(TARsave[i,2]) & as.numeric(ClinVar[,3])<as.numeric(TARsave[i,3]) & del==T)
      j2=which(ClinVar[,1]==TARsave[i,1] & ((as.numeric(ClinVar[,2])>as.numeric(TARsave[i,2]) & as.numeric(ClinVar[,3])<as.numeric(TARsave[i,3])) | (as.numeric(ClinVar[,2])<as.numeric(TARsave[i,2]) & as.numeric(ClinVar[,3])>as.numeric(TARsave[i,2]))  | (as.numeric(ClinVar[,2])<as.numeric(TARsave[i,3]) & as.numeric(ClinVar[,3])>as.numeric(TARsave[i,3])) | (as.numeric(ClinVar[,2])<as.numeric(TARsave[i,2]) & as.numeric(ClinVar[,3])>as.numeric(TARsave[i,3]))) & del==T)
    }
    if(length(j1)>0){
        TARsave[i,27]=paste(ClinVar[j1,4],collapse=",")
    } else {
        TARsave[i,27]="-"
    }
    if(length(j2)>0){
        TARsave[i,30]=paste(ClinVar[j2,4],collapse=",")
    } else {
        TARsave[i,30]="-"
    }
  }

  # annotation of ncRNA
  for (i in 1:dim(TARsave)[1]){
    j1=which(ncRNA[,1]==TARsave[i,1] & as.numeric(ncRNA[,2])>as.numeric(TARsave[i,2]) & as.numeric(ncRNA[,3])<as.numeric(TARsave[i,3]))
    j2=which(ncRNA[,1]==TARsave[i,1] & as.numeric(ncRNA[,2])<as.numeric(TARsave[i,2]) & as.numeric(ncRNA[,3])>as.numeric(TARsave[i,2]))
    j3=which(ncRNA[,1]==TARsave[i,1] & as.numeric(ncRNA[,2])<as.numeric(TARsave[i,3]) & as.numeric(ncRNA[,3])>as.numeric(TARsave[i,3]))
    j4=which(ncRNA[,1]==TARsave[i,1] & as.numeric(ncRNA[,2])<as.numeric(TARsave[i,2]) & as.numeric(ncRNA[,3])>as.numeric(TARsave[i,3]))
    j5=c(j1,j2,j3,j4)
    if(length(j5)>0){
        TARsave[i,28]=paste(ncRNA[j5,4],collapse=",")
      } else {
        TARsave[i,28]="-"
      }
  }

  # annotation of all gnomAD-CNV
  for (i in 1:dim(TARsave)[1]){
    j1=which(gnomAD2[,1]==TARsave[i,1] & as.numeric(gnomAD2[,2])<as.numeric(TARsave[i,2]) & as.numeric(gnomAD2[,3])>as.numeric(TARsave[i,3]))
    if(length(j1)>0){
        TARsave[i,29]=paste(gnomAD2[j1,4],collapse=",")
      } else {
        TARsave[i,29]="-"
      }
  }

  # annotation of number of targets from chromosome ends
  ex=unique(exons[,1:3])
  for (i in 1:dim(TARsave)[1]){
    TARsave[i,38]=length(which(ex[,1]==TARsave[i,1] & as.numeric(ex[,2])<as.numeric(TARsave[i,2])))
    TARsave[i,39]=length(which(ex[,1]==TARsave[i,1] & as.numeric(ex[,3])>as.numeric(TARsave[i,3])))
  }

  # RefSeq Functional Element Features
  for (i in 1:dim(TARsave)[1]){
    j1=which(FEfeats[,1]==TARsave[i,1] & as.numeric(FEfeats[,2])>as.numeric(TARsave[i,2]) & as.numeric(FEfeats[,3])<as.numeric(TARsave[i,3]))
    j2=which(FEfeats[,1]==TARsave[i,1] & as.numeric(FEfeats[,2])<as.numeric(TARsave[i,2]) & as.numeric(FEfeats[,3])>as.numeric(TARsave[i,2]))
    j3=which(FEfeats[,1]==TARsave[i,1] & as.numeric(FEfeats[,2])<as.numeric(TARsave[i,3]) & as.numeric(FEfeats[,3])>as.numeric(TARsave[i,3]))
    j4=which(FEfeats[,1]==TARsave[i,1] & as.numeric(FEfeats[,2])<as.numeric(TARsave[i,2]) & as.numeric(FEfeats[,3])>as.numeric(TARsave[i,3]))
    j5=c(j1,j2,j3,j4)
    if(length(j5)>0){
        TARsave[i,40]=paste(FEfeats[j5,4],collapse=",")
      } else {
        TARsave[i,40]="-"
      }
  }

  # QUAL
  TARsave[,32]=round((1/as.numeric(TARsave[,32])),digits=2)
  TARsave[which(TARsave[,32]=="Inf"),32]=10000
  TARsave[,33]=round(abs(as.numeric(TARsave[,10]))/as.numeric(TARsave[,16]),digits=2)
  TARsave[,25]=round(as.numeric(TARsave[,25]),digits=2)

  # sorting by QUAL
  TARsave=TARsave[sort(abs(as.numeric(TARsave[,10])),index.return=T,decreasing=T)$ix,]

  # writing outputs
  cols=c(5,1:3,34:35,23,24,16,4,22,31,36,28,40,21,26,29,11:14,10,25,32,33,20,18,27,30,38:39)
  sel=which(TARsave[,23]!="normal")
  write.table(TARsave[sel,cols],file=paste(folder,"/04_CNVs-results/CNVs-targets-only.tsv",sep=""),quote=F,sep="\t",row.names=F)
  sel=which(abs(as.numeric(TARsave[,10]))>5 & TARsave[,24]!="2" & as.numeric(TARsave[,25])>2 & as.numeric(TARsave[,32])>2.5 & as.numeric(TARsave[,33])>3)
  write.table(TARsave[sel,cols],file=paste(folder,"/04_CNVs-results/CNVs-targets-only.HQ.tsv",sep=""),quote=F,sep="\t",row.names=F)
  sel=which(abs(as.numeric(TARsave[,10]))>5 & TARsave[,24]!="2" & as.numeric(TARsave[,25])>2 & as.numeric(TARsave[,32])>2.5 & as.numeric(TARsave[,33])>3 & as.numeric(TARsave[,21])==0)
  write.table(TARsave[sel,cols],file=paste(folder,"/04_CNVs-results/CNVs-targets-only.HQ.unique.tsv",sep=""),quote=F,sep="\t",row.names=F)

  varav=mean(as.numeric(sampleinfo[,6]))
  varsd=sd(as.numeric(sampleinfo[,6]))
  lowqual=which(abs((as.numeric(sampleinfo[,6])-varav)/varsd)>2)
  if(length(lowqual)>0){
      good=as.character(sampleinfo[-lowqual,1])
  } else {
      good=as.character(sampleinfo[,1])
  }
  sel=which(abs(as.numeric(TARsave[,10]))>5 & TARsave[,24]!="2" & as.numeric(TARsave[,25])>2 & as.numeric(TARsave[,32])>2.5 & as.numeric(TARsave[,33])>3 & as.numeric(TARsave[,21])==0 & is.element(TARsave[,5],good))
  write.table(TARsave[sel,cols],file=paste(folder,"/04_CNVs-results/CNVs-targets-only.HQ.unique.HQ-sample.tsv",sep=""),quote=F,sep="\t",row.names=F)
  sel=which(abs(as.numeric(TARsave[,10]))>5 & TARsave[,24]!="2" & as.numeric(TARsave[,25])>2 & as.numeric(TARsave[,32])>2.5 & as.numeric(TARsave[,33])>3 & is.element(TARsave[,5],good))
  write.table(TARsave[sel,cols],file=paste(folder,"/04_CNVs-results/CNVs-targets-only.HQ.HQ-sample.tsv",sep=""),quote=F,sep="\t",row.names=F)
  sel=which(TARsave[,23]!="normal" & is.element(TARsave[,5],good))
  write.table(TARsave[sel,cols],file=paste(folder,"/04_CNVs-results/CNVs-targets-only.HQ-sample.tsv",sep=""),quote=F,sep="\t",row.names=F)

  # BED file for IGV
  bed=TARsave[,c(1:3,5,23)]
  for(i in 1:dim(bed)[1]){
    bed[i,4]=paste(bed[i,4],"_",bed[i,5],sep="")
  }
  write.table(bed[,1:4],file=paste(folder,"/02_BED-files/CNVs-targets-only-IGV.bed",sep=""),quote=F,sep="\t",row.names=F,col.names=F)

  # BED file for AnnotSV
  bed=TARsave[,c(1:3,5,23)]
  bed[which(bed[,5]=="duplication"),5]="DUP"
  bed[which(bed[,5]=="deletion"),5]="DEL"
  write.table(bed[which(bed[,5]=="DUP" | bed[,5]=="DEL"),c(1,2,3,5,4)],file=paste(folder,"/02_BED-files/CNVs-targets-only-AnnotSV.bed",sep=""),quote=F,sep="\t",row.names=F,col.names=F)


  # writing individual outputs and doing plots for best CNVs
  for (pat in colnames(dataALL)[5:num]){
     write.table(TARsave[which(TARsave[,5]==pat),cols],file=paste(folder,"/",pat,"/",pat,"-merged_targets-only.tsv",sep=""),quote=F,sep="\t",row.names=F)
    allBOTH=TARsave[which(TARsave[,5]==pat),]
    allBOTH=allBOTH[which(abs(as.numeric(allBOTH[,10]))>5),]
    allBOTH=allBOTH[sort(abs(as.numeric(allBOTH[,10])),index.return=T,decreasing=T)$ix,]

    #plotting best CNVs
    if (dim(allBOTH)[1]>0){
      nb=length(which(abs(as.numeric(allBOTH[,10]))>5))
      if (nb>0){
        for (i in 1:min(dim(allBOTH)[1],nb,nbPlots)){
          if(as.numeric(allBOTH[i,16])<1000){
            chr=allBOTH[i,1]
            begin=as.numeric(allBOTH[i,2])
            end=as.numeric(allBOTH[i,3])
            data=paste(folder,"/05_RData-files/data-plot-",pat,".RData",sep="")
            side=20
            z=allBOTH[i,10]

            pdf=paste(folder,"/",pat,"/plots_on-targets-only/",pat,"-",chr,"-",begin,"-",end,"_on-targets.HQ.pdf",sep="")
            plotcnv(chr,begin,end,pat,data,pdf,side,FALSE,FALSE)

            pdf=paste(folder,"/",pat,"/plots_on-targets-only/",pat,"-",chr,"-",begin,"-",end,"_on-targets.all.pdf",sep="")
            plotcnvall(chr,begin,end,pat,data,pdf,side,FALSE,FALSE)
          }
        }
      }
    }

  }

}

