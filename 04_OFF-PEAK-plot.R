
#!/usr/bin/env Rscript

library(optparse)

option_list = list(
  make_option(c("--ID"), type="character", default="NA", 
              help="ID of sample", metavar="character"),

  make_option(c("--chr"), type="character", default="NA", 
              help="chromosome", metavar="character"),

  make_option(c("--begin"), type="numeric", default=NA, 
              help="begining of CNV", metavar="numeric"),

  make_option(c("--end"), type="numeric", default=NA, 
              help="end of CNV", metavar="numeric"),

  make_option(c("--side"), type="numeric", default=20, 
              help="number of side targets to plot", metavar="numeric"),

  make_option(c("--batch"), type="character", default="NA", 
              help="RData file with data", metavar="character"),

  make_option(c("--out"), type="character", default="NA", 
              help="output folder", metavar="character"),

  make_option(c("--databasefile"), type="character", default="NA", 
              help="RData file for CNV annotation", metavar="character"),

  make_option(c("--UseCano"), type="character", default="FALSE", 
              help="If TRUE, only plot canonical isoforms", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
args = parse_args(opt_parser);

ID=as.character(args[1])
chr=as.character(args[2])
begin=as.numeric(args[3])
end=as.numeric(args[4])
side=as.numeric(args[5])
data=paste(as.character(args[6]),"/05_RData-files/data-plot-",ID,".RData",sep="")
out=as.character(args[7])
databasefile=as.character(args[8])
UseCano=as.character(args[9])

if(ID=="NA"){
  stop("You need to include the ID with the --ID option. Exit.")
}
if(chr=="NA"){
  stop("You need to include the chromosome with the --chr option. Exit.")
}
if(is.na(begin)){
  stop("You need to include the begin position with the --begin option. Exit.")
}
if(is.na(end)){
  stop("You need to include the end position with the --end option. Exit.")
}
if(data=="NA"){
  stop("You need to include the Rdata file for the sample with the --data option. Exit.")
}
if(out=="NA"){
  stop("You need to include the output directory with the --out option. Exit.")
}
if(databasefile=="NA"){
  stop("You need to include the Rdata file with information from databases with the --databasefile option. Exit.")
}

load(databasefile)

# creating output directories
dir.create(out, showWarnings = FALSE)
out1=paste(out,"plots_on-targets_off-targets",sep="/")
out2=paste(out,"plots_on-targets-only",sep="/")
dir.create(out1, showWarnings = FALSE)
dir.create(out2, showWarnings = FALSE)


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


#plotting targets and antitargets
print("Doing plot with on-targets and offtargets")
pdf=paste(out1,"/",ID,"-",chr,"-",begin,"-",end,".HQ.pdf",sep="")
plotcnv(chr,begin,end,ID,data,pdf,side,TRUE,UseCano)
pdf=paste(out1,"/",ID,"-",chr,"-",begin,"-",end,".all.pdf",sep="")
plotcnvall(chr,begin,end,ID,data,pdf,side,TRUE,UseCano)

#plotting targets and antitargets
print("Doing plot with on-targets only")
pdf=paste(out2,"/",ID,"-",chr,"-",begin,"-",end,"_on-targets.HQ.pdf",sep="")
plotcnv(chr,begin,end,ID,data,pdf,side,FALSE,UseCano)
pdf=paste(out2,"/",ID,"-",chr,"-",begin,"-",end,"_on-targets.all.pdf",sep="")
plotcnvall(chr,begin,end,ID,data,pdf,side,FALSE,UseCano)


