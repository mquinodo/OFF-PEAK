# Written by Mathieu Quinodoz, Basel, Switzerland

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

a=read.table(file=args[1],header=F)
b=read.table(file=args[2],header=F)

a=cbind(a,a[,4])
a[,5]="No"

for (i in 1:dim(a)[1]){
	if(i==1){
		temp=b[which(b[,1]==a[i,1]),]
		} else {
			if(a[i,1]!=a[i-1,1]){
				temp=b[which(b[,1]==a[i,1]),]
			}
		}
	temp2=temp[which(abs(temp[,2]-a[i,2])<100000),]

	done=0
	# completely covered
	overlap=which(temp2[,2]<=a[i,2] & temp2[,3]>=a[i,3])
	if(length(overlap)>0){
		overlap=overlap[which.max(temp2[overlap,3]-temp2[overlap,2])]
	}
	if(length(overlap)>0){
		a[i,5]="Yes"
		done=1
	}

	# covered only on the left
	overlap=which(temp2[,2]<=a[i,2] & temp2[,3]>=a[i,2] & temp2[,3]<a[i,3] & done==0) 
	if(length(overlap)>0){
		overlap=overlap[which.max(temp2[overlap,3]-temp2[overlap,2])]
	}
	if(length(overlap)>0){
		a[i,5]="Yes"
		done=1
	}

	# covered only on the right
	overlap=which(temp2[,2]<=a[i,3] & temp2[,3]>=a[i,3] & temp2[,2]>a[i,2] & done==0) 
	if(length(overlap)>0){
		overlap=overlap[which.max(temp2[overlap,3]-temp2[overlap,2])]
	}
	if(length(overlap)>0){
		a[i,5]="Yes"
		done=1
	}

	# covered only in the center
	overlap=which(temp2[,2]>a[i,2] & temp2[,3]<a[i,3] & done==0) 
	if(length(overlap)>0){
		overlap=overlap[which.max(temp2[overlap,3]-temp2[overlap,2])]
	}
	if(length(overlap)>0){
		a[i,5]="Yes"
	}
}

write.table(a[which(a[,5]=="Yes"),1:4],file=args[3],quote=F,sep="\t",row.names=F, col.names=F)
write.table(a[which(a[,5]=="No"),1:4],file=args[4],quote=F,sep="\t",row.names=F, col.names=F)
