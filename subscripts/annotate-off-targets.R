
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

a=read.table(file=args[1],header=F)
b=read.table(file=args[2],header=F)

com=which(is.element(a[,4],b[,1]))

for (i in com){
	sel=which(b[,1]==a[i,4])
	if(length(sel)>0){
		for(j in 1:length(sel)){
			a[i,4]=paste(a[i,4],b[sel[j],2],sep=",")
		}
	}
}

write.table(a,file=args[3],quote=F,sep="\t",row.names=F, col.names=F)

