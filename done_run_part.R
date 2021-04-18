library(tidyverse)


input = commandArgs(trailingOnly=T)[1]
color = commandArgs(trailingOnly=T)[2]


## Begining of PART ##
## This part of source code was provided by Dr. Atsushi Niida ##

library("poibin")

PART<-function(M){
Mk<-apply(M,1,sum)
Mpp<-apply(M,2,mean);
#Mpv<-1-ppoibin(kk=Mk, pp=Mpp, method = "DFT-CF")
Mpv<-1-ppoibin(kk=Mk, pp=Mpp, method = "RNA")
Mpv<-p.adjust(Mpv,method="fdr")
Mpv<- -log10(Mpv)
names(Mk)<-rownames(M)
names(Mpv)<-rownames(M)
									 
#return(Mpv)
										  
###extrapolate p-values because ppoibin cannot calculate very small pvalues###
data<-cbind(Mk, Mpv)
colnames(data)<-c("k","p")

data<-as.data.frame(data[which(data[,2]<14&data[,2]>0.5),]) #pvalue 0.5~13?
result<-lm(p~1+k+I(k^2),data=data)
a<-sort(Mk)
x<-result$coefficients[1]+a*result$coefficients[2]+(a^2)*result$coefficients[3]

if(any(is.na(x))){
	x<-Mpv
	x[x==Inf]<-max(x[x!=Inf])
}else{
	tmp<-names(which(x<13))
	x[tmp]<-Mpv[tmp]
}

###comfirm extrapolation###
#plot(cbind(Mk, Mpv),xlim=c(0,max(Mk)),ylim=c(0,max(x)),xlab="k",ylab="p")
#par(new=T)
#plot(cbind(a,x),xlim=c(0,max(Mk)),ylim=c(0,max(x)),col="red",type="l",xlab="",ylab="")
#legend("topleft", legend = c("original","estimate"), col = c("black","red"),lty=1)  

return (rev(x))
}

## End of PART ##
## This part of source code was provided by Dr. Atsushi Niida ##


## marker density ##
N = 200

probe_file = "array_marker.txt"
probe_dat = read.table(probe_file,header=T,stringsAsFactor=F,sep="\t",quote="")
probe_dat = probe_dat[(1:round(nrow(probe_dat)/N))*N,]
sep = cbind(rep(0,23),1:23,rep(0,23)); colnames(sep) = colnames(probe_dat)
probe_dat = probe_dat %>% arrange(Chrom,BasePair)
probe_labels = paste0(probe_dat[,2],":",probe_dat[,3])
zero_labels = probe_labels[which(probe_dat[,3] == 0)]

cytoinfo <- read.table(file="cytoBand/cytoBand.txt")
arminfo <- read.table(file="cytoBand/arm_chr_length.txt")
cytoinfo$col <- ifelse(cytoinfo$V5=="gneg","white",ifelse(cytoinfo$V5=="gpos25","gray25",ifelse(cytoinfo$V5=="gpos50","gray",
ifelse(cytoinfo$V5=="gpos75","darkgray",ifelse(cytoinfo$V5=="gpos100","black",ifelse(cytoinfo$V5=="acen","red","blue4")))
)))
cytoinfo$col = as.character(cytoinfo$col)
cumpos = c(0,cumsum(arminfo$V2/(10^6)))
cumpos = cumpos[1:(length(cumpos)-1)]


cna_dat = read.table(input,header=F,stringsAsFactor=F,sep="\t",quote="")
colnames(cna_dat) = c("ID","chrom","loc.start","loc.end","num.mark","seg.mean")
cases = as.character(unique(cna_dat$ID))
cna_dat = cna_dat %>% filter(seg.mean != 0)

mat = matrix(0,nrow=length(probe_labels),ncol=length(cases))
colnames(mat) = cases; rownames(mat) = probe_labels

for(i in 1:nrow(mat)){
	chr = probe_dat[i,"Chrom"]
	bp = probe_dat[i,"BasePair"]
	cur_ID = unlist(cna_dat %>% filter(chrom == chr & loc.start <= bp & loc.end >= bp) %>% select(ID))
	mat[i,as.character(cur_ID)] = 1
}
mat[zero_labels,] = 0

res = PART(mat)
pos_sig = res[probe_labels]
pos_sig[which(pos_sig > 20)] = 20

sar = matrix(NA,ncol=5,nrow=length(which(pos_sig > -log10(0.25))))
sar[,1] = which(pos_sig > -log10(0.25))
tmp = unlist(strsplit(names(which(pos_sig > -log10(0.25))),":"))
sar[,2] = as.numeric(tmp[(1:(length(tmp)/2))*2-1])
sar[,3] = as.numeric(tmp[(1:(length(tmp)/2))*2])
sar[,4] = c(0,sar[2:nrow(sar),1]-sar[1:(nrow(sar)-1),1]-1)
sar[,5] = c(0,sar[2:nrow(sar),2]-sar[1:(nrow(sar)-1),2])

idx = sort(c(1,which(sar[,4]!=0)-1,which(sar[,4]!=0),nrow(sar)),decreasing=F)
idx1 = idx[(1:(length(idx)/2))*2-1]
idx2 = idx[(1:(length(idx)/2))*2]
res = cbind(sar[idx1,2],sar[idx1,3],sar[idx2,3])

output_file = paste0(input,"_output.txt")
write.table(res,file=output_file,col.name=F,row.name=F,sep="\t",quote=F)

output_pdf = paste0(input,"_output.pdf")
pdf(output_pdf,width=10,height=5)

plot(NULL,NULL,
xlim=c(-200,max(cumpos)+200),ylim=c(-max(pos_sig)/5,max(pos_sig)*1.1),
axes=FALSE,xlab="",ylab="")

rect(cumpos[1:(length(cumpos)-1)],-1,cumpos[2:length(cumpos)],-2,col=rep(c("white","gray"),20))
p1 = cumpos[as.numeric(probe_dat[1:(nrow(probe_dat)-1),"Chrom"])] + probe_dat[1:(nrow(probe_dat)-1),"BasePair"]/(10^6)
p2 = cumpos[as.numeric(probe_dat[2:nrow(probe_dat),"Chrom"])] + probe_dat[2:nrow(probe_dat),"BasePair"]/(10^6)
segments(p1,pos_sig[1:(length(pos_sig)-1)],p2,pos_sig[2:length(pos_sig)],col=color)
segments(0,-log10(0.25),cumpos[length(cumpos)],-log10(0.25),lwd=0.5,lty="dashed")
text((cumpos[1:(length(cumpos)-1)]+cumpos[2:length(cumpos)])/2,rep(c(-3.5,-4),20),1:22,cex=0.5,srt=90)

dev.off()


## End of analysis ##
q()
