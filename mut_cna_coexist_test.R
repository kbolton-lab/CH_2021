
## draw Supplementary Fig.6a ##

# library
library(tidyverse)

genes = c("DNMT3A","TET2","JAK2","TP53","NRAS","EZH2","CBL")
chr = c(2,4,9,17,1,7,11)
starts = c(25455830,106034562,4949510,7563450,115247085,148504464,119076986)
ends = c(25565459,106234240,5163918,7592173,115259515,148581441,119178859)
type = c("LOH","LOH","UPD","LOH","UPD","LOH","LOH")
gene_reg = data.frame(genes,chr,starts,ends,type)

# CH
mut_path = "MatMed_2021_mut_table.txt"
cna_path = "MatMed_2021_cna_table.csv"

mut = read.table(mut_path,header=T,stringsAsFactor=F,sep="\t",quote="",colClass="character")
mut = mut[which(!is.na(mut$Chr)),]
mut$Gene.refGene[which(mut$Gene.refGene=="TET2;TET2")] = "TET2"
mut$Gene.refGene[which(mut$Gene.refGene=="U2AF1;U2AF1L5 ")] = "U2AF1"

cna = read.csv(cna_path,header=T,stringsAsFactor=F,quote="",comment.char="")
cna = cna[which(!is.na(cna$chr)),]

id = as.character(1:11234)
starts = starts/(10^6)
ends = ends/(10^6)
type = c("LOH","LOH","neutral","LOH","neutral","LOH","LOH")

# g = as.numeric(commandArgs(trailingOnly=T)[1])

pdf(paste0("cf_plot.pdf"),width=10,height=10) 
par(mfrow=c(2,2),mar=c(.5,.5,.5,.5),omi=c(0,0,0,0))

n_coex = rep(NA,4)
xisq_list = xisq_null_list = list()
for(g in 1:4){
res = matrix(ncol=2,nrow=length(id))
rownames(res) = id

for(i in 1:length(id)){
	cur_id = id[i]
	cur_mut = mut[which(mut$id == cur_id & mut$Gene.refGene == genes[g]),]
	cur_cna = cna[which(cna$id == cur_id & cna$COPY_CHANGE == "neutral" & cna$chr==chr[g] & (as.numeric(cna$start)-20-ends[g])*(as.numeric(cna$end)+20-starts[g])<=0 ),]
	
	if(nrow(cur_mut) != 0){
		res[i,1] = max(as.numeric(cur_mut$misRate))
	}else{
		res[i,1] = 0
	}
	if(nrow(cur_cna) != 0){
		res[i,2] = max(as.numeric(cur_cna$CELL_FRAC))
	}else{
		res[i,2] = 0
	}
}

res_coex = res[which(res[,1] != 0 & res[,2] != 0),]
n_coex[g] = nrow(res_coex)

d_mut = density(res[which(res[,1]!=0),1],bw=.1,from=0,to=0.5)
d_cna = density(res[which(res[,2]!=0),2],bw=.1,from=0,to=0.5)
dist_mat_mut = matrix(nrow=length(d_mut$y),ncol=length(d_mut$y)) 
for(i in 1:ncol(dist_mat_mut)) dist_mat_mut[,i] = d_mut$y
dist_mat_cna = matrix(nrow=length(d_cna$y),ncol=length(d_cna$y)) 
for(i in 1:nrow(dist_mat_cna)) dist_mat_cna[i,] = d_cna$y
dist_mat = dist_mat_mut * dist_mat_cna
dist_mat = dist_mat / sum(dist_mat)

obsv_mat = matrix(0,ncol=ncol(dist_mat),nrow=nrow(dist_mat))
hig = d_mut$x; low = c(0,hig[-length(hig)])

xisq = 0
for(i in 1:nrow(res_coex)){
	if(res_coex[i,1] >= res_coex[i,2]){ xisq = xisq + 1 }
}

print(xisq)

N_sim = 100000
xisq_null = rep(NA,N_sim)
for(k in 1:N_sim){
	xi = sample(1:512,nrow(res_coex),prob=d_mut$y)
	yi = sample(1:512,nrow(res_coex),prob=d_cna$y)
	xisq_null[k] = sum(xi >= yi)
}

p = sum(xisq_null>=xisq)/length(xisq_null)


plot(NULL,NULL,axes=F,xlab="CF of CNA",ylab="VAF of SNV/indel",xlim=c(-.1,1.6),ylim=c(-.1,1.6))
r = 0.008
theta = seq(-pi,pi,length=100)
for(i in 1:nrow(res_coex)){
	xs = res_coex[i,2] + r * cos(theta)
	ys = res_coex[i,1] + r * sin(theta)
	if(res_coex[i,2]<res_coex[i,1]){
		polygon(xs,ys,col=rgb(0.9,0.05,0.05,alpha=0.8),border="transparent")
	}else{
		polygon(xs,ys,col=rgb(0.05,0.05,0.9,alpha=0.8),border="transparent")
	}
}
#par(new=T); contour(d_cna$x,d_mut$x,t(dist_mat),col="gray",axes=F,xlab="",ylab="",xlim=c(-.1,.6),ylim=c(-.1,.6))
segments(0,0,.5,0,col="black"); segments(0,0,0,.5,col="black")
segments(0,.5,.5,.5,col="black"); segments(.5,0,.5,.5,col="black")
segments(0,0,.5,.5,col="black",lty="dashed")
text(.25,.55,genes[g],cex=1.3,font=3)
text(.3,.45,paste0("P = ",round(p,5)))
text(-.07,.25,"VAF of SNV/indel",srt=90,adj=0.5)
text(.25,-.07,"Cell fraction of UPD",srt=0,adj=0.5)

lat1 = c(0,.1,.2,.3,.4,.5)
lat2 = c(.05,.15,.25,.35,.45)
for(i in 1:length(lat1)){
	segments(lat1[i],0,lat1[i],-.01)
	segments(0,lat1[i],-.01,lat1[i])
#	if(i >= 2){
		text(lat1[i],-.025,lat1[i],cex=.6)
		text(-.025,lat1[i],lat1[i],srt=90,cex=.6)
#	}
}
for(i in 1:length(lat2)){
	segments(lat2[i],0,lat2[i],-.01)
	segments(0,lat2[i],-.01,lat2[i])
}
#dev.off()

xisq_list[[g]] = xisq
xisq_null_list[[g]] = xisq_null
}

dev.off()


