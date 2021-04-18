
mut_path = ## import mutation file on JGA
cna_path = ## import cna file on JGA

id = as.character(1:11234)

mut = read.table(mut_path,header=T,stringsAsFactor=F,quote="",sep="\t",comment.char="")
mut = mut[which(is.element(mut$id,id)),]
mut_vaf = mut$misRate

cna = read.table(cna_path,header=T,stringsAsFactor=F,quote="",sep="\t",comment.char="")
cna_cf = cna$CELL_FRAC

median(na.omit(cna_cf))
max(na.omit(cna_cf))
min(na.omit(cna_cf))

baf = unlist(strsplit(cna$BAF_SE," "))
baf = as.numeric(baf[(1:(length(baf)/2))*2-1])
cna$mudiff = 2*0.01*baf

cna$AFgain = 2*cna$mudiff / (1-cna$mudiff)
cna$AFupd = cna$mudiff
cna$AFloss = 2*cna$mudiff / (1+cna$mudiff)

cna$cf =(cna$COPY_CHANGE=="gain")*cna$AFgain+
(cna$COPY_CHANGE=="loss")*cna$AFloss+
(cna$COPY_CHANGE=="neutral")*cna$AFupd+
(cna$COPY_CHANGE=="unknown")*cna$AFgain

cna_cf = cna$cf

print("mut vaf")
print("n"); length(mut_vaf)
print("median"); median(mut_vaf)
print("range"); min(mut_vaf); max(mut_vaf)

print("cna cf")
print("n"); length(cna_cf)
print("median"); median(na.omit(cna_cf))
print("range"); min(na.omit(cna_cf)); max(na.omit(cna_cf))

h1 = hist(log10(mut_vaf),xlim=c(-3,0),breaks=seq(-3,0,0.05),col="skyblue",border="transparent",ylim=c(0,250),main="",xaxt="n", yaxt="n")
dev.off()
h2 = hist(log10(cna_cf),xlim=c(-3,0),breaks=seq(-3,0,0.05),col="skyblue",border="transparent",ylim=c(0,250),main="",xaxt="n", yaxt="n")
dev.off()



pdf("figure_mut_cna_cf_hist_2.pdf")

plot(
NULL,NULL,
xlab="",ylab="",
xaxt="n",yaxt="n",
bty="n",
xlim=c(-4,1),
ylim=c(-100,800)
)

lat1 = c(0.001,0.01,0.1,1)
lat2 = c(0.001*(2:9),0.01*(2:9),0.1*(2:9))

for(i in 1:length(lat1)){
	segments(log10(lat1),0,log10(lat1),-10)
	text(log10(lat1),-30,lat1,cex=0.7)
	segments(log10(lat1)+log10(2),0+400,log10(lat1)+log10(2),-10+400)
	text(log10(lat1)+log10(2),-30+400,lat1,cex=0.7)
}

for(i in 1:length(lat2)){
	segments(log10(lat2),0,log10(lat2),-5)
	segments(log10(lat2)+log10(2),0+400,log10(lat2)+log10(2),-5+400)
}

lat1 = c(50,100,150,200,250)
lat2 = setdiff(c(1:25)*10,lat1)

for(i in 1:length(lat1)){
	segments(-3,lat1[i],-3-0.05,lat1[i])
	segments(-3+log10(2),lat1[i]+400,-3-0.05+log10(2),lat1[i]+400)
	text(-3-0.1,lat1[i],lat1[i],cex=0.7,adj=1)
	text(-3-0.1+log10(2),lat1[i]+400,lat1[i],cex=0.7,adj=1)
}
for(i in 1:length(lat2)){
	segments(-3,lat2[i],-3-0.025,lat2[i])
	segments(-3+log10(2),lat2[i]+400,-3-0.025+log10(2),lat2[i]+400)
}

segments(-3,0,0,0); text(log10(0.04),-70,"Cell fraction of CNA",cex=0.8)
segments(-3,0,-3,250); text(-3-0.5,125,"Count",cex=0.8,srt=90)
segments(-3+log10(2),400,0+log10(2),400); text(log10(0.04)+log10(2),400-70,"VAF of mutation",cex=0.8)
segments(-3+log10(2),400,-3+log10(2),650); text(-3-0.5+log10(2),125+400,"Count",cex=0.8,srt=90)

for(i in 1:length(h1$counts)){
	rect(h1$breaks[i]+log10(2),400,h1$breaks[i+1]+log10(2),h1$count[i]+400,col="skyblue",border="transparent")
	#rect(h1$breaks[i]+log10(2),300,h1$breaks[i+1]+log10(2),h1$count[i]+300,col="skyblue",border="transparent")
	rect(h2$breaks[i],0,h2$breaks[i+1],h2$count[i],col="skyblue",border="transparent")
}
segments(-2,0,-2,650,col="red") 

dev.off()


