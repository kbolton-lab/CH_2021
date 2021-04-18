
## import data
mut_path = "MatMed_2021_mut_table.txt"
cna_path = "MatMed_2021_cna_table.csv"  
id = as.character(1:11234)


# get info of chr length
arminfo = read.table(file="cytoBand/arm_chr_length.txt")
cent = arminfo[1:22,4]


mut = read.table(mut_path,header=T,stringsAsFactor=F,sep="\t",quote="",colClass="character")
mut = mut[which(!is.na(mut$Chr)),]
mat_mut = matrix(0,ncol=length(unique(mut$Gene.refGene)),nrow=length(id))
colnames(mat_mut) = unique(mut$Gene.refGene); rownames(mat_mut) = id



# get info of chr length
arminfo = read.table("cytoband/arm_chr_length.txt")
cent = arminfo[1:22,4]

mat_mut = matrix(0,ncol=length(unique(mut$Gene.refGene)),nrow=length(id))
colnames(mat_mut) = unique(mut$Gene.refGene); rownames(mat_mut) = id

mut_count = matrix(0,ncol=7,nrow=length(id))
rownames(mut_count) = id
n_mut = rep(0,length(id)); names(n_mut) = id

for(i in 1:nrow(mut)){
	cur_id = mut$id[i]
	cur_gene = mut$Gene.refGene[i]

	if(! is.element(cur_id,id)){
		next
	}
	
	n_mut[cur_id] = n_mut[cur_id] + 1
	mat_mut[cur_id,cur_gene] = max(mat_mut[cur_id,cur_gene],mut$misRate[i])
}



cna = read.csv(cna_path,header=T,stringsAsFactor=F,quote="",comment.char="")
cna = cna[which(!is.na(cna$chr)),]
label_cna = paste0(sort(rep(1:22,4*2)),":",sort(rep(c("CNN-LOH","Duplication","Deletion","Unknown"),2)),":",c("p","q"))
mat_cna =  matrix(0,ncol=length(label_cna),nrow=length(id))
colnames(mat_cna) = label_cna; rownames(mat_cna) = id
n_cna = rep(0,length(id)); names(n_cna) = rownames(mat_cna)

for(i in 1:nrow(cna)){
	cur_id = cna$id[i]
	
	cur_type = "unknown"
	if(cna$COPY_CHANGE[i] == "loss"){
		cur_type = "Deletion"
	}else if(cna$COPY_CHANGE[i] == "neutral"){
		cur_type = "CNN-LOH"
	}else if(cna$COPY_CHANGE[i] == "gain"){
		cur_type = "Duplication"
	}else{
		cur_type = "Unknown"
	}
	cur_cna = paste0(cna$chr[i],":",cur_type)
	if(! is.element(cur_id,id)){
		next
	}
	n_cna[cur_id] = n_cna[cur_id] + 1
	
	cur_cent = cent[as.numeric(cna$chr[i])]
	if(cna$start[i]*(10^6) <= cur_cent & cna$end[i]*(10^6) <= cur_cent){
		cur_cna = paste0(cur_cna,":p")
	}else if(cna$end[i]*(10^6) >= cur_cent & cna$start[i]*(10^6) >= cur_cent){
		cur_cna = paste0(cur_cna,":q")
	}else if(cna$start[i]*(10^6) <= cur_cent & cna$end[i]*(10^6) >= cur_cent)
	if(! is.element(cur_cna,label_cna)){
		next
	}
	
	mat_cna[cur_id,cur_cna] = 1
	# mat_cna[cur_id,cur_cna] = mat_cna[cur_id,cur_cna] + 1
}


fac_tp53 = as.numeric(mat_mut[,"TP53"] != 0)
fac_tet2 = as.numeric(mat_mut[,"TET2"] != 0)
fac_jak2 = as.numeric(mat_mut[,"JAK2"] != 0)
fac_dnmt3a = as.numeric(mat_mut[,"DNMT3A"] != 0)
fac_gnb1 = as.numeric(mat_mut[,"GNB1"] != 0)
fac_gnas = as.numeric(mat_mut[,"GNAS"] != 0)
fac_sf3b1 = as.numeric(mat_mut[,"SF3B1"] != 0)
fac_srsf2 = as.numeric(mat_mut[,"SRSF2"] != 0)
fac_asxl1 = as.numeric(mat_mut[,"ASXL1"] != 0)
fac_ppm1d = as.numeric(mat_mut[,"PPM1D"] != 0)
fac_cbl = as.numeric(mat_mut[,"CBL"] != 0)
fac_srsf2 = as.numeric(mat_mut[,"SRSF2"] != 0)
fac_u2af1 = as.numeric(mat_mut[,"U2AF1;U2AF1L5"] != 0)
fac_gnas = as.numeric(mat_mut[,"GNAS"] != 0)

n_cna = rep(0,length(id))
for(i in 1:nrow(cna)){
n_cna[which(id==cna$id[i])] = n_cna[which(id==cna$id[i])] + 1
}

n_cna = n_cna + n_mut
n_cna[which(n_mut!=0)] = n_cna[which(n_mut!=0)] - 1

cna_tp53 = numeric(max(n_cna)+1)
for(i in 0:max(n_cna)){cna_tp53[i+1] = sum(n_cna[which(fac_tp53==1)]==i)}
cs_cna_tp53 = c(0,cumsum(cna_tp53/sum(cna_tp53)))

cna_dnmt3a = numeric(max(n_cna)+1)
for(i in 0:max(n_cna)){cna_dnmt3a[i+1] = sum(n_cna[which(fac_dnmt3a==1)]==i)}
cs_cna_dnmt3a = c(0,cumsum(cna_dnmt3a/sum(cna_dnmt3a)))

cna_tet2 = numeric(max(n_cna)+1)
for(i in 0:max(n_cna)){cna_tet2[i+1] = sum(n_cna[which(fac_tet2==1)]==i)}
cs_cna_tet2 = c(0,cumsum(cna_tet2/sum(cna_tet2)))

cna_jak2 = numeric(max(n_cna)+1)
for(i in 0:max(n_cna)){cna_jak2[i+1] = sum(n_cna[which(fac_jak2==1)]==i)}
cs_cna_jak2 = c(0,cumsum(cna_jak2/sum(cna_jak2)))

cna_asxl1 = numeric(max(n_cna)+1)
for(i in 0:max(n_cna)){cna_asxl1[i+1] = sum(n_cna[which(fac_asxl1==1)]==i)}
cs_cna_asxl1 = c(0,cumsum(cna_asxl1/sum(cna_asxl1)))

cna_ppm1d = numeric(max(n_cna)+1)
for(i in 0:max(n_cna)){cna_ppm1d[i+1] = sum(n_cna[which(fac_ppm1d==1)]==i)}
cs_cna_ppm1d = c(0,cumsum(cna_ppm1d/sum(cna_ppm1d)))

cna_sf3b1 = numeric(max(n_cna)+1)
for(i in 0:max(n_cna)){cna_sf3b1[i+1] = sum(n_cna[which(fac_sf3b1==1)]==i)}
cs_cna_sf3b1 = c(0,cumsum(cna_sf3b1/sum(cna_sf3b1)))

cna_gnb1 = numeric(max(n_cna)+1)
for(i in 0:max(n_cna)){cna_gnb1[i+1] = sum(n_cna[which(fac_gnb1==1)]==i)}
cs_cna_gnb1 = c(0,cumsum(cna_gnb1/sum(cna_gnb1)))

cna_cbl = numeric(max(n_cna)+1)
for(i in 0:max(n_cna)){cna_cbl[i+1] = sum(n_cna[which(fac_cbl==1)]==i)}
cs_cna_cbl = c(0,cumsum(cna_cbl/sum(cna_cbl)))

cna_srsf2 = numeric(max(n_cna)+1)
for(i in 0:max(n_cna)){cna_srsf2[i+1] = sum(n_cna[which(fac_srsf2==1)]==i)}
cs_cna_srsf2 = c(0,cumsum(cna_srsf2/sum(cna_srsf2)))

cna_u2af1 = numeric(max(n_cna)+1)
for(i in 0:max(n_cna)){cna_u2af1[i+1] = sum(n_cna[which(fac_u2af1==1)]==i)}
cs_cna_u2af1 = c(0,cumsum(cna_u2af1/sum(cna_u2af1)))

cna_gnas = numeric(max(n_cna)+1)
for(i in 0:max(n_cna)){cna_gnas[i+1] = sum(n_cna[which(fac_gnas==1)]==i)}
cs_cna_gnas = c(0,cumsum(cna_gnas/sum(cna_gnas)))




pdf("figure_gene_n_alt_compare.pdf",width=20,height=15)

par(mfrow=c(1,1))

plot(
NULL,NULL,
xlab="",ylab="",
xaxt="n",yaxt="n",
bty="n",
xlim=c(-3,20),
ylim=c(-3,10)
)

y_mag = 10

colPal2 <- colorRampPalette(c("white", "red"))
col_p = colPal2(5)   ### edit here
col_p = col_p[-1]
col_p = c(col_p,rep(col_p[length(col_p)],10))
for(i in 2:length(cs_cna_tp53)){
	rect(1.5-0.3,(1-cs_cna_tp53[i])*y_mag,1.5+0.3,(1-cs_cna_tp53[i+1])*y_mag,border="transparent",col=col_p[i-1])
	rect(2.5-0.3,(1-cs_cna_dnmt3a[i])*y_mag,2.5+0.3,(1-cs_cna_dnmt3a[i+1])*y_mag,border="transparent",col=col_p[i-1])
	rect(3.5-0.3,(1-cs_cna_tet2[i])*y_mag,3.5+0.3,(1-cs_cna_tet2[i+1])*y_mag,border="transparent",col=col_p[i-1])
	rect(4.5-0.3,(1-cs_cna_jak2[i])*y_mag,4.5+0.3,(1-cs_cna_jak2[i+1])*y_mag,border="transparent",col=col_p[i-1])
	rect(5.5-0.3,(1-cs_cna_asxl1[i])*y_mag,5.5+0.3,(1-cs_cna_asxl1[i+1])*y_mag,border="transparent",col=col_p[i-1])
	rect(6.5-0.3,(1-cs_cna_ppm1d[i])*y_mag,6.5+0.3,(1-cs_cna_ppm1d[i+1])*y_mag,border="transparent",col=col_p[i-1])
	rect(7.5-0.3,(1-cs_cna_sf3b1[i])*y_mag,7.5+0.3,(1-cs_cna_sf3b1[i+1])*y_mag,border="transparent",col=col_p[i-1])
	rect(8.5-0.3,(1-cs_cna_gnb1[i])*y_mag,8.5+0.3,(1-cs_cna_gnb1[i+1])*y_mag,border="transparent",col=col_p[i-1])
	rect(9.5-0.3,(1-cs_cna_cbl[i])*y_mag,9.5+0.3,(1-cs_cna_cbl[i+1])*y_mag,border="transparent",col=col_p[i-1])
	rect(10.5-0.3,(1-cs_cna_srsf2[i])*y_mag,10.5+0.3,(1-cs_cna_srsf2[i+1])*y_mag,border="transparent",col=col_p[i-1])
	rect(11.5-0.3,(1-cs_cna_u2af1[i])*y_mag,11.5+0.3,(1-cs_cna_u2af1[i+1])*y_mag,border="transparent",col=col_p[i-1])
	rect(12.5-0.3,(1-cs_cna_gnas[i])*y_mag,12.5+0.3,(1-cs_cna_gnas[i+1])*y_mag,border="transparent",col=col_p[i-1])
}

ts = c(
1-cs_cna_tp53[2],
1-cs_cna_dnmt3a[2],
1-cs_cna_tet2[2],
1-cs_cna_jak2[2],
1-cs_cna_asxl1[2],
1-cs_cna_ppm1d[2],
1-cs_cna_sf3b1[2],
1-cs_cna_gnb1[2],
1-cs_cna_cbl[2],
1-cs_cna_srsf2[2],
1-cs_cna_u2af1[2],
1-cs_cna_gnas[2]
)

text(1.5,-0.5,"TP53",srt=90,font=4,cex=1.5,adj=1)
text(2.5,-0.5,"DNMT3A",srt=90,font=4,cex=1.5,adj=1)
text(3.5,-0.5,"TET2",srt=90,font=4,cex=1.5,adj=1)
text(4.5,-0.5,"JAK2",srt=90,font=4,cex=1.5,adj=1)
text(5.5,-0.5,"ASXL1",srt=90,font=4,cex=1.5,adj=1)
text(6.5,-0.5,"PPM1D",srt=90,font=4,cex=1.5,adj=1)
text(7.5,-0.5,"SF3B1",srt=90,font=4,cex=1.5,adj=1)
text(8.5,-0.5,"GNB1",srt=90,font=4,cex=1.5,adj=1)
text(9.5,-0.5,"CBL",srt=90,font=4,cex=1.5,adj=1)
text(10.5,-0.5,"SRSF2",srt=90,font=4,cex=1.5,adj=1)
text(11.5,-0.5,"U2AF1",srt=90,font=4,cex=1.5,adj=1)
text(12.5,-0.5,"GNAS",srt=90,font=4,cex=1.5,adj=1)

segments(0,0,0,1*y_mag)
lat1 = 0.1*(0:10)
lat2 = 0.1*(1:10) - 0.05
for(i in 1:length(lat1)){
	segments(0,lat1[i]*y_mag,0-0.1,lat1[i]*y_mag)
	text(0-0.3,lat1[i]*y_mag,lat1[i]*100,adj=0.5,cex=1.5,font=2,srt=90)
}
for(i in 1:length(lat2)){
	segments(0,lat2[i]*y_mag,0-0.05,lat2[i]*y_mag)
}
text(-1.3,0.25*y_mag,"Proportion of subjects (%)",srt=90,font=2,cex=1.5)

## color annotation

rect(14+0.5*(1-1),2,14.5+0.5*(1-1),2.5,col=col_p[1],border="transparent")
text(14.25+0.5*(1-1),2.8,1,adj=0.5,cex=1.2,font=2)
rect(14+0.5*(2-1),2,14.5+0.5*(2-1),2.5,col=col_p[2],border="transparent")
text(14.25+0.5*(2-1),2.8,2,adj=0.5,cex=1.2,font=2)
rect(14+0.5*(3-1),2,14.5+0.5*(3-1),2.5,col=col_p[3],border="transparent")
text(14.25+0.5*(3-1),2.8,">=3",adj=0.5,cex=1.2,font=2)

text(14,4.3,"Number of CNA",adj=0,cex=1.3,font=2)
#text(14,3.8,"(Other than 2pLOH, 4qLOH, 9pUPD, 17pLOH)",adj=0,cex=1.3,font=2)

dev.off()


## End of analysis ##
q()



