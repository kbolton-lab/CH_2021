
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
id_tp53 = mut$id[which(mut$Gene.refGene == "TP53")]
id_tp53_LOH = cna$id[which(cna$chr == 17 & cna$start<=7.661779 & cna$end>=7.687550 & (cna$COPY_CHANGE == "loss" | cna$COPY_CHANGE == "neutral"))]
fac_tp53_wloh = as.numeric(is.element(id,intersect(id_tp53,id_tp53_LOH)))
fac_tp53_woloh = as.numeric(is.element(id,setdiff(id_tp53,id_tp53_LOH)))

fac_17p = as.numeric(apply(mat_cna[,c(129,131,133,135)],1,sum) != 0)
# fac_4q = as.numeric(apply(mat_cna[,c(26,28,30,32)],1,sum) != 0)

fac_tet2 = as.numeric(mat_mut[,"TET2"] != 0)
id_tet2 = mut$id[which(mut$Gene.refGene == "TET2")]
id_tet2_LOH = cna$id[which(cna$chr == 4 & cna$start<=105.145875 & cna$end>=105.279816 & (cna$COPY_CHANGE == "loss" | cna$COPY_CHANGE == "neutral"))]
fac_tet2_wloh = as.numeric(is.element(id,intersect(id_tet2,id_tet2_LOH)))
fac_tet2_woloh = as.numeric(is.element(id,setdiff(id_tet2,id_tet2_LOH)))

fac_jak2 = as.numeric(mat_mut[,"JAK2"] != 0)
id_jak2 = mut$id[which(mut$Gene.refGene == "JAK2")]
id_jak2_LOH = cna$id[which(cna$chr == 9 & cna$start<=4.984390 & cna$end>=5.128183 & cna$COPY_CHANGE == "neutral")]
fac_jak2_wloh = as.numeric(is.element(id,intersect(id_jak2,id_jak2_LOH)))
fac_jak2_woloh = as.numeric(is.element(id,setdiff(id_jak2,id_jak2_LOH)))

fac_dnmt3a = as.numeric(mat_mut[,"DNMT3A"] != 0)
id_dnmt3a = mut$id[which(mut$Gene.refGene == "DNMT3A")]
id_dnmt3a_LOH = cna$id[which(cna$chr == 2 & cna$start<=25.227855 & cna$end>=25.342590 & (cna$COPY_CHANGE == "loss" | cna$COPY_CHANGE == "neutral"))]
fac_dnmt3a_wloh = as.numeric(is.element(id,intersect(id_dnmt3a,id_dnmt3a_LOH)))
fac_dnmt3a_woloh = as.numeric(is.element(id,setdiff(id_dnmt3a,id_dnmt3a_LOH)))

is_17p = as.numeric(cna$chr == 17 & cna$start<=7.661779 & cna$end>=7.687550 & (cna$COPY_CHANGE == "loss" | cna$COPY_CHANGE == "neutral"))
is_4q = as.numeric(cna$chr == 4 & cna$start<=105.145875 & cna$end>=105.279816 & (cna$COPY_CHANGE == "loss" | cna$COPY_CHANGE == "neutral"))
is_9p = as.numeric(cna$chr == 9 & cna$start<=4.984390 & cna$end>=5.128183 & cna$COPY_CHANGE == "neutral")
is_2p = as.numeric(cna$chr == 2 & cna$start<=25.227855 & cna$end>=25.342590 & (cna$COPY_CHANGE == "loss" | cna$COPY_CHANGE == "neutral"))

cna_2 = cna[which(is_17p!=1 & is_4q!=1 & is_9p!=1 & is_2p!=1),]
#cna_2 = cna

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

fac_none = as.numeric(apply(mat_mut,1,function(x){sum(as.numeric(x))}) == 0)

n_cna = rep(0,length(id))
for(i in 1:nrow(cna)){
n_cna[which(id==cna$id[i])] = n_cna[which(id==cna$id[i])] + 1
}

n_cna_2 = rep(0,length(id))
for(i in 1:nrow(cna_2)){
n_cna_2[which(id==cna_2$id[i])] = n_cna_2[which(id==cna_2$id[i])] + 1
}

cna_all = numeric(max(n_cna)+1)
for(i in 0:max(n_cna)){cna_all[i+1] = sum(n_cna[which(fac_none==1)]==i)}
cs_cna_all = c(0,cumsum(cna_all/sum(cna_all)))

cna_tp53 = numeric(max(n_cna)+1)
for(i in 0:max(n_cna)){cna_tp53[i+1] = sum(n_cna[which(fac_tp53==1)]==i)}
cs_cna_tp53 = c(0,cumsum(cna_tp53/sum(cna_tp53)))

cna_tp53_wloh = numeric(max(n_cna_2)+1)
for(i in 0:max(n_cna_2)){cna_tp53_wloh[i+1] = sum(n_cna_2[which(fac_tp53_wloh==1)]==i)}
cs_cna_tp53_wloh = c(0,cumsum(cna_tp53_wloh/sum(cna_tp53_wloh)))

cna_tp53_woloh = numeric(max(n_cna_2)+1)
for(i in 0:max(n_cna_2)){cna_tp53_woloh[i+1] = sum(n_cna_2[which(fac_tp53_woloh==1)]==i)}
cs_cna_tp53_woloh = c(0,cumsum(cna_tp53_woloh/sum(cna_tp53_woloh)))

cna_dnmt3a = numeric(max(n_cna)+1)
for(i in 0:max(n_cna)){cna_dnmt3a[i+1] = sum(n_cna[which(fac_dnmt3a==1)]==i)}
cs_cna_dnmt3a = c(0,cumsum(cna_dnmt3a/sum(cna_dnmt3a)))

cna_dnmt3a_wloh = numeric(max(n_cna_2)+1)
for(i in 0:max(n_cna_2)){cna_dnmt3a_wloh[i+1] = sum(n_cna_2[which(fac_dnmt3a_wloh==1)]==i)}
cs_cna_dnmt3a_wloh = c(0,cumsum(cna_dnmt3a_wloh/sum(cna_dnmt3a_wloh)))

cna_dnmt3a_woloh = numeric(max(n_cna_2)+1)
for(i in 0:max(n_cna_2)){cna_dnmt3a_woloh[i+1] = sum(n_cna_2[which(fac_dnmt3a_woloh==1)]==i)}
cs_cna_dnmt3a_woloh = c(0,cumsum(cna_dnmt3a_woloh/sum(cna_dnmt3a_woloh)))

cna_tet2 = numeric(max(n_cna)+1)
for(i in 0:max(n_cna)){cna_tet2[i+1] = sum(n_cna[which(fac_tet2==1)]==i)}
cs_cna_tet2 = c(0,cumsum(cna_tet2/sum(cna_tet2)))

cna_tet2_wloh = numeric(max(n_cna_2)+1)
for(i in 0:max(n_cna_2)){cna_tet2_wloh[i+1] = sum(n_cna_2[which(fac_tet2_wloh==1)]==i)}
cs_cna_tet2_wloh = c(0,cumsum(cna_tet2_wloh/sum(cna_tet2_wloh)))

cna_tet2_woloh = numeric(max(n_cna_2)+1)
for(i in 0:max(n_cna_2)){cna_tet2_woloh[i+1] = sum(n_cna_2[which(fac_tet2_woloh==1)]==i)}
cs_cna_tet2_woloh = c(0,cumsum(cna_tet2_woloh/sum(cna_tet2_woloh)))

cna_jak2 = numeric(max(n_cna)+1)
for(i in 0:max(n_cna)){cna_jak2[i+1] = sum(n_cna[which(fac_jak2==1)]==i)}
cs_cna_jak2 = c(0,cumsum(cna_jak2/sum(cna_jak2)))

cna_jak2_wloh = numeric(max(n_cna_2)+1)
for(i in 0:max(n_cna_2)){cna_jak2_wloh[i+1] = sum(n_cna_2[which(fac_jak2_wloh==1)]==i)}
cs_cna_jak2_wloh = c(0,cumsum(cna_jak2_wloh/sum(cna_jak2_wloh)))

cna_jak2_woloh = numeric(max(n_cna_2)+1)
for(i in 0:max(n_cna_2)){cna_jak2_woloh[i+1] = sum(n_cna_2[which(fac_jak2_woloh==1)]==i)}
cs_cna_jak2_woloh = c(0,cumsum(cna_jak2_woloh/sum(cna_jak2_woloh)))


#n_cna = apply(mat_cna,1,sum)
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

p = rep(NA,8)
p[1] = wilcox.test(n_cna[which(fac_tp53==1)],n_cna[which(fac_none==1)])$ p.value
p[2] = wilcox.test(n_cna[which(fac_dnmt3a==1)],n_cna[which(fac_none==1)])$ p.value
p[3] = wilcox.test(n_cna[which(fac_tet2==1)],n_cna[which(fac_none==1)])$ p.value
p[4] = wilcox.test(n_cna[which(fac_jak2==1)],n_cna[which(fac_none==1)])$ p.value
p[5] = wilcox.test(n_cna[which(fac_asxl1==1)],n_cna[which(fac_none==1)])$ p.value
p[6] = wilcox.test(n_cna[which(fac_ppm1d==1)],n_cna[which(fac_none==1)])$ p.value
p[7] = wilcox.test(n_cna[which(fac_sf3b1==1)],n_cna[which(fac_none==1)])$ p.value
p[8] = wilcox.test(n_cna[which(fac_gnb1==1)],n_cna[which(fac_none==1)])$ p.value
p[9] = wilcox.test(n_cna[which(fac_cbl==1)],n_cna[which(fac_none==1)])$ p.value
p[10] = wilcox.test(n_cna[which(fac_srsf2==1)],n_cna[which(fac_none==1)])$ p.value
p[11] = wilcox.test(n_cna[which(fac_u2af1==1)],n_cna[which(fac_none==1)])$ p.value
p[12] = wilcox.test(n_cna[which(fac_gnas==1)],n_cna[which(fac_none==1)])$ p.value

n_cna[which(fac_tp53_woloh==1)]
n_cna[which(fac_tp53_wloh==1)]
n_cna[which(fac_dnmt3a_woloh==1)]
n_cna[which(fac_dnmt3a_wloh==1)]
n_cna[which(fac_tet2_woloh==1)]
n_cna[which(fac_tet2_wloh==1)]
n_cna[which(fac_jak2_woloh==1)]
n_cna[which(fac_jak2_wloh==1)]

p2 = rep(NA,5)
p2[1] = wilcox.test(n_cna[which(fac_tp53_woloh==1)],n_cna[which(fac_tp53_wloh==1)]-1)$ p.value
p2[2] = wilcox.test(n_cna[which(fac_dnmt3a_woloh==1)],n_cna[which(fac_dnmt3a_wloh==1)]-1)$ p.value
p2[3] = wilcox.test(n_cna[which(fac_tet2_woloh==1)],n_cna[which(fac_tet2_wloh==1)]-1)$ p.value
p2[4] = wilcox.test(n_cna[which(fac_jak2_woloh==1)],n_cna[which(fac_jak2_wloh==1)]-1)$ p.value
p2[5] = wilcox.test(
c(n_cna[which(fac_tp53_woloh==1)],n_cna[which(fac_dnmt3a_woloh==1)],n_cna[which(fac_tet2_woloh==1)],n_cna[which(fac_jak2_woloh==1)]),
c(n_cna[which(fac_tp53_wloh==1)],n_cna[which(fac_dnmt3a_wloh==1)],n_cna[which(fac_tet2_wloh==1)],n_cna[which(fac_jak2_wloh==1)])-1
)$p.value

q = p.adjust(p,"BH")
q2 = p.adjust(p2,"BH")


pdf("figure_gene_n_cna_compare_all.pdf",width=20,height=15)

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
#col_p = cm.colors((length(cs_cna_all)-1))
#col_p = colPal2((length(cs_cna_all)))
col_p = colPal2(5) # edit here
col_p = col_p[-1]
col_p = c(col_p,rep(col_p[length(col_p)],10))
#col_p = col_p[((length(col_p)/2)+1):length(col_p)]
for(i in 2:length(cs_cna_all)){
	rect(0.5-0.3,(1-cs_cna_all[i])*y_mag,0.5+0.3,(1-cs_cna_all[i+1])*y_mag,border="transparent",col=col_p[i-1])
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

for(i in 1:length(q)){
	if(q[i]<0.001){
		text(i+0.5,ts[i]*y_mag+0.3,"***",adj=0.5,cex=2)
	}else if(q[i]<0.01){
		text(i+0.5,ts[i]*y_mag+0.3,"**",adj=0.5,cex=2)
	}else if(q[i]<0.1){
		text(i+0.5,ts[i]*y_mag+0.3,"*",adj=0.5,cex=2)
	}
}

text(0.5,-0.5,"No mutation",srt=90,font=2,cex=1.5,adj=1)
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

## significance annotation
text(14,8.8,"q value",adj=0,cex=1.5,font=2)
text(14.5,8,"***",adj=0.5,cex=2); text(15,8,"< 0.001",adj=0,cex=1.5,font=2)
text(14.5,7.5,"**",adj=0.5,cex=2); text(15,7.5,"< 0.01",adj=0,cex=1.5,font=2)
text(14.5,7,"*",adj=0.5,cex=2); text(15,7,"< 0.1",adj=0,cex=1.5,font=2)

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



pdf("figure_tp53_cna_compare_loh.pdf",width=20,height=15)

par(mfrow=c(1,1))

plot(
NULL,NULL,
xlab="",ylab="",
xaxt="n",yaxt="n",
bty="n",
xlim=c(-3,10),
ylim=c(-3,6.5)
)

y_mag = 10

colPal2 <- colorRampPalette(c("white", "red"))
#col_p = cm.colors((length(cs_cna_all)-1))
#col_p = colPal2((length(cs_cna_all)))
col_p = colPal2(4)
col_p = col_p[-1]
col_p = c(col_p,rep(col_p[length(col_p)],10))
#col_p = col_p[((length(col_p)/2)+1):length(col_p)]
for(i in 2:length(cs_cna_all)){
	rect(1.5-0.3,(1-cs_cna_tp53_woloh[i])*y_mag,1.5,(1-cs_cna_tp53_woloh[i+1])*y_mag,border="transparent",col=col_p[i-1])
	rect(1.5,(1-cs_cna_tp53_wloh[i])*y_mag,1.5+0.3,(1-cs_cna_tp53_wloh[i+1])*y_mag,border="transparent",col=col_p[i-1])
	rect(2.5-0.3,(1-cs_cna_dnmt3a_woloh[i])*y_mag,2.5,(1-cs_cna_dnmt3a_woloh[i+1])*y_mag,border="transparent",col=col_p[i-1])
	rect(2.5,(1-cs_cna_dnmt3a_wloh[i])*y_mag,2.5+0.3,(1-cs_cna_dnmt3a_wloh[i+1])*y_mag,border="transparent",col=col_p[i-1])
	rect(3.5-0.3,(1-cs_cna_tet2_woloh[i])*y_mag,3.5,(1-cs_cna_tet2_woloh[i+1])*y_mag,border="transparent",col=col_p[i-1])
	rect(3.5,(1-cs_cna_tet2_wloh[i])*y_mag,3.5+0.3,(1-cs_cna_tet2_wloh[i+1])*y_mag,border="transparent",col=col_p[i-1])
	rect(4.5-0.3,(1-cs_cna_jak2_woloh[i])*y_mag,4.5,(1-cs_cna_jak2_woloh[i+1])*y_mag,border="transparent",col=col_p[i-1])
	rect(4.5,(1-cs_cna_jak2_wloh[i])*y_mag,4.5+0.3,(1-cs_cna_jak2_wloh[i+1])*y_mag,border="transparent",col=col_p[i-1])
}

tet2_max = max((1-cs_cna_tet2_woloh[-1])*y_mag,(1-cs_cna_tet2_wloh[-1])*y_mag)
jak2_max = max((1-cs_cna_jak2_woloh[-1])*y_mag,(1-cs_cna_jak2_wloh[-1])*y_mag)

#segments(3.5-0.3,tet2_max+y_mag/60,3.5+0.3,tet2_max+y_mag/60,lwd=2)
#segments(4.5-0.3,jak2_max+y_mag/60,4.5+0.3,jak2_max+y_mag/60,lwd=2)

#theta = seq(pi/2,pi/2+2*pi*4/5,length=5)
#r = y_mag/200
#x0 = 3.5; y0 = tet2_max+y_mag*2/60
#segments(x0,y0,x0+r*cos(theta),y0+r*sin(theta),lwd=6)
#x0 = 4.5; y0 = jak2_max+y_mag*2/60
#segments(x0,y0,x0+r*cos(theta),y0+r*sin(theta),lwd=6)



ts = c(
1-cs_cna_tp53[2],
1-cs_cna_dnmt3a[2],
1-cs_cna_tet2[2],
1-cs_cna_asxl1[2],
1-cs_cna_ppm1d[2],
1-cs_cna_jak2[2],
1-cs_cna_sf3b1[2],
1-cs_cna_gnb1[2],
1-cs_cna_cbl[2],
1-cs_cna_srsf2[2],
1-cs_cna_u2af1[2],
1-cs_cna_gnas[2]
)



text(-.3,-.5,"LOH",srt=0,font=2,cex=1.5,adj=0.5)
text(-.3,-1,"Gene",srt=0,font=2,cex=1.5,adj=0.5)

#text(0.5,-1,"No mutation",srt=0,font=2,cex=1.5,adj=0.5)

text(1.5,-1,"TP53",srt=0,font=4,cex=1.5,adj=0.5)
text(1.3,-.5,"-",srt=0,font=2,cex=1.5,adj=0.5)
text(1.7,-.5,"+",srt=0,font=2,cex=1.5,adj=0.5)
text(2.5,-1,"DNMT3A",srt=0,font=4,cex=1.5,adj=0.5)
text(2.3,-.5,"-",srt=0,font=2,cex=1.5,adj=0.5)
text(2.7,-.5,"+",srt=0,font=2,cex=1.5,adj=0.5)
text(3.5,-1,"TET2",srt=0,font=4,cex=1.5,adj=0.5)
text(3.3,-.5,"-",srt=0,font=2,cex=1.5,adj=0.5)
text(3.7,-.5,"+",srt=0,font=2,cex=1.5,adj=0.5)
text(4.5,-1,"JAK2",srt=0,font=4,cex=1.5,adj=0.5)
text(4.3,-.5,"-",srt=0,font=2,cex=1.5,adj=0.5)
text(4.7,-.5,"+",srt=0,font=2,cex=1.5,adj=0.5)


segments(0,0,0,0.6*y_mag)
lat1 = 0.1*(0:6)
lat2 = 0.1*(1:6) - 0.05
for(i in 1:length(lat1)){
	segments(0,lat1[i]*y_mag,0-0.1,lat1[i]*y_mag)
	text(0-0.3,lat1[i]*y_mag,lat1[i]*100,adj=0.5,cex=1.5,font=2,srt=90)
}
for(i in 1:length(lat2)){
	segments(0,lat2[i]*y_mag,0-0.05,lat2[i]*y_mag)
}
text(-1.3,0.25*y_mag,"Proportion of subjects (%)",srt=90,font=2,cex=1.8)

dev.off()


## End of analysis ##
q()
