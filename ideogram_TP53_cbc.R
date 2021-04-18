
mut_path = "MatMed_2021_mut_table.txt"
cna_path = "MatMed_2021_cna_table.csv"  
id = as.character(1:11234)

# get info of chr length
arminfo = read.table(file="cytoBand/arm_chr_length.txt")
cent = arminfo[1:22,4]

mut = read.table(mut_path,header=T,stringsAsFactor=F,sep="\t",quote="",colClass="character")
mut = mut[which(!is.na(mut$Chr)),]
mut$Gene.refGene[which(mut$Gene.refGene=="TET2;TET2")] = "TET2"
mut$Gene.refGene[which(mut$Gene.refGene=="U2AF1;U2AF1L5 ")] = "U2AF1"


cna = read.csv(cna_path,header=T,stringsAsFactor=F,quote="",comment.char="")
cna = cna[which(!is.na(cna$chr)),]


## 
library(tidyverse)

cytoinfo <- read.table("cytoBand/cytoBand.txt")
arminfo <- read.table("cytoBand/arm_chr_length.txt")
cytoinfo$col <- ifelse(cytoinfo$V5=="gneg","white",ifelse(cytoinfo$V5=="gpos25","gray25",ifelse(cytoinfo$V5=="gpos50","gray",
ifelse(cytoinfo$V5=="gpos75","darkgray",ifelse(cytoinfo$V5=="gpos100","black",ifelse(cytoinfo$V5=="acen","red","blue4")))
)))
cytoinfo$col = as.character(cytoinfo$col)

cum_chr_pos = c(0,cumsum(arminfo$V2/1000000/10))

cna = cna %>% 
mutate(TYPE = if_else(COPY_CHANGE == "loss","Deletion",if_else(COPY_CHANGE == "neutral","CNN-LOH",if_else(COPY_CHANGE == "gain","Duplication","unknown")))) %>%
mutate(BEG_GRCh37 = as.numeric(start) * (10^6)) %>% mutate(END_GRCh37 = as.numeric(end) * (10^6)) %>% mutate(CHROM = as.numeric(chr))
cna$COLOR = ifelse(cna$TYPE=="Duplication","firebrick3",ifelse(cna$TYPE=="Deletion","dodgerblue3",ifelse(cna$TYPE=="CNN-LOH","chartreuse3","tan")))
cna$DENSE = as.numeric(cna$CELL_FRAC)

pos = rep(Inf,length(id))
typ = rep("ZZ",length(id))
names(pos) = names(typ) = id

for(i in 1:nrow(cna)){
cur_id = cna$id[i]
cur_pos = as.numeric(cna$start[i]) + cum_chr_pos[as.numeric(cna$chr[i])]
if(!is.na(pos[cur_id]) & as.numeric(pos[cur_id])>cur_pos){
pos[cur_id] = cur_pos; typ[cur_id] = cna$COPY_CHANGE[i]
}
}

TP53_id = mut$id[which(mut$Gene.refGene == "TP53")]
LOH17p_id = cna$id[which(cna$chr==17 & (as.numeric(cna$end)-7571720/1000000)*(as.numeric(cna$start)-7590868/1000000)<=0)]



#### positive in ddPCR
tp53_ddpcr_id = c(
"2159",
"1027",
"6204",
"19",
"6201",
"10885",
"11228",
"10932",
"2708",
"4304",
"2118",
"10922",
"7914",
"10235",
"9758",
"9208",
"10790",
"6017",
"11075",
"3474",
"5543",
"32"
)



TP53_id = union(TP53_id,tp53_ddpcr_id)

TP53_or_LOH17p_id = union(union(TP53_id,LOH17p_id),tp53_ddpcr_id)
TP53_and_LOH17p_id = tp53_ddpcr_id #intersect(TP53_id,LOH17p_id)

pos = rep(Inf,length(TP53_or_LOH17p_id))
typ = rep("ZZ",length(TP53_or_LOH17p_id))
names(pos) = names(typ) = TP53_or_LOH17p_id

for(i in 1:nrow(cna)){
cur_id = cna$id[i]
if(!is.element(cur_id,TP53_or_LOH17p_id)){ next }
cur_pos = as.numeric(cna$start[i]) + cum_chr_pos[as.numeric(cna$chr[i])]
if(as.numeric(pos[as.character(cur_id)])>cur_pos){
pos[as.character(cur_id)] = cur_pos; typ[as.character(cur_id)] = cna$COPY_CHANGE[i]
}
}

#TP53_or_LOH17p_id = TP53_or_LOH17p_id[order(typ,decreasing=T)]
TP53_or_LOH17p_id = TP53_or_LOH17p_id[order(pos,decreasing=F)]

TP53_or_LOH17p_id = TP53_or_LOH17p_id[order(is.element(TP53_or_LOH17p_id,TP53_id),decreasing=T)]
TP53_or_LOH17p_id = TP53_or_LOH17p_id[order(is.element(TP53_or_LOH17p_id,LOH17p_id),decreasing=T)]

n_tp53 = rep(0,length(TP53_or_LOH17p_id))
names(n_tp53) = TP53_or_LOH17p_id
for(i in 1:length(n_tp53)){
n_tp53[i] = sum(mut$id[which(mut$Gene.refGene == "TP53")] == TP53_or_LOH17p_id[i])
}
for(i in 1:length(n_tp53)){
if(is.element(names(n_tp53)[i],setdiff(tp53_ddpcr_id,unique(mut$id[which(mut$Gene.refGene == "TP53")])))){
	n_tp53[i] = n_tp53[i] + 1
}
}

n_cna = rep(0,length(TP53_or_LOH17p_id))
names(n_cna) = TP53_or_LOH17p_id
for(i in 1:length(n_cna)){
n_cna[i] = sum(cna$id == TP53_or_LOH17p_id[i])
}


xmag = 1; ymag = 1

x_wid = sum(arminfo$V2[-23])/1000000/10
y_wid = length(TP53_or_LOH17p_id)


## import cbc data for all cases in TP53_or_LOH17p_id ##
## please access to clinical data provided by BBJ ##
cbc = as.data.frame(matrix(0,ncol=4,nrow=length(TP53_or_LOH17p_id))) ## data.frame of Hb, Ht, WBC, and PLT for each subject
rownames(cbc) = TP53_or_LOH17p_id; colnames(cbc) = c("Hb","Ht","WBC","PLT") ## data.frame of Hb, Ht, WBC, and PLT for each subject
gender = rep(1,length((TP53_or_LOH17p_id)))

Hb_l_thld = 11; Hb_h_m_thld = 16.5; Hb_h_f_thld = 16
Ht_h_thld = 49
WBC_l_thld = 3000; WBC_h_thld = 10000
Plt_l_thld = 10; Plt_h_thld = 45

cbc$Hb[which(cbc$Hb <= Hb_l_thld)] = -1
cbc$Hb[which(cbc$Hb >= Hb_h_m_thld & gender == 1)] = 1
cbc$Hb[which(cbc$Hb >= Hb_h_f_thld & gender == 2)] = 1
cbc$Hb[which(cbc$Hb != 1 & cbc$Hb != -1)] = 0
cbc$Hb[which(is.na(cbc$Hb))] = 2

cbc$Ht[which(cbc$Ht >= Ht_h_thld)] = 1
cbc$Ht[which(cbc$Ht != 1)] = 0
cbc$Ht[which(is.na(cbc$Ht))] = 2

cbc$WBC[which(cbc$WBC <= WBC_l_thld)] = -1
cbc$WBC[which(cbc$WBC >= WBC_h_thld)] = 1
cbc$WBC[which(cbc$WBC != 1 & cbc$WBC != -1)] = 0
cbc$WBC[which(is.na(cbc$WBC))] = 2

cbc$PLT[which(cbc$PLT <= Plt_l_thld)] = -1
cbc$PLT[which(cbc$PLT >= Plt_h_thld)] = 1
cbc$PLT[which(cbc$PLT!= 1 & cbc$PLT!= -1)] = 0
cbc$PLT[which(is.na(cbc$PLT))] = 2



pdf("TP53_subjects_copy_num_profile.pdf")

plot(NULL,NULL,
xlim=c(-50,x_wid + 50),
ylim=c(-10-y_wid,10),
axes=FALSE,
xlab="",ylab=""
) 

text(-15,1.5,"TP53 mutation",srt=-45,cex=0.4,adj=1)
text(-6,1.5,"17p alteration",srt=-45,cex=0.4,adj=1)


for(i in 1:length(TP53_or_LOH17p_id)){

# text(-20,-(i-1/2),i,cex=0.15)

if(i%%2 == 1){
rect(0*xmag,-(i-1)*ymag,x_wid*xmag,-i*ymag,col="gray90",border="transparent")
}

if(n_tp53[i]==1){
polygon(c(-16,-11,-11,-16),c(-i,-i,-(i-1),-(i-1)),col="orange",border="transparent")
}else if(n_tp53[i]>=2){
polygon(c(-16,-11,-11,-16),c(-i,-i,-(i-1),-(i-1)),col="orange4",border="transparent")
}

#polygon(c(-19,-17,-17,-19),c(-i,-i,-(i-1),-(i-1)),col=ifelse(cbc[i,"Ht"]==1,"firebrick2",ifelse(cbc[i,"Ht"]==0,"grey70","grey90")),border="transparent")
#polygon(c(-21,-19,-19,-21),c(-i,-i,-(i-1),-(i-1)),col=ifelse(cbc[i,"Hb"]==1,"firebrick2",ifelse(cbc[i,"Hb"]==-1,"skyblue",ifelse(cbc[i,"Hb"]==0,"grey70","grey90"))),border="transparent")
#polygon(c(-23,-21,-21,-23),c(-i,-i,-(i-1),-(i-1)),col=ifelse(cbc[i,"PLT"]==1,"firebrick2",ifelse(cbc[i,"PLT"]==-1,"skyblue",ifelse(cbc[i,"PLT"]==0,"grey70","grey90"))),border="transparent")
#polygon(c(-25,-23,-23,-25),c(-i,-i,-(i-1),-(i-1)),col=ifelse(cbc[i,"WBC"]==1,"firebrick2",ifelse(cbc[i,"WBC"]==-1,"skyblue",ifelse(cbc[i,"WBC"]==0,"grey70","grey90"))),border="transparent")

cur_cna = cna[which(cna$id == TP53_or_LOH17p_id[i]),]
if(nrow(cur_cna) == 0){ next }

for(j in 1:nrow(cur_cna)){
	
	chr_pos =  cum_chr_pos[as.numeric(cur_cna$chr[j])]
	stt_pos = as.numeric(cur_cna$start[j])/10 + chr_pos
	end_pos = as.numeric(cur_cna$end[j])/10 + chr_pos
	type_col = cur_cna$COLOR[j]
	
	rect(stt_pos*xmag,-(i-1)*ymag,end_pos*xmag,-i*ymag,col=type_col,border="transparent")
	
	if((as.numeric(cur_cna$chr[j]) == 17) & ((as.numeric(cur_cna$end[j])-7571720/1000000)*(as.numeric(cur_cna$start[j])-7590868/1000000)<=0)){
		polygon(c(-10,-5,-5,-10),c(-i,-i,-(i-1),-(i-1)),col=type_col,border="transparent")	
	}
}

col_cna = "gray"; if(n_cna[i]>=3){col_cna = "violet"}else{col_cna = "skyblue"}
rect(x_wid*xmag*1.01,-(i-1)*ymag,x_wid*xmag*1.01+n_cna[i]*3.5,-i*ymag,col=col_cna,border="transparent")
}

segments(x_wid*xmag*1.01,-(length(TP53_or_LOH17p_id)-1)*ymag,x_wid*xmag*1.01+n_cna[i]*3.5,-length(TP53_or_LOH17p_id)*ymag)


for(i in 1:22){
	text((cum_chr_pos[i]+cum_chr_pos[i+1])/2,1.5+(i%%2)*1.5,i,cex=0.4)
	segments(cum_chr_pos[i+1],0,cum_chr_pos[i+1],-length(TP53_or_LOH17p_id),col="white")
}

rect(x_wid + 10,-y_wid+9,x_wid + 20,-y_wid+8,col="orange",border="transparent"); text(x_wid + 30,-y_wid+8.5,"TP53 mutation",cex=0.3,adj=0)
rect(x_wid + 10,-y_wid+7,x_wid + 20,-y_wid+6,col="firebrick3",border="transparent"); text(x_wid + 30,-y_wid+6.5,"Duplication",cex=0.3,adj=0)
rect(x_wid + 10,-y_wid+5,x_wid + 20,-y_wid+4,col="dodgerblue3",border="transparent"); text(x_wid + 30,-y_wid+4.5,"Deletion",cex=0.3,adj=0)
rect(x_wid + 10,-y_wid+3,x_wid + 20,-y_wid+2,col="chartreuse3",border="transparent"); text(x_wid + 30,-y_wid+2.5,"UPD",cex=0.3,adj=0)
rect(x_wid + 10,-y_wid+0,x_wid + 20,-y_wid+1,col="tan",border="transparent"); text(x_wid + 30,-y_wid+0.5,"Unclassifiable",cex=0.3,adj=0)

dev.off()



