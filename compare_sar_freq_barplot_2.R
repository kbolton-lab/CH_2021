
#
cf_thld = 0.0

# get info of chr length
arminfo = read.table("cytoBand/arm_chr_length.txt")
cent = arminfo[1:22,4]

# sar
sar_upd_all_file = "part/upd_all_for_part.txt_output.txt"
sar_gain_all_file = "part/gain_all_for_part.txt_output.txt"
sar_loss_all_file = "part/loss_all_for_part.txt_output.txt"

sar_upd_all = as.matrix(read.table(sar_upd_all_file,header=F,quote="",sep="\t"))
sar_gain_all = as.matrix(read.table(sar_gain_all_file,header=F,quote="",sep="\t"))
sar_loss_all = as.matrix(read.table(sar_loss_all_file,header=F,quote="",sep="\t"))

sar = rbind(
cbind(sar_upd_all,rep(0,nrow(sar_upd_all))),
cbind(sar_gain_all,rep(1,nrow(sar_gain_all))),
cbind(sar_loss_all,rep(-1,nrow(sar_loss_all)))
)

labels = sar[,4]
sar = sar[,-4]
sar = apply(sar,c(1,2),as.numeric)

labels = c(

paste0(c(
"11p","12q","13q","14q",
"15q","16p","16q","22q",
"1p","1q","4q","6p",
"9p","9q","11q",
"17q"),"UPD"),

paste0("+",
c("3q","12","14q",
"8","15","18","21","22")),

paste0("del(",c(
"4q23-24","10q25-q26","14q24",
"14q32","16p11","2p23","3p13",
"5q14-32","6q16-24","7q32-36","8p23",
"9q31","11q14-23","13q13-31","14q11",
"17p13-11","17q11","20q11-13","21q22",
"22q12"),
")")

)
rownames(sar) = labels 


cna_path = "MatMed_2021_cna_table.csv"  
cna_bbj = read.csv(cna_path,header=T,stringsAsFactor=F,quote="",colClass="character")
cna_bbj = cna_bbj[which(!is.na(cna_bbj$chr)),]
cna_bbj = cna_bbj[which(as.numeric(cna_bbj$CELL_FRAC)>=cf_thld),]

## import clinical imformation ##

## non-HM cases
id_path = ~~  ## 10623 randomly selected cases
id_bbj = unlist(read.table(id_path,header=F,stringsAsFactor=F,quote="",sep="\t",comment.char=""))

## age and gender
sex_bbj = ## sex
age_bbj = ## age

id_bbj = id_bbj[which(age_bbj>=60 & age_bbj<=75)]

cna_bbj = cna_bbj[which(is.element(cna_bbj$id,id_bbj)),]

label_cna_bbj = paste0(sort(rep(1:22,4*2)),":",sort(rep(c("CNN-LOH","Duplication","Deletion","Unknown"),2)),":",c("p","q"))

mat_cna_bbj =  matrix(0,ncol=nrow(sar),nrow=length(id_bbj))
colnames(mat_cna_bbj) = rownames(sar); rownames(mat_cna_bbj) = id_bbj

for(i in 1:nrow(cna_bbj)){

	cur_id = cna_bbj$id[i]
	if(!is.element(cur_id,id_bbj)) next

	cur_type = NA
	if(cna_bbj$COPY_CHANGE[i] == "loss"){
		cur_type = -1
	}else if(cna_bbj$COPY_CHANGE[i] == "neutral"){
		cur_type = 0
	}else if(cna_bbj$COPY_CHANGE[i] == "gain"){
		cur_type = 1
	}else{
		cur_type = NA
	}
	if(is.na(cur_type)) next 
	
	xx = which(
	as.numeric(cna_bbj$chr[i]) == sar[,1] & 
	(as.numeric(cna_bbj$start[i])*(10^6) - sar[,3])*(as.numeric(cna_bbj$end[i])*(10^6) - sar[,2]) <= 0 & 
	cur_type == sar[,4]
	)
	
	mat_cna_bbj[cur_id,xx] = max(cna_bbj$CELL_FRAC[i],mat_cna_bbj[cur_id,xx])
}

sort(apply(mat_cna_bbj!=0,2,sum))


jacobs_path = paste0(path,"previous_study/Jacobs_2012_table.txt")
cna_jacobs = read.table(jacobs_path,header=T,stringsAsFactor=F,quote="",sep="\t",comment.char="")

cna_jacobs$type = cna_jacobs$STATE
cna_jacobs$type[which(cna_jacobs$type=="NEUTRAL")] = "CN-LOH"
cna_jacobs$type[which(cna_jacobs$type=="GAIN")] = "Gain"
cna_jacobs$type[which(cna_jacobs$type=="LOSS")] = "Loss"
cna_jacobs$cf = cna_jacobs$PROPORTION_ABNORMAL
cna_jacobs$CHR = cna_jacobs$CHROM
cna_jacobs$beg_GRCh37 = as.numeric(cna_jacobs$SEG_START)
cna_jacobs$end_GRCh37 = as.numeric(cna_jacobs$SEG_END)
cna_jacobs$START_MB =  as.numeric(cna_jacobs$SEG_START)/1000000
cna_jacobs$END_MB = as.numeric(cna_jacobs$SEG_END)/1000000
cna_jacobs$sample_id = cna_jacobs$ID
cna_jacobs$id = paste0("jacobs_",cna_jacobs$ID)

cna_jacobs$COPY_CHANGE = cna_jacobs$type

cna_jacobs$CELL_FRAC = cna_jacobs$cf

cna_jacobs = cna_jacobs[which(cna_jacobs$CHR!="X"),]
cna_jacobs = cna_jacobs[which(as.numeric(cna_jacobs$CELL_FRAC)>=cf_thld),]

cna_jacobs_status = matrix(NA,ncol=2,nrow=max(cna_jacobs$ID))

for(i in 1:nrow(cna_jacobs)){
	cna_jacobs_status[cna_jacobs$ID[i],1] = cna_jacobs$CANCER_CATEGORY[i]
	cna_jacobs_status[cna_jacobs$ID[i],2] = cna_jacobs$AGE5_DNA[i]
}
cna_jacobs_status[which(cna_jacobs_status[,1]!="cancer-free"),1] = "cancer"

c_1 = table(cna_jacobs_status[which(cna_jacobs_status[,1]!="cancer"),2])["60-65"]/0.008
c_2 = table(cna_jacobs_status[which(cna_jacobs_status[,1]!="cancer"),2])["65-70"]/0.0113
c_3 = table(cna_jacobs_status[which(cna_jacobs_status[,1]!="cancer"),2])["70-75"]/0.0166

f_1 = table(cna_jacobs_status[which(cna_jacobs_status[,1]!="cancer-free"),2])["60-65"]/0.006
f_2 = table(cna_jacobs_status[which(cna_jacobs_status[,1]!="cancer-free"),2])["65-70"]/0.0075
f_3 = table(cna_jacobs_status[which(cna_jacobs_status[,1]!="cancer-free"),2])["70-75"]/0.0133

cna_jacobs = cna_jacobs[which(is.element(cna_jacobs$AGE5_DNA,c("60-65","65-70","70-75"))),]
id_jacobs = paste0("jacobs_",1:(c_1+c_2+c_3+f_1+f_2+f_3))

mat_cna_jacobs =  matrix(0,ncol=nrow(sar),nrow=length(id_jacobs))
colnames(mat_cna_jacobs) = rownames(sar); rownames(mat_cna_jacobs) = id_jacobs
for(i in 1:nrow(cna_jacobs)){

	cur_id = cna_jacobs$id[i]
	if(!is.element(cur_id,id_jacobs)) next
	
	cur_type = "unknown"
	if(cna_jacobs[i,"COPY_CHANGE"] == "CN-LOH"){
		cur_type = 0
	}else if(cna_jacobs[i,"COPY_CHANGE"] == "Gain"){
		cur_type = 1
	}else if(cna_jacobs[i,"COPY_CHANGE"] == "Loss"){
		cur_type = -1
	}else{
		cur_type = NA
	}
	if(is.na(cur_type)) next 

	xx = which(
	as.numeric(as.numeric(cna_jacobs[i,"CHR"])) == sar[,1] & 
	(as.numeric(cna_jacobs[i,"START_MB"])*(10^6) - sar[,3])*(as.numeric(cna_jacobs[i,"END_MB"])*(10^6) - sar[,2]) <= 0 & 
	cur_type == sar[,4]
	)
	
	mat_cna_jacobs[cur_id,xx] = max(cna_jacobs$cf[i],mat_cna_jacobs[cur_id,xx])
}
sort(apply(mat_cna_jacobs!=0,2,sum))


laurie_path = paste0(path,"previous_study/Laurie_2012_table.txt")
cna_laurie = read.table(laurie_path,header=T,stringsAsFactor=F,quote="",sep="\t",comment.char="")
cna_laurie = cna_laurie[which(cna_laurie$age.tissue.collect >= 60 & cna_laurie$age.tissue.collect <= 75),]

cna_laurie$CHR = cna_laurie$chromosome

cna_laurie$type = cna_laurie$mosaic.type
cna_laurie$type[which(cna_laurie$type=="aupd")] = "CN-LOH"
cna_laurie$type[which(cna_laurie$type=="gain")] = "Gain"
cna_laurie$type[which(cna_laurie$type=="loss")] = "Loss"

cna_laurie$COPY_CHANGE = cna_laurie$type

cna_laurie$muDiff = 2*as.numeric(cna_laurie$anom.baf.dev.med)
cna_laurie$cf = 
(cna_laurie$type == "Loss")*(2*cna_laurie$muDiff/(1+cna_laurie$muDiff))+
(cna_laurie$type == "Gain")*(2*cna_laurie$muDiff/(1-cna_laurie$muDiff))+
(cna_laurie$type == "CN-LOH")*cna_laurie$muDiff

cna_laurie$START_MB = as.numeric(cna_laurie$Start..GRCh37.)/(10^6)
cna_laurie$END_MB = as.numeric(cna_laurie$End..GRCh37.)/(10^6)
cna_laurie$sample_id = cna_laurie$subject.id
cna_laurie$ID = cna_laurie$subject.id
cna_laurie$CELL_FRAC = cna_laurie$cf

cna_laurie = cna_laurie[which(cna_laurie$CHR!="X"),]
cna_laurie = cna_laurie[which(as.numeric(cna_laurie$CELL_FRAC)>=cf_thld),]

cna_laurie = cna_laurie[which(!is.na(cna_laurie$age.tissue.collect)),]

id_laurie = paste0("laurie_",1:(3772+3845+3607))

label_cna_laurie = paste0(sort(rep(1:22,4*2)),":",sort(rep(c("CNN-LOH","Duplication","Deletion","Unknown"),2)),":",c("p","q"))
mat_cna_laurie =  matrix(0,ncol=nrow(sar),nrow=length(id_laurie))
colnames(mat_cna_laurie) = rownames(sar); rownames(mat_cna_laurie) = id_laurie
for(i in 1:nrow(cna_laurie)){
	
	cur_id = paste0("laurie_",cna_laurie[i,"ID"])
	if(!is.element(cur_id,id_laurie)) next
	
	cur_type = NA
	if(cna_laurie[i,"COPY_CHANGE"] == "CN-LOH"){
		cur_type = 0
	}else if(cna_laurie[i,"COPY_CHANGE"] == "Gain"){
		cur_type = 1
	}else if(cna_laurie[i,"COPY_CHANGE"] == "Loss"){
		cur_type = -1
	}else{
		cur_type = NA
	}
	if(is.na(cur_type)) next
	
	xx = which(
	as.numeric(cna_laurie[i,"CHR"]) == sar[,1] & 
	(as.numeric(cna_laurie[i,"START_MB"])*(10^6) - sar[,3])*(as.numeric(cna_laurie[i,"END_MB"])*(10^6) - sar[,2]) <= 0 & 
	cur_type == sar[,4]
	)
	
	mat_cna_laurie[cur_id,xx] = max(cna_laurie$cf[i],mat_cna_laurie[cur_id,xx])
}
sort(apply(mat_cna_laurie!=0,2,sum))


loh_path =  paste0(path,"previous_study/Loh_2020.txt")
cna_loh_2 = read.table(loh_path,header=T,stringsAsFactor=F,quote="",sep="\t",comment.char="")
cna_loh = cna_loh_2[which(cna_loh_2$CHR!="X"),]
cna_loh = cna_loh[which(is.element(cna_loh$AGE,c("60-65","65-70","70-75"))),]


m1 = round(length(unique(cna_loh_2[which(cna_loh_2$AGE == "60-65" & cna_loh_2$SEX == "M"),"ID"]))/0.047)
f1 = round(length(unique(cna_loh_2[which(cna_loh_2$AGE == "60-65" & cna_loh_2$SEX == "F"),"ID"]))/0.04)
m2 = round(length(unique(cna_loh_2[which((cna_loh_2$AGE == "65-70" | cna_loh_2$AGE == "70-75") & cna_loh_2$SEX == "M"),"ID"]))/0.060)
f2 = round(length(unique(cna_loh_2[which((cna_loh_2$AGE == "65-70" | cna_loh_2$AGE == "70-75") & cna_loh_2$SEX == "F"),"ID"]))/0.049)
n_loh_1 = m1 + f1
n_loh_2 = m2 + f2
n_loh = n_loh_1 + n_loh_2

id_loh = paste0("loh_",1:n_loh)

label_cna_loh = paste0(sort(rep(1:22,4*2)),":",sort(rep(c("CNN-LOH","Duplication","Deletion","Unknown"),2)),":",c("p","q"))
mat_cna_loh =  matrix(0,ncol=nrow(sar),nrow=length(id_loh))
colnames(mat_cna_loh) = rownames(sar); rownames(mat_cna_loh) = id_loh

for(i in 1:nrow(cna_loh)){

	cur_id = paste0("loh_",cna_loh[i,"ID"])
	if(!is.element(cur_id,id_loh)) next

	cur_type = NA
	if(cna_loh[i,"COPY_CHANGE"] == "neutral"){
		cur_type = 0
	}else if(cna_loh[i,"COPY_CHANGE"] == "gain"){
		cur_type = 1
	}else if(cna_loh[i,"COPY_CHANGE"] == "loss"){
		cur_type = -1
	}else{
		cur_type = NA
	}
	if(is.na(cur_type)) next
	
	xx = which(
	as.numeric(cna_loh[i,"CHR"]) == sar[,1] & 
	(as.numeric(cna_loh[i,"START_MB"])*(10^6) - sar[,3])*(as.numeric(cna_loh[i,"END_MB"])*(10^6) - sar[,2]) <= 0 & 
	cur_type == sar[,4]
	)
	
	mat_cna_loh[cur_id,xx] = max(as.numeric(cna_loh$CELL_FRAC[i]),mat_cna_loh[cur_id,xx])
}
sort(apply(mat_cna_loh!=0,2,sum))



## make frequency table ##
## stratified by cell fractions ##
thld1 = .05
thld2 = .00001
thld3 = .02
thld4 = .001

bbj = cbind(
apply((mat_cna_bbj>=thld1),2,sum)/nrow(mat_cna_bbj),
apply((mat_cna_bbj>=thld2 & mat_cna_bbj<thld1),2,sum)/nrow(mat_cna_bbj),
apply((mat_cna_bbj>=thld3)&(mat_cna_bbj<thld2),2,sum)/nrow(mat_cna_bbj),
apply((mat_cna_bbj>=thld4)&(mat_cna_bbj<thld3),2,sum)/nrow(mat_cna_bbj)
)

jacobs = cbind(
apply(mat_cna_jacobs>=thld1,2,sum)/nrow(mat_cna_jacobs),
apply((mat_cna_jacobs>=thld2 & mat_cna_jacobs<thld1),2,sum)/nrow(mat_cna_jacobs),
apply(mat_cna_jacobs>=thld3,2,sum)/nrow(mat_cna_jacobs),
apply(mat_cna_jacobs>=thld4,2,sum)/nrow(mat_cna_jacobs)
)

laurie = cbind(
apply(mat_cna_laurie>=thld1,2,sum)/nrow(mat_cna_laurie),
apply((mat_cna_laurie>=thld2 & mat_cna_laurie<thld1),2,sum)/nrow(mat_cna_laurie),
apply(mat_cna_laurie>=thld3,2,sum)/nrow(mat_cna_laurie),
apply(mat_cna_laurie>=thld4,2,sum)/nrow(mat_cna_laurie)
)

loh = cbind(
apply(mat_cna_loh>=thld1,2,sum)/nrow(mat_cna_loh),
apply((mat_cna_loh>=thld2 & mat_cna_loh<thld1),2,sum)/nrow(mat_cna_loh),
apply(mat_cna_loh>=thld3,2,sum)/nrow(mat_cna_loh),
apply(mat_cna_loh>=thld4,2,sum)/nrow(mat_cna_loh)
)

col = sar[order(bbj[,4],decreasing=T),4]
jacobs = jacobs[order(bbj[,4],decreasing=T),]
laurie = laurie[order(bbj[,4],decreasing=T),]
loh = loh[order(bbj[,4],decreasing=T),]
bbj = bbj[order(bbj[,4],decreasing=T),]

mat = round(cbind(bbj,loh,laurie,jacobs)*100,3)

col = col[order(mat[,2],decreasing=T)]
col2 = rep(NA,length(col)); names(col2) = rownames(mat)
for(i in 1:length(col)){
if(col[i]==0){
col2[i] = "chartreuse3"
}else if(col[i]==1){
col2[i] = "firebrick3"
}else if(col[i]==-1){
col2[i] = "dodgerblue3"
}
}

mat = mat[order(mat[,1],decreasing=T),]

incl = c(
"del(20q11-13)","del(14q11)","+15","del(2p23)","1pUPD","del(9q31)","14qUPD",
"del(13q13-31)","del(5q14-32)","del(6q16-24)","+14q","9pUPD","+21","1qUPD",
"11qUPD","16pUPD","12qUPD","del(11q14-23)","del(7q32-36)","del(3p13)","del(22q12)",
"9qUPD","+22","22qUPD","del(17p13-11)","17qUPD","6pUPD","+18","11pUPD","+12","del(21q22)",
"4qUPD","16pUPD","13qUPD","15qUPD","+8","+3q","del(10q25-q26)"
)

col2 = col2[incl]
mat = mat[incl,]
mat = mat + 0.0002


## Draw ED Fig.4d ##
## cell fraction < 5% ##
pdf("barplot_1.pdf",width=5,height=5)

plot(
NULL,NULL,
xlab="",ylab="",
axes=F,
xlim=c(-5,50),
ylim=c(-0.2,2)
)

segments(-1,0,-1,0.7,lwd=.5)
lat = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7)
lat2= 0.01*(1:70)
for(i in 1:length(lat)){
segments(-1,lat[i],-1-0.2,lat[i],lwd=.5)
text(-1.4,lat[i],lat[i],adj=1,cex=0.15,font=2)
}
for(i in 1:length(lat2)){
segments(-1,lat2[i],-1-0.1,lat2[i],lwd=.5)
}

col_001_bbj = "white" #"thistle1" #"grey65"
col_01_bbj =  "lightpink"  #"white" #"grey45"
col_05_bbj = "firebrick2" #"grey25"
col_1_bbj = "firebrick4" #"grey5"

col_001_loh = "white" #"darkseagreen1" #"grey65"
col_01_loh = "palegreen" #"grey45"
col_05_loh = "palegreen3" #"grey25"
col_1_loh = "forestgreen" #"grey5"

col_001_laurie = "white" #"lightcyan" #"grey65"
col_01_laurie = "lightskyblue1" #"grey45"
col_05_laurie = "royalblue2" #"grey25"
col_1_laurie = "royalblue4" #"grey5"

col_001_jacobs = "white" #"burlywood1" #"grey65"
col_01_jacobs = "papayawhip" #"grey45"
col_05_jacobs = "navajowhite2" #"grey25"
col_1_jacobs = "navajowhite4" #"grey5"

for(i in 1:nrow(mat)){
	rect(i-1,0,i-4/5,mat[i,2],col=col_05_bbj,border="transparent")
	rect(i-4/5,0,i-3/5,mat[i,6],col=col_05_loh,border="transparent")
	rect(i-3/5,0,i-2/5,mat[i,10],col=col_05_laurie,border="transparent")
	rect(i-2/5,0,i-1/5,mat[i,14],col=col_05_jacobs,border="transparent")
	
	text(i-3/5,-1/100,rownames(mat)[i],srt=45,adj=1,cex=0.15,font=2)
}


i=i+1
rect(i-1,0,i-4/5,0.1,col=col_1_bbj,border="transparent")
rect(i-1,0.1,i-4/5,0.2,col=col_05_bbj,border="transparent")
rect(i-1,0.2,i-4/5,0.3,col=col_01_bbj,border="transparent")
rect(i-1,0.3,i-4/5,0.4,col=col_001_bbj,border="transparent")
	
rect(i-4/5,0,i-3/5,0.1,col=col_1_loh,border="transparent")
rect(i-4/5,0.1,i-3/5,0.2,col=col_05_loh,border="transparent")
rect(i-4/5,0.2,i-3/5,0.3,col=col_01_loh,border="transparent")
rect(i-4/5,0.3,i-3/5,0.4,col=col_001_loh,border="transparent")

rect(i-3/5,0,i-2/5,0.1,col=col_1_laurie,border="transparent")
rect(i-3/5,0.1,i-2/5,0.2,col=col_05_laurie,border="transparent")
rect(i-3/5,0.2,i-2/5,0.3,col=col_01_laurie,border="transparent")
rect(i-3/5,0.3,i-2/5,0.4,col=col_001_laurie,border="transparent")

rect(i-2/5,0,i-1/5,0.1,col=col_1_jacobs,border="transparent")
rect(i-2/5,0.1,i-1/5,0.2,col=col_05_jacobs,border="transparent")
rect(i-2/5,0.2,i-1/5,0.3,col=col_01_jacobs,border="transparent")
rect(i-2/5,0.3,i-1/5,0.4,col=col_001_jacobs,border="transparent")

dev.off()


## Draw ED Fig.4e ##
## cell fraction >= 5% ##
pdf("barplot_2.pdf",width=5,height=5)

plot(
NULL,NULL,
xlab="",ylab="",
axes=F,
xlim=c(-5,50),
ylim=c(-0.2,0.8)
)

segments(-1,0,-1,0.7,lwd=.5)
lat = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7)
#lat = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5)
lat2= 0.01*(1:70)
for(i in 1:length(lat)){
segments(-1,lat[i],-1-0.2,lat[i],lwd=.5)
text(-1.4,lat[i],lat[i],adj=1,cex=0.15,font=2)
}
for(i in 1:length(lat2)){
segments(-1,lat2[i],-1-0.1,lat2[i],lwd=.5)
}

col_001_bbj = "white" #"thistle1" #"grey65"
col_01_bbj =  "lightpink"  #"white" #"grey45"
col_05_bbj = "firebrick2" #"grey25"
col_1_bbj = "firebrick2"#"firebrick4" #"grey5"

col_001_loh = "white" #"darkseagreen1" #"grey65"
col_01_loh = "palegreen" #"grey45"
col_05_loh = "palegreen3" #"grey25"
col_1_loh = "palegreen3"#"forestgreen" #"grey5"

col_001_laurie = "white" #"lightcyan" #"grey65"
col_01_laurie = "lightskyblue1" #"grey45"
col_05_laurie = "royalblue2" #"grey25"
col_1_laurie = "royalblue2" #"grey5"

col_001_jacobs = "white" #"burlywood1" #"grey65"
col_01_jacobs = "papayawhip" #"grey45"
col_05_jacobs = "navajowhite2" #"grey25"
col_1_jacobs = "navajowhite2" #"grey5"

for(i in 1:nrow(mat)){
	rect(i-1,0,i-4/5,mat[i,1],col=col_1_bbj,border="transparent")
	rect(i-4/5,0,i-3/5,mat[i,5],col=col_1_loh,border="transparent")
	rect(i-3/5,0,i-2/5,mat[i,9],col=col_1_laurie,border="transparent")
	rect(i-2/5,0,i-1/5,mat[i,13],col=col_1_jacobs,border="transparent")
	text(i-3/5,-1/100,rownames(mat)[i],srt=45,adj=1,cex=0.15,font=2)
}

i=i+1
rect(i-1,0,i-4/5,0.1,col=col_1_bbj,border="transparent")
rect(i-1,0.1,i-4/5,0.2,col=col_05_bbj,border="transparent")
rect(i-1,0.2,i-4/5,0.3,col=col_01_bbj,border="transparent")
rect(i-1,0.3,i-4/5,0.4,col=col_001_bbj,border="transparent")
	
rect(i-4/5,0,i-3/5,0.1,col=col_1_loh,border="transparent")
rect(i-4/5,0.1,i-3/5,0.2,col=col_05_loh,border="transparent")
rect(i-4/5,0.2,i-3/5,0.3,col=col_01_loh,border="transparent")
rect(i-4/5,0.3,i-3/5,0.4,col=col_001_loh,border="transparent")

rect(i-3/5,0,i-2/5,0.1,col=col_1_laurie,border="transparent")
rect(i-3/5,0.1,i-2/5,0.2,col=col_05_laurie,border="transparent")
rect(i-3/5,0.2,i-2/5,0.3,col=col_01_laurie,border="transparent")
rect(i-3/5,0.3,i-2/5,0.4,col=col_001_laurie,border="transparent")

rect(i-2/5,0,i-1/5,0.1,col=col_1_jacobs,border="transparent")
rect(i-2/5,0.1,i-1/5,0.2,col=col_05_jacobs,border="transparent")
rect(i-2/5,0.2,i-1/5,0.3,col=col_01_jacobs,border="transparent")
rect(i-2/5,0.3,i-1/5,0.4,col=col_001_jacobs,border="transparent")

dev.off()

