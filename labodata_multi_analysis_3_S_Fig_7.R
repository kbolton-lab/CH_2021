

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
label_cna = paste0(sort(rep(1:22,3*2)),":",sort(rep(c("CNN-LOH","Duplication","Deletion"),2)),":",c("p","q"))
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
	
	if(! is.element(cur_id,id)) next
	
	n_cna[cur_id] = n_cna[cur_id] + 1
	
	if(cur_type == "Unknown") next
	
	cur_cent = cent[as.numeric(cna$chr[i])]
	if(cna$start[i]*(10^6) <= cur_cent & cna$end[i]*(10^6) <= cur_cent){
		cur_cna = paste0(cur_cna,":p")
	}else if(cna$end[i]*(10^6) >= cur_cent & cna$start[i]*(10^6) >= cur_cent){
		cur_cna = paste0(cur_cna,":q")
	}else if(cna$start[i]*(10^6) <= cur_cent & cna$end[i]*(10^6) >= cur_cent)
	if(! is.element(cur_cna,label_cna)){
		next
	}
	
	mat_cna[cur_id,cur_cna] = max(mat_cna[cur_id,cur_cna],cna$CELL_FRAC[i])
}


all_mat_mut = mat_mut
all_mat_cna = mat_cna

is_mut = as.numeric(apply(all_mat_mut!=0,1,sum)!=0)
is_cna = as.numeric(apply(all_mat_cna!=0,1,sum)!=0)
is_all = as.numeric((is_mut + is_cna) != 0)
is_both = is_mut * is_cna
is_cis = xxx
is_dis = is_both - is_cis
is_none = (1-is_mut) * (1-is_cna)

all_mat_mut = all_mat_mut[,which(apply(all_mat_mut != 0,2,sum) >= 20)]
all_mat_cna = all_mat_cna[,which(apply(all_mat_cna != 0,2,sum) >= 15)]
all_mat_cna = all_mat_cna[,setdiff(1:ncol(all_mat_cna),grep("Unknown",colnames(all_mat_cna)))]


## import clinical data ##
## data.frame of blood test values indicated in Supplementary Fig.7 for 11,234 cases ##
all_dat = ~~ 


###################
### Correlation ###
###################

transform = function(clinical){
	library(car)
	for(i in 1:ncol(clinical)){
		nna = which(!is.na(clinical[,i]))
		v = na.omit(clinical[,i])
		v1 = powerTransform(v+0.0001)
		v2 = na.omit(bcPower(v+0.0001, v1$lambda))
		clinical[nna,i] = v2
	}
	return(clinical)
}


mat_p = mat_t = matrix(NA,nrow=length(lab_all),ncol=4+ncol(all_mat_mut)+ncol(all_mat_cna))

for(i in 1:length(lab_all)){

cur_dat = all_dat[,3*i]
if(sum(!is.na(cur_dat))<=1){next}

cur_dat_raw = cur_dat

if(type[i]==0){
next
}else if(type[i]==1){
cur_dat[which(cur_dat=="99")] = NA
if(length(table(cur_dat))<=1){next}
if(sum(!is.na(cur_dat))<=50){next}
}else if(type[i] == 2){
if(sum(!is.na(cur_dat))<=50){next}
cur_dat = transform(as.matrix(as.numeric(cur_dat),ncol=1))
}
cur_data = as.numeric(cur_dat)

qt = quantile(na.omit(cur_data))
ol_u = qt[4]+1.5*(qt[4]-qt[2])
ol_l = qt[2]-1.5*(qt[4]-qt[2])
cur_data[which(cur_data > ol_u)] = NA
cur_data[which(cur_data < ol_l)] = NA

###########
### mut ###
###########
cmp_data = data.frame(cur_dat,is_mut,all_age,all_sex)
cmp_data = cmp_data[which(is_mut==1 | is_none==1),]
s = summary(lm(cur_dat ~ ., data=cmp_data))
mat_p[i,1] = s$coef[2,4]; mat_t[i,1] = s$coef[2,1]/sd(na.omit(cur_dat))


###########
### cna ###
###########
cmp_data = data.frame(cur_dat,is_cna,all_age,all_sex,all_version)
cmp_data = cmp_data[which(is_cna==1 | is_none==1),]
s = summary(lm(cur_dat ~ ., data=cmp_data))
mat_p[i,2] = s$coef[2,4]; mat_t[i,2] = s$coef[2,1]/sd(na.omit(cur_dat))


###########
### cis ###
###########
cmp_data = data.frame(cur_dat,is_cis,all_age,all_sex,all_version)
cmp_data = cmp_data[which(is_cis==1 | is_none==1),]
s = summary(lm(cur_dat ~ ., data=cmp_data))
mat_p[i,3] = s$coef[2,4]; mat_t[i,3] = s$coef[2,1]/sd(na.omit(cur_dat))


###########
### dis ###
###########
cmp_data = data.frame(cur_dat,is_dis,all_age,all_sex,all_version)
cmp_data = cmp_data[which(is_dis==1 | is_none==1),]
s = summary(lm(cur_dat ~ ., data=cmp_data))
mat_p[i,4] = s$coef[2,4]; mat_t[i,4] = s$coef[2,1]/sd(na.omit(cur_dat))


######################
### individual SNV ###
######################
for(j in 1:ncol(all_mat_mut)){
cur_alt = as.numeric(all_mat_mut[,j]!=0)
cmp_data = data.frame(cur_dat,cur_alt,all_age,all_sex)
cmp_data = cmp_data[which(cur_alt==1 | is_none==1),]
s = summary(lm(cur_dat ~ ., data=cmp_data))
if(rownames(s$coefficients)[2] == "cur_alt"){
	cat(colnames(all_mat_mut)[j],"\t")
	cat(lab_all[i],"\t")
	cat(sum((all_mat_mut[,j]!=0)*(!is.na(cur_dat))),"\n")
	mat_p[i,j+4] = s$coef[2,4]; mat_t[i,j+4] = s$coef[2,1]/sd(na.omit(cur_dat))
}
}



######################
### individual CNA ###
######################
for(j in 1:ncol(all_mat_cna)){
cur_alt = as.numeric(all_mat_cna[,j]!=0)
cmp_data = data.frame(cur_dat,cur_alt,all_age,all_sex,all_version)
cmp_data = cmp_data[which(cur_alt==1 | is_none==1),]
s = summary(lm(cur_dat ~ ., data=cmp_data))
if(rownames(s$coefficients)[2] == "cur_alt"){
	cat(colnames(all_mat_cna)[j],"\t")
	cat(lab_all[i],"\t")
	cat(sum((all_mat_cna[,j]!=0)*(!is.na(cur_dat))),"\n")
	mat_p[i,j+ncol(all_mat_mut)+4] = s$coef[2,4]; mat_t[i,j+ncol(all_mat_mut)+4] = s$coef[2,1]/sd(na.omit(cur_dat))
}
}


mat_p = mat_p[,c(1:4,which(apply(all_mat_mut!=0,2,sum)>=30)+4,which(apply(all_mat_cna!=0,2,sum)>=30)+ncol(all_mat_mut)+4)]
mat_t = mat_t[,c(1:4,which(apply(all_mat_mut!=0,2,sum)>=30)+4,which(apply(all_mat_cna!=0,2,sum)>=30)+ncol(all_mat_mut)+4)]
mat_q = matrix(p.adjust(as.numeric(mat_q),"BH"),nrow=nrow(mat_q))



##############################
## draw Supplementary Fig.7 ##
##############################

getCnaLabel = function(s){
splited = unlist(strsplit(s,":"))
type=NA
prf = ifelse(splited[2]=="CNN-LOH","",ifelse(splited[2]=="Duplication","+",ifelse(splited[2]=="Deletion","del(","Unknown")))
suf = ifelse(splited[2]=="CNN-LOH","UPD",ifelse(splited[2]=="Duplication","",ifelse(splited[2]=="Deletion",")","Unknown")))
return(paste0(prf,splited[1],splited[3],suf))
}

xwid = nrow(mat_p)
ywid = ncol(mat_p)
xlab = lab_all[23:length(lab_all)]
ylab = c("Mut","CNA","Cis","Dis",colnames(all_mat_mut)[which(apply(all_mat_mut!=0,2,sum)>=30)],sapply(colnames(all_mat_cna)[which(apply(all_mat_cna!=0,2,sum)>=30)],getCnaLabel))
fonts = c(rep(1,8),rep(3,sum(apply(all_mat_mut!=0,2,sum)>=30)),rep(1,sum(apply(all_mat_cna!=0,2,sum)>=30)))

r = 0.15; theta = c(18,18+72,18+72*2,18+72*3,18+72*4)*pi/180

pdf("rect_cor_clinical.pdf",w=ywid/5,h=ywid/4.5)
plot(NULL,NULL,
xlim=c(-ywid/10,ywid + ywid/10),
#xlim=c(-xwid/10,xwid + xwid/10),
ylim=c(-ywid/10,ywid + ywid/10),
xlab="",ylab="",axes=F
)
rect(-.05,-.05,xwid+.05,ywid+.05,col="gray",border="transparent")

for(i in 1:xwid){
for(j in 1:ywid){
p = mat_p[i,j]; t = mat_t[i,j]

if(is.na(p)){s=0}else{
if(p < 10^(-10)){p = 10^(-10)}
s = min(.45,.45*abs(-log10(p))/abs(-log10(0.001)))
}
if(is.na(t)){c="gray"}else{
if(t>0){c = rgb(1,0,0,alpha=min(1,t))}
if(t<0){c = rgb(0,0,1,alpha=min(1,-t))}
}

rect(i-1,ywid-(j-1),i,ywid-j,col="gray",border="transparent")
rect(i-1+.05,ywid-(j-1+.05),i-.05,ywid-(j-.05),col="white",border="transparent")
rect(i-.5-s,ywid-(j-.5-s),i-.5+s,ywid-(j-.5+s),col=c,border="transparent")
text(i-.5,ywid+ywid/50,xlab[i],srt=45,adj=0,cex=0.7)
text(-xwid/100,ywid-(j-.5),ylab[j],adj=1,cex=0.7,font=fonts[j])

if(!is.na(mat_p[i,j]) &  mat_p[i,j] < 0.05 & mat_q[i,j] < 0.1){
x = min(1,abs(t)/1.3)
if(x<0.5){x = 0}else{x = 1}
segments(i-.5,ywid-(j-.5),(i-.5)+r*cos(theta),(ywid-(j-.5))+r*sin(theta),col=rgb(x,x,x,alpha=1),lwd=1.5)
}

}
}

s.05 = min(.45,.45*abs(-log10(.05))/abs(-log10(0.001)))
s.01 = min(.45,.45*abs(-log10(.01))/abs(-log10(0.001)))
s.001 = min(.45,.45*abs(-log10(.001))/abs(-log10(0.001)))

polygon(c(10,0,10),c(-5.5+s.001,-5.5,-5.5-s.001),col="gray",border="transparent")
segments(10*s.05/s.001,-5.5+s.05,10*s.05/s.001,-5.5-s.05,col="white")
segments(10*s.01/s.001,-5.5+s.01,10*s.01/s.001,-5.5-s.01,col="white")
segments(10*s.001/s.001,-5.5+s.001,10*s.001/s.001,-5.5-s.001,col="white")

for(i in 1:100){
t = seq(-1,1,length=100)[i]
if(t>0){c = rgb(1,0,0,alpha=min(1,t))}
if(t<0){c = rgb(0,0,1,alpha=min(1,-t))}
rect(3+(i-1)*3/100,-7.5,3+i*3/100,-7,col=c,border="transparent")
}

dev.off()

}


