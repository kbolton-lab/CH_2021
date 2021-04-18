
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

is_mut = as.numeric(apply(all_mat_mut!=0,1,sum) != 0)
is_cna = as.numeric(apply(all_mat_cna!=0,1,sum) != 0)
is_both = as.numeric(apply(all_mat_mut!=0,1,sum) != 0 & apply(all_mat_cna!=0,1,sum) != 0)

vaf = as.numeric(apply(all_mat_mut,1,function(x){max(as.numeric(x))}))
cf = as.numeric(apply(all_mat_cna,1,function(x){max(as.numeric(x))}))

all_n_mut = n_mut
all_n_cna = n_cna




### please import clinical data provided by BBJ ##
all_dat = ~~  ### matrix of WBC, Hb, MCV, MCHC, and PLT
all_age = ~~
all_sex = ~~
all_HM = ~~


## Examine correlation ##

## cox-box transformation##
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

mat_p = mat_t = matrix(NA,nrow=length(lab_all),ncol=12+ncol(all_mat_mut)+ncol(all_mat_cna))
colnames(mat_p) = colnames(mat_t) = c("SNV+","CNA+","nSNV","VAF","nCNA","CF","SNV+CNA",colnames(all_mat_mut),colnames(all_mat_cna))
rownames(mat_p) = rownames(mat_t) = lab_all

mat_cf_p = mat_cf_t = matrix(NA,nrow=length(lab_all),ncol=ncol(all_mat_mut)+ncol(all_mat_cna))
colnames(mat_cf_p) = colnames(mat_cf_t) = c(colnames(all_mat_mut),colnames(all_mat_cna))
rownames(mat_cf_p) = rownames(mat_cf_t) = lab_all

for(i in 1:length(lab_all)){

cur_dat = all_dat[,i]
if(sum(!is.na(cur_dat))<=1){next}

cur_dat_raw = cur_dat
cur_dat = transform(as.matrix(as.numeric(cur_dat),ncol=1))
cur_data = as.numeric(cur_dat)

qt = quantile(na.omit(cur_data))
ol_u = qt[4]+1.5*(qt[4]-qt[2])
ol_l = qt[2]-1.5*(qt[4]-qt[2])
print(ol_u)
print(ol_l)
print(max(na.omit(cur_data)))
print(min(na.omit(cur_data)))
cur_data[which(cur_data > ol_u)] = NA
cur_data[which(cur_data < ol_l)] = NA

###########
### mut ###
###########
cmp_data = data.frame(cur_dat,is_mut,all_age,all_sex,all_HM)
cmp_data = cmp_data[which(is_any==1 | is_none==1),]
s = summary(lm(cur_dat ~ ., data=cmp_data))
cat(lab_all[i])
print(s)
mat_p[i,1] = s$coef[2,4]; mat_t[i,1] = s$coef[2,1]/sd(na.omit(cur_dat))


###########
### cna ###
###########
cmp_data = data.frame(cur_dat,is_cna,all_age,all_sex,all_version,all_HM)
cmp_data = cmp_data[which(is_cna==1 | is_none==1),]
s = summary(lm(cur_dat ~ ., data=cmp_data))
cat(lab_all[i])
print(s)
mat_p[i,2] = s$coef[2,4]; mat_t[i,2] = s$coef[2,1]/sd(na.omit(cur_dat))


############
### mut_n ##
############
cmp_data = data.frame(cur_dat,all_n_mut,all_age,all_sex,all_HM)
cmp_data = cmp_data[which(is_mut==1),]
s = summary(lm(cur_dat ~ ., data=cmp_data))
cat(lab_all[i])
print(s)
mat_p[i,3] = s$coef[2,4]; mat_t[i,3] = s$coef[2,1]/sd(na.omit(cur_dat))


###############
### mut_vaf ###
###############
cmp_data = data.frame(cur_dat,vaf,all_age,all_sex,all_HM)
cmp_data = cmp_data[which(is_mut==1),]
s = summary(lm(cur_dat ~ ., data=cmp_data))
cat(lab_all[i])
print(s)
mat_p[i,4] = s$coef[2,4]; mat_t[i,4] = s$coef[2,1]/sd(na.omit(cur_dat))


#############
### cna_n ###
#############
cmp_data = data.frame(cur_dat,all_n_cna,is_cna,all_age,all_sex,all_version,all_HM)
cmp_data = cmp_data[which(is_cna==1),]
s = summary(lm(cur_dat ~ ., data=cmp_data))
cat(lab_all[i])
print(s)
mat_p[i,5] = s$coef[2,4]; mat_t[i,5] = s$coef[2,1]/sd(na.omit(cur_dat))
print(s$coef[2,])


##############
### cna_cf ###
##############
cmp_data = data.frame(cur_dat,cf,is_cna,all_age,all_sex,all_version,all_HM)
cmp_data = cmp_data[which(is_cna==1),]
s = summary(lm(cur_dat ~ ., data=cmp_data))
cat(lab_all[i])
print(s)
mat_p[i,6] = s$coef[2,4]; mat_t[i,6] = s$coef[2,1]/sd(na.omit(cur_dat))


###############
### SNV+CNA ###
###############
cmp_data = data.frame(cur_dat,is_both,all_age,all_sex,all_version,all_HM)
cmp_data = cmp_data[which(is_both==1 | is_none==1),]
s = summary(lm(cur_dat ~ ., data=cmp_data))
cat(lab_all[i])
print(s)
mat_p[i,7] = s$coef[2,4]; mat_t[i,12] = s$coef[2,1]/sd(na.omit(cur_dat))
print(s$coef[2,])


######################
### individual SNV ###
######################
for(j in 1:ncol(all_mat_mut)){
cur_alt = as.numeric(all_mat_mut[,j]!=0)
cur_vaf = as.numeric(all_mat_mut[,j])
cmp_data = data.frame(cur_dat,cur_alt,all_age,all_sex,all_HM)
cmp_data = cmp_data[which(cur_alt==1 | is_none==1),]
s = summary(lm(cur_dat ~ ., data=cmp_data))
cmp_data_2 = data.frame(cur_dat,cur_vaf,all_age,all_sex,all_HM)
cmp_data_2 = cmp_data_2[which(cur_vaf>0),]
s2 = summary(lm(cur_dat ~ ., data=cmp_data_2))
if(rownames(s$coefficients)[2] == "cur_alt"){
	cat(colnames(all_mat_mut)[j],"\t")
	cat(lab_all[i],"\t")
	cat(sum((all_mat_mut[,j]!=0)*(!is.na(cur_dat))),"\n")
	
	mat_p[i,j+7] = s$coef[2,4]; mat_t[i,j+7] = s$coef[2,1]/sd(na.omit(cur_dat))
	mat_cf_p[i,j] = s2$coef[2,4]; mat_cf_t[i,j] = s2$coef[2,1]/sd(na.omit(cur_dat))
	
	pdf(paste0(colnames(all_dat)[i]," - ",colnames(all_mat_mut)[j],"_cf.pdf"))
	plot(cur_dat_raw[which(cur_alt==1)] ~ cur_vaf[which(cur_alt==1)])
	dev.off()
}
}


######################
### individual cna ###
######################
for(j in 1:ncol(all_mat_cna)){
cur_alt = as.numeric(all_mat_cna[,j]!=0)
cur_cf = as.numeric(all_mat_cna[,j])
cmp_data = data.frame(cur_dat,cur_alt,all_age,all_sex,all_version,all_HM)
cmp_data = cmp_data[which(cur_alt==1 | is_none==1),]
s = summary(lm(cur_dat ~ ., data=cmp_data))
cmp_data_2 = data.frame(cur_dat,cur_cf,all_age,all_sex,all_HM)
cmp_data_2 = cmp_data_2[which(cur_cf>0),]
s2 = summary(lm(cur_dat ~ ., data=cmp_data_2))
if(rownames(s$coefficients)[2] == "cur_alt"){
	cat(colnames(all_mat_cna)[j],"\t")
	cat(lab_all[i],"\t")
	cat(sum((all_mat_cna[,j]!=0)*(!is.na(cur_dat))),"\n")

	mat_p[i,j+ncol(all_mat_mut)+7] = s$coef[2,4]; mat_t[i,j+ncol(all_mat_mut)+12] = s$coef[2,1]/sd(na.omit(cur_dat))
	mat_cf_p[i,j+ncol(all_mat_mut)] = s2$coef[2,4]; mat_cf_t[i,j+ncol(all_mat_mut)] = s2$coef[2,1]/sd(na.omit(cur_dat))

	pdf(paste0(colnames(all_dat)[i]," - ",colnames(all_mat_cna)[j],"_cf.pdf"))
	plot(cur_dat_raw[which(cur_alt==1)] ~ cur_cf[which(cur_alt==1)])
	dev.off()
}
}

}

mat_p = mat_p[,c(1:7,which(apply(all_mat_mut!=0,2,sum)>=30)+7,which(apply(all_mat_cna!=0,2,sum)>=30)+ncol(all_mat_mut)+7)]
mat_t = mat_t[,c(1:7,which(apply(all_mat_mut!=0,2,sum)>=30)+7,which(apply(all_mat_cna!=0,2,sum)>=30)+ncol(all_mat_mut)+7)]
mat_q = matrix(p.adjust(as.numeric(mat_q),"BH"),nrow=nrow(mat_q))



## plot  of correlation ##
## Draw Fig.3b ##

getCnaLabel = function(s){
splited = unlist(strsplit(s,":"))
type=NA
prf = ifelse(splited[2]=="CNN-LOH","",ifelse(splited[2]=="Duplication","+",ifelse(splited[2]=="Deletion","del(","Unknown")))
suf = ifelse(splited[2]=="CNN-LOH","UPD",ifelse(splited[2]=="Duplication","",ifelse(splited[2]=="Deletion",")","Unknown")))
return(paste0(prf,splited[1],splited[3],suf))
}

xwid = nrow(mat_p)
ywid = ncol(mat_p)
xlab = lab_all
ylab = colnames(mat_p)

fonts = c(rep(1,8),rep(3,sum(apply(all_mat_mut!=0,2,sum)>=30)),rep(1,sum(apply(all_mat_cna!=0,2,sum)>=30)))

r = 0.15; theta = c(18,18+72,18+72*2,18+72*3,18+72*4)*pi/180

pdf("Fig_3b_CH_cbc_correlation.pdf",w=ywid/5,h=ywid/4.5)
plot(NULL,NULL,
xlim=c(-ywid/10,ywid + ywid/10),
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

if(mat_p[i,j] >= 0.05){
c = "white"
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

j = j + 2
i=1
s = s.05
rect(i-.5-s,ywid-(j-.5-s),i-.5+s,ywid-(j-.5+s),col="gray",border="transparent")
i=3
s = s.01
rect(i-.5-s,ywid-(j-.5-s),i-.5+s,ywid-(j-.5+s),col="gray",border="transparent")
i=5
s = s.001
rect(i-.5-s,ywid-(j-.5-s),i-.5+s,ywid-(j-.5+s),col="gray",border="transparent")

for(i in 1:100){
t = seq(-1,1,length=100)[i]
if(t>0){c = rgb(1,0,0,alpha=min(1,t))}
if(t<0){c = rgb(0,0,1,alpha=min(1,-t))}
rect(3+(i-1)*3/100,-7.5,3+i*3/100,-7,col=c,border="transparent")
}

dev.off()



## Draw ED Fig.7c,d ##

###########################
### Hb plot n_mut n_cna ###
###########################

is_none = (1-is_mut) * (1-is_cna)

hb_none = na.omit(all_dat[which(is_none==1),"Hb"])
hb_mut = na.omit(all_dat[which(is_mut==1),"Hb"])
hb_mut_1 = na.omit(all_dat[which(all_n_alt==1),"Hb"])
hb_mut_2 = na.omit(all_dat[which(all_n_alt==2),"Hb"])
hb_mut_3 = na.omit(all_dat[which(all_n_alt>=3),"Hb"])
hb_cna = na.omit(all_dat[which(is_cna==1),"Hb"])
hb_cna_1 = na.omit(all_dat[which(all_n_cna==1),"Hb"])
hb_cna_2 = na.omit(all_dat[which(all_n_cna==2),"Hb"])
hb_cna_3 = na.omit(all_dat[which(all_n_cna>=3),"Hb"])

box = function(x){
	q = quantile(x)
	return(c(q[2]-1.5*(q[4]-q[2]),q[2],q[3],q[4],q[4]+1.5*(q[4]-q[2])))
}

box_none = box(hb_none)
box_mut = box(hb_mut)
box_mut_1 = box(hb_mut_1)
box_mut_2 = box(hb_mut_2)
box_mut_3 = box(hb_mut_3)
box_cna = box(hb_cna)
box_cna_1 = box(hb_cna_1)
box_cna_2 = box(hb_cna_2)
box_cna_3 = box(hb_cna_3)

x_none = 1*1.7 -1 + runif(length(hb_none),-0.2,0.2)
x_mut = 2*1.7 -1 + runif(length(hb_none),-0.2,0.2)
x_mut_1 = 3*1.7 -1 + runif(length(hb_none),-0.2,0.2)
x_mut_2 = 4*1.7 -1 + runif(length(hb_none),-0.2,0.2)
x_mut_3 = 5*1.7 -1 + runif(length(hb_none),-0.2,0.2)
x_cna = 6*1.7 -1 + runif(length(hb_none),-0.2,0.2)
x_cna_1 = 7*1.7 -1 + runif(length(hb_none),-0.2,0.2)
x_cna_2 = 8*1.7 -1 + runif(length(hb_none),-0.2,0.2)
x_cna_3 = 9*1.7 -1 + runif(length(hb_none),-0.2,0.2)

m = max(
hb_none,
hb_mut,hb_mut_1,hb_mut_2,hb_mut_3,
hb_cna,hb_cna_1,hb_cna_2,hb_cna_3
)

l = min(
hb_none,
hb_mut,hb_mut_1,hb_mut_2,hb_mut_3,
hb_cna,hb_cna_1,hb_cna_2,hb_cna_3
)

dat_list = list(
hb_none,
hb_mut,hb_mut_1,hb_mut_2,hb_mut_3,
hb_cna,hb_cna_1,hb_cna_2,hb_cna_3
)

x_list = list(
x_none,
x_mut,x_mut_1,x_mut_2,x_mut_3,
x_cna,x_cna_1,x_cna_2,x_cna_3
)

box_list = list(
box_none,
box_mut,box_mut_1,box_mut_2,box_mut_3,
box_cna,box_cna_1,box_cna_2,box_cna_3
)

labels = c(
"No alteration",
"SNV","1 SNV","2 SNVs","3 SNVs",
"CNA","1 CNA","2 CNAs","3 CNAs"
)

s = 7/(m-2)
print(m)
print(l)

pdf("ED_Fig_7c_Hb_plot_n_mut_cna.pdf",width=12,height=9)
plot(NULL,NULL,xlim = c(-1,15),ylim = c(2,m*(6/5)),
xlab = "",ylab = "",axes = F)

# median in normal
segments(0,box_list[[1]][3],length(box_list)*1.7+1,box_list[[1]][3],col="gray",lty="dashed",lwd=2)

segments(0,4,0,20)
text(-0.8,12,"Hemoglobin in male subjects (mg/dl)",srt=90,adj=0.5,font=2)
lats = (2:10)*2
for(i in 1:length(lats)){
	segments(0,lats[i],0-0.05,lats[i])
	text(-0.2,lats[i],lats[i],adj=1,font=2)
}
lat2 = setdiff((20:100)/5,lats)
for(i in 1:length(lat2)){
	segments(0,lat2[i],0-0.025,lat2[i])
}

theta = seq(-pi,pi,length=100)
r = 0.1
for(l in 1:length(dat_list)){

cur_dat = dat_list[[l]]
cur_x = x_list[[l]]
for(i in 1:length(cur_dat)){
	polygon(cur_x[i]+s*cos(theta)*r,cur_dat[i]+sin(theta)*r,col=rgb(0.5,0,0.5, alpha=0.3),border="transparent")
}

b = box_list[[l]]
x0 = l*1.7-0.5; w0 = 0.15; xs = c(x0-w0,x0-w0,x0,x0-w0,x0-w0,x0+w0,x0+w0,x0,x0+w0,x0+w0)
ys = c(b[2],(b[2]+3*b[3])/4,b[3],(b[4]+3*b[3])/4,b[4])
ys = c(ys,rev(ys))
polygon(xs,ys,border="black",lwd=2,col="gray")
arrows(x0,b[4],x0,max(cur_dat[which(cur_dat <= b[5])]),angle=90,lwd=2,length=w0)
arrows(x0,b[2],x0,min(cur_dat[which(cur_dat >= b[1])]),angle=90,lwd=2,length=w0)

text(l*1.7-0.75,3,labels[l],adj=0.5,font=2)
}

dev.off()



#################
### hb  n_alt ###
#################

is_none = (1-is_mut) * (1-is_cna)
all_n_alt = all_n_mut + all_n_cna

hb_none = na.omit(all_dat[which(is_none==1),"Hb"])
hb_alt_1 = na.omit(all_dat[which(all_n_alt==1),"Hb"])
hb_alt_2 = na.omit(all_dat[which(all_n_alt==2),"Hb"])
hb_alt_3 = na.omit(all_dat[which(all_n_alt==3),"Hb"])
hb_alt_4 = na.omit(all_dat[which(all_n_alt>=4),"Hb"])

wilcox.test(hb_none,hb_alt_1)
wilcox.test(hb_alt_1,hb_alt_2)
wilcox.test(hb_alt_2,hb_alt_3)
wilcox.test(hb_alt_3,hb_alt_4)

wilcox.test(hb_none,hb_alt_2)
wilcox.test(hb_alt_1,hb_alt_3)
wilcox.test(hb_alt_2,hb_alt_4)

wilcox.test(hb_none,hb_alt_3)
wilcox.test(hb_alt_1,hb_alt_4)

wilcox.test(hb_none,hb_alt_4)


hb_ = all_dat[,"Hb"]
lm_ = lm(hb_ ~ all_n_alt + all_sex + all_age + all_HM)
summary(lm_)

hb_alt_12 = all_dat[which(all_n_alt==4 | all_n_alt==2),"Hb"]
sex_alt_12 = all_sex[which(all_n_alt==4 | all_n_alt==2)]
age_alt_12 = all_age[which(all_n_alt==4 | all_n_alt==2)]
all_n_alt_12 = all_n_alt[which(all_n_alt==4 | all_n_alt==2)]
lm_12 = lm(hb_alt_12 ~ all_n_alt_12 + sex_alt_12 + age_alt_12)
summary(lm_12)

box = function(x){
	q = quantile(x)
	return(c(q[2]-1.5*(q[4]-q[2]),q[2],q[3],q[4],q[4]+1.5*(q[4]-q[2])))
}

box_none = box(hb_none)
box_alt_1 = box(hb_alt_1)
box_alt_2 = box(hb_alt_2)
box_alt_3 = box(hb_alt_3)
box_alt_4 = box(hb_alt_4)

x_none = 1*1.7 -1 + runif(length(hb_none),-0.2,0.2)
x_alt_1 = 2*1.7 -1 + runif(length(hb_none),-0.2,0.2)
x_alt_2 = 3*1.7 -1 + runif(length(hb_none),-0.2,0.2)
x_alt_3 = 4*1.7 -1 + runif(length(hb_none),-0.2,0.2)
x_alt_4 = 5*1.7 -1 + runif(length(hb_none),-0.2,0.2)

m = max(hb_none,hb_alt_1,hb_alt_2,hb_alt_3,hb_alt_4)
l = min(hb_none,hb_alt_1,hb_alt_2,hb_alt_3,hb_alt_4)
dat_list = list(hb_none,hb_alt_1,hb_alt_2,hb_alt_3,hb_alt_4)
x_list = list(x_none,x_alt_1,x_alt_2,x_alt_3,x_alt_4)
box_list = list(box_none,box_alt_1,box_alt_2,box_alt_3,box_alt_4)
labels = c("0","1","2","3",">=4")

s = 7/(m-2)
print(m)
print(l)

pdf("ED_Fig_7c_Hb_plot_n_alt.pdf",width=12,height=10)
plot(NULL,NULL,
xlim = c(-1,10),ylim = c(2,m*(6/5)),
xlab = "",ylab = "",axes = F)

segments(0,4,0,20)
lats = (2:10)*2
for(i in 1:length(lats)){
	segments(0,lats[i],0-0.05,lats[i])
	text(-0.1,lats[i],lats[i],adj=1)
}
lat2 = setdiff((20:100)/5,lats)
for(i in 1:length(lat2)){
	segments(0,lat2[i],0-0.025,lat2[i])
}

theta = seq(-pi,pi,length=100)
r = 0.1
for(l in 1:length(dat_list)){

cur_dat = dat_list[[l]]
cur_x = x_list[[l]]
for(i in 1:length(cur_dat)){
	polygon(cur_x[i]+s*cos(theta)*r,cur_dat[i]+sin(theta)*r,col=rgb(0.5,0,0.5, alpha=0.3),border="transparent")
}

b = box_list[[l]]
x0 = l*1.7-0.5; w0 = 0.15; xs = c(x0-w0,x0-w0,x0,x0-w0,x0-w0,x0+w0,x0+w0,x0,x0+w0,x0+w0)
ys = c(b[2],(b[2]+3*b[3])/4,b[3],(b[4]+3*b[3])/4,b[4])
ys = c(ys,rev(ys))
polygon(xs,ys,border="black",lwd=2,col="gray")
arrows(x0,b[4],x0,max(cur_dat[which(cur_dat <= b[5])]),angle=90,lwd=2,length=w0)
arrows(x0,b[2],x0,min(cur_dat[which(cur_dat >= b[1])]),angle=90,lwd=2,length=w0)

text(l*1.7-0.75,3,labels[l],adj=0.5,font=2)
}

# median in normal
segments(0,box_list[[1]][3],length(box_list)*1.7+1,box_list[[1]][3],col="gray",lty="dashed",lwd=2)
dev.off()



##########
### Hb ###
##########

sex = 1
is_none = (1-is_mut) * (1-is_cna)
is_none = as.numeric(is_none==1 & all_sex == sex)
is_del20q = as.numeric(all_mat_cna[,"20:Deletion:q"]!=0 & all_sex == sex)

is_ppm1d = as.numeric(all_mat_mut[,"PPM1D"]!=0 & all_sex == sex)
is_tp53 = as.numeric(all_mat_mut[,"TP53"]!=0 & all_sex == sex)
is_sf3b1 = as.numeric(all_mat_mut[,"SF3B1"]!=0 & all_sex == sex)
is_u2af1 = as.numeric(all_mat_mut[,"U2AF1"]!=0 & all_sex == sex)
is_4qupd = as.numeric(all_mat_cna[,"4:CNN-LOH:q"]!=0 & all_sex == sex)
#is_del20q = as.numeric(all_mat_cna[,"20:Deletion:q"]!=0)

hb_none = na.omit(all_dat[which(is_none==1),"Hb"])
hb_mut = na.omit(all_dat[which(is_mut==1),"Hb"])
hb_cna = na.omit(all_dat[which(is_cna==1),"Hb"])
hb_ppm1d = na.omit(all_dat[which(is_ppm1d==1),"Hb"])
hb_tp53 = na.omit(all_dat[which(is_tp53==1),"Hb"])
hb_sf3b1 = na.omit(all_dat[which(is_sf3b1==1),"Hb"])
hb_u2af1 = na.omit(all_dat[which(is_u2af1==1),"Hb"])
hb_4qupd = na.omit(all_dat[which(is_4qupd==1),"Hb"])
hb_del20q = na.omit(all_dat[which(is_del20q==1),"Hb"])

box = function(x){
	q = quantile(x)
	return(c(q[2]-1.5*(q[4]-q[2]),q[2],q[3],q[4],q[4]+1.5*(q[4]-q[2])))
}

box_none = box(hb_none)
box_mut = box(hb_mut)
box_cna = box(hb_cna)
box_ppm1d = box(hb_ppm1d)
box_tp53 = box(hb_tp53)
box_sf3b1 = box(hb_sf3b1)
box_u2af1 = box(hb_u2af1)
box_4qupd = box(hb_4qupd)
box_del20q = box(hb_del20q)

x_none = 1*1.7 -1 + runif(length(hb_none),-0.2,0.2)
x_ppm1d = 2*1.7 -1 + runif(length(hb_ppm1d),-0.2,0.2)
x_tp53 = 3*1.7 -1 + runif(length(hb_tp53),-0.2,0.2)
x_sf3b1 = 4*1.7 -1 + runif(length(hb_sf3b1),-0.2,0.2)
x_u2af1 = 5*1.7 -1 + runif(length(hb_u2af1),-0.2,0.2)
x_4qupd = 6*1.7 -1 + runif(length(hb_4qupd),-0.2,0.2)
x_del20q = 7*1.7 -1 + runif(length(hb_del20q),-0.2,0.2)

m = max(hb_none,
hb_ppm1d,
hb_tp53,
hb_sf3b1,
hb_u2af1,
hb_4qupd,
hb_del20q)

l = min(hb_none,
hb_ppm1d,
hb_tp53,
hb_sf3b1,
hb_u2af1,
hb_4qupd,
hb_del20q)

dat_list = list(
hb_none,
hb_ppm1d,
hb_tp53,
hb_sf3b1,
hb_u2af1,
hb_4qupd,
hb_del20q)

x_list = list(
x_none,
x_ppm1d,
x_tp53,
x_sf3b1,
x_u2af1,
x_4qupd,
x_del20q)

box_list = list(
box_none,
box_ppm1d,
box_tp53,
box_sf3b1,
box_u2af1,
box_4qupd,
box_del20q)

labels = c(
"No alteration",
"PPM1D",
"TP53",
"SF3B1",
"U2AF1",
"4qUPD",
"del(20q)"
)

s = 7/(m-2)
print(m)
print(l)

pdf("ED_Fig_7c_Hb_plot_male.pdf",width=12,height=10)
plot(NULL,NULL,xlim = c(-1,15),ylim = c(2,m*(6/5)),
xlab = "",ylab = "",axes = F)

segments(0,4,0,20)
lats = (2:10)*2
for(i in 1:length(lats)){
	segments(0,lats[i],0-0.05,lats[i])
	text(-0.1,lats[i],lats[i],adj=1)
}
lat2 = setdiff((20:100)/5,lats)
for(i in 1:length(lat2)){
	segments(0,lat2[i],0-0.025,lat2[i])
}

theta = seq(-pi,pi,length=100)
r = 0.1
for(l in 1:length(dat_list)){

cur_dat = dat_list[[l]]
cur_x = x_list[[l]]
for(i in 1:length(cur_dat)){
	polygon(cur_x[i]+s*cos(theta)*r,cur_dat[i]+sin(theta)*r,col=rgb(0.5,0,0.5, alpha=0.3),border="transparent")
}

b = box_list[[l]]
x0 = l*1.7-0.5; w0 = 0.15; xs = c(x0-w0,x0-w0,x0,x0-w0,x0-w0,x0+w0,x0+w0,x0,x0+w0,x0+w0)
ys = c(b[2],(b[2]+3*b[3])/4,b[3],(b[4]+3*b[3])/4,b[4])
ys = c(ys,rev(ys))
polygon(xs,ys,border="black",lwd=2,col="gray")
arrows(x0,b[4],x0,max(cur_dat[which(cur_dat <= b[5])]),angle=90,lwd=2,length=w0)
arrows(x0,b[2],x0,min(cur_dat[which(cur_dat >= b[1])]),angle=90,lwd=2,length=w0)

text(l*1.7-0.75,3,labels[l],adj=0.5,font=2)
}

# median in normal
segments(0,box_list[[1]][3],length(box_list)*1.7+1,box_list[[1]][3],col="gray",lty="dashed",lwd=2)
dev.off()



sex = 2
is_none = (1-is_mut) * (1-is_cna)
is_none = as.numeric(is_none==1 & all_sex == sex)
is_del20q = as.numeric(all_mat_cna[,"20:Deletion:q"]!=0 & all_sex == sex)

is_ppm1d = as.numeric(all_mat_mut[,"PPM1D"]!=0 & all_sex == sex)
is_tp53 = as.numeric(all_mat_mut[,"TP53"]!=0 & all_sex == sex)
is_sf3b1 = as.numeric(all_mat_mut[,"SF3B1"]!=0 & all_sex == sex)
is_u2af1 = as.numeric(all_mat_mut[,"U2AF1"]!=0 & all_sex == sex)
is_4qupd = as.numeric(all_mat_cna[,"4:CNN-LOH:q"]!=0 & all_sex == sex)
#is_del20q = as.numeric(all_mat_cna[,"20:Deletion:q"]!=0)

hb_none = na.omit(all_dat[which(is_none==1),"Hb"])
hb_mut = na.omit(all_dat[which(is_mut==1),"Hb"])
hb_mut_1 = na.omit(all_dat[which(all_n_mut==1),"Hb"])
hb_mut_2 = na.omit(all_dat[which(all_n_mut==2),"Hb"])
hb_mut_3 = na.omit(all_dat[which(all_n_mut==3),"Hb"])
hb_mut_4 = na.omit(all_dat[which(all_n_mut>=4),"Hb"])
hb_cna = na.omit(all_dat[which(is_cna==1),"Hb"])
hb_cna_1 = na.omit(all_dat[which(all_n_cna==1),"Hb"])
hb_cna_2 = na.omit(all_dat[which(all_n_cna==2),"Hb"])
hb_cna_3 = na.omit(all_dat[which(all_n_cna==3),"Hb"])
hb_cna_4 = na.omit(all_dat[which(all_n_cna>=4),"Hb"])
hb_ppm1d = na.omit(all_dat[which(is_ppm1d==1),"Hb"])
hb_tp53 = na.omit(all_dat[which(is_tp53==1),"Hb"])
hb_sf3b1 = na.omit(all_dat[which(is_sf3b1==1),"Hb"])
hb_u2af1 = na.omit(all_dat[which(is_u2af1==1),"Hb"])
hb_4qupd = na.omit(all_dat[which(is_4qupd==1),"Hb"])
hb_del20q = na.omit(all_dat[which(is_del20q==1),"Hb"])

box = function(x){
	q = quantile(x)
	return(c(q[2]-1.5*(q[4]-q[2]),q[2],q[3],q[4],q[4]+1.5*(q[4]-q[2])))
}

box_none = box(hb_none)
box_mut = box(hb_mut)
box_mut_1 = box(hb_mut_1)
box_mut_2 = box(hb_mut_2)
box_mut_3 = box(hb_mut_3)
box_mut_4 = box(hb_mut_4)
box_cna = box(hb_cna)
box_cna_1 = box(hb_cna_1)
box_cna_2 = box(hb_cna_2)
box_cna_3 = box(hb_cna_3)
box_cna_4 = box(hb_cna_4)
box_ppm1d = box(hb_ppm1d)
box_tp53 = box(hb_tp53)
box_sf3b1 = box(hb_sf3b1)
box_u2af1 = box(hb_u2af1)
box_4qupd = box(hb_4qupd)
box_del20q = box(hb_del20q)

x_none = 1*1.7 -1 + runif(length(hb_none),-0.2,0.2)
x_ppm1d = 2*1.7 -1 + runif(length(hb_ppm1d),-0.2,0.2)
x_tp53 = 3*1.7 -1 + runif(length(hb_tp53),-0.2,0.2)
x_sf3b1 = 4*1.7 -1 + runif(length(hb_sf3b1),-0.2,0.2)
x_u2af1 = 5*1.7 -1 + runif(length(hb_u2af1),-0.2,0.2)
x_4qupd = 6*1.7 -1 + runif(length(hb_4qupd),-0.2,0.2)
x_del20q = 7*1.7 -1 + runif(length(hb_del20q),-0.2,0.2)

m = max(hb_none,
hb_ppm1d,
hb_tp53,
hb_sf3b1,
hb_u2af1,
hb_4qupd,
hb_del20q)

l = min(hb_none,
hb_ppm1d,
hb_tp53,
hb_sf3b1,
hb_u2af1,
hb_4qupd,
hb_del20q)

dat_list = list(
hb_none,
hb_ppm1d,
hb_tp53,
hb_sf3b1,
hb_u2af1,
hb_4qupd,
hb_del20q)

x_list = list(
x_none,
x_ppm1d,
x_tp53,
x_sf3b1,
x_u2af1,
x_4qupd,
x_del20q)

box_list = list(
box_none,
box_ppm1d,
box_tp53,
box_sf3b1,
box_u2af1,
box_4qupd,
box_del20q)

labels = c(
"No alteration",
"PPM1D",
"TP53",
"SF3B1",
"U2AF1",
"4qUPD",
"del(20q)"
)

s = 7/(m-2)
print(m)
print(l)

pdf("ED_Fig_7c_Hb_plot_female.pdf",width=12,height=10)
plot(NULL,NULL,xlim = c(-1,15),ylim = c(2,m*(6/5)),
xlab = "",ylab = "",axes = F)

segments(0,4,0,20)
lats = (2:10)*2
for(i in 1:length(lats)){
	segments(0,lats[i],0-0.05,lats[i])
	text(-0.1,lats[i],lats[i],adj=1)
}
lat2 = setdiff((20:100)/5,lats)
for(i in 1:length(lat2)){
	segments(0,lat2[i],0-0.025,lat2[i])
}

theta = seq(-pi,pi,length=100)
r = 0.1
for(l in 1:length(dat_list)){

cur_dat = dat_list[[l]]
cur_x = x_list[[l]]
for(i in 1:length(cur_dat)){
	polygon(cur_x[i]+s*cos(theta)*r,cur_dat[i]+sin(theta)*r,col=rgb(0.5,0,0.5, alpha=0.3),border="transparent")
}

b = box_list[[l]]
x0 = l*1.7-0.5; w0 = 0.15; xs = c(x0-w0,x0-w0,x0,x0-w0,x0-w0,x0+w0,x0+w0,x0,x0+w0,x0+w0)
ys = c(b[2],(b[2]+3*b[3])/4,b[3],(b[4]+3*b[3])/4,b[4])
ys = c(ys,rev(ys))
polygon(xs,ys,border="black",lwd=2,col="gray")
arrows(x0,b[4],x0,max(cur_dat[which(cur_dat <= b[5])]),angle=90,lwd=2,length=w0)
arrows(x0,b[2],x0,min(cur_dat[which(cur_dat >= b[1])]),angle=90,lwd=2,length=w0)

text(l*1.7-0.75,3,labels[l],adj=0.5,font=2)
}

# median in normal
segments(0,box_list[[1]][3],length(box_list)*1.7+1,box_list[[1]][3],col="gray",lty="dashed",lwd=2)
dev.off()



###########
### MCV ###
###########

is_none = (1-is_mut) * (1-is_cna)
is_tp53 = as.numeric(all_mat_mut[,"TP53"]!=0)
is_sf3b1 = as.numeric(all_mat_mut[,"SF3B1"]!=0)
is_11qupd = as.numeric(all_mat_cna[,"11:CNN-LOH:q"]!=0)

mcv_none = na.omit(all_dat[which(is_none==1),"MCV"])
mcv_mut = na.omit(all_dat[which(is_mut==1),"MCV"])
mcv_cna = na.omit(all_dat[which(is_cna==1),"MCV"])
mcv_tp53 = na.omit(all_dat[which(is_tp53==1),"MCV"])
mcv_sf3b1 = na.omit(all_dat[which(is_sf3b1==1),"MCV"])
mcv_11qupd = na.omit(all_dat[which(is_11qupd==1),"MCV"])

box = function(x){
	q = quantile(x)
	return(c(q[2]-1.5*(q[4]-q[2]),q[2],q[3],q[4],q[4]+1.5*(q[4]-q[2])))
}

box_none = box(mcv_none)
box_mut = box(mcv_mut)
box_cna = box(mcv_cna)
box_tp53 = box(mcv_tp53)
box_sf3b1 = box(mcv_sf3b1)
box_11qupd = box(mcv_11qupd)

x_none = 1*1.7 -1 + runif(length(mcv_none),-0.2,0.2)
x_mut = 2*1.7 -1 + runif(length(mcv_mut),-0.2,0.2)
x_cna = 3*1.7 -1 + runif(length(mcv_cna),-0.2,0.2)
x_tp53 = 4*1.7 -1 + runif(length(mcv_tp53),-0.2,0.2)
x_sf3b1 = 5*1.7 -1 + runif(length(mcv_sf3b1),-0.2,0.2)
x_11qupd = 6*1.7 -1 + runif(length(mcv_11qupd),-0.2,0.2)

m = max(mcv_none,
mcv_mut,
mcv_cna,
mcv_tp53,
mcv_sf3b1,
mcv_11qupd
)

l = min(mcv_none,
mcv_mut,
mcv_cna,
mcv_tp53,
mcv_sf3b1,
mcv_11qupd
)

dat_list = list(
mcv_none,
mcv_mut,
mcv_cna,
mcv_tp53,
mcv_sf3b1,
mcv_11qupd
)

x_list = list(
x_none,
x_mut,
x_cna,
x_tp53,
x_sf3b1,
x_11qupd
)

box_list = list(
box_none,
box_mut,
box_cna,
box_tp53,
box_sf3b1,
box_11qupd
)

labels = c(
"No alteration",
"SNV/indels",
"CNAs",
"TP53",
"SF3B1",
"11qUPD"
)

s = 7/(m-2)
print(m)
print(l)

pdf("ED_Fig_7c_MCV_plot.pdf",width=13,height=8)
plot(NULL,NULL,xlim = c(-1,10),ylim = c(50,150),
xlab = "",ylab = "",axes = F)

# median in normal
segments(0,box_list[[1]][3],length(box_list)*1.7+1,box_list[[1]][3],col="gray",lty="dashed",lwd=3)

segments(0,60,0,140)
lats = (6:14)*10
for(i in 1:length(lats)){
	segments(0,lats[i],0-0.05,lats[i])
	text(-0.1,lats[i],lats[i],adj=1)
}
lat2 = setdiff((12:28)*5,lats)
for(i in 1:length(lat2)){
	segments(0,lat2[i],0-0.025,lat2[i])
}

theta = seq(-pi,pi,length=100)
r = 1
for(l in 1:length(dat_list)){

cur_dat = dat_list[[l]]
cur_x = x_list[[l]]
for(i in 1:length(cur_dat)){
	polygon(cur_x[i]+s*cos(theta)*r,cur_dat[i]+sin(theta)*r,col=rgb(0.5,0,0.5, alpha=0.3),border="transparent")
}

b = box_list[[l]]
x0 = l*1.7-0.5; w0 = 0.15; xs = c(x0-w0,x0-w0,x0,x0-w0,x0-w0,x0+w0,x0+w0,x0,x0+w0,x0+w0)
ys = c(b[2],(b[2]+3*b[3])/4,b[3],(b[4]+3*b[3])/4,b[4])
ys = c(ys,rev(ys))
polygon(xs,ys,border="black",lwd=2,col="gray")
arrows(x0,b[4],x0,max(cur_dat[which(cur_dat <= b[5])]),angle=90,lwd=2,length=w0)
arrows(x0,b[2],x0,min(cur_dat[which(cur_dat >= b[1])]),angle=90,lwd=2,length=w0)

text(l*1.7-0.75,60,labels[l],adj=0.5,font=2)
}

dev.off()



###########
### Plt ###
###########

is_none = (1-is_mut) * (1-is_cna)
is_ppm1d = as.numeric(all_mat_mut[,"PPM1D"]!=0)
is_u1af1 = as.numeric(all_mat_mut[,"U2AF1"]!=0)
is_jak2= as.numeric(all_mat_mut[,"JAK2"]!=0)
is_6pupd = as.numeric(all_mat_cna[,"6:CNN-LOH:p"]!=0)
is_9pupd = as.numeric(all_mat_cna[,"9:CNN-LOH:p"]!=0)
is_11qupd = as.numeric(all_mat_cna[,"11:CNN-LOH:q"]!=0)
is_del20q = as.numeric(all_mat_cna[,"20:Deletion:q"]!=0)

plt_none = na.omit(all_dat[which(is_none==1),"PLT"])
plt_mut = na.omit(all_dat[which(is_mut==1),"PLT"])
plt_cna = na.omit(all_dat[which(is_cna==1),"PLT"])
plt_jak2 = na.omit(all_dat[which(is_jak2==1),"PLT"])
plt_9pupd = na.omit(all_dat[which(is_9pupd==1),"PLT"])
plt_ppm1d = na.omit(all_dat[which(is_ppm1d==1),"PLT"])
plt_u2af1 = na.omit(all_dat[which(is_u2af1==1),"PLT"])
plt_6pupd = na.omit(all_dat[which(is_6pupd==1),"PLT"])
plt_11qupd = na.omit(all_dat[which(is_11qupd==1),"PLT"])
plt_del20q = na.omit(all_dat[which(is_del20q==1),"PLT"])

box = function(x){
	q = quantile(x)
	return(c(q[2]-1.5*(q[4]-q[2]),q[2],q[3],q[4],q[4]+1.5*(q[4]-q[2])))
}

box_none = box(plt_none)
box_mut = box(plt_mut)
box_cna = box(plt_cna)
box_jak2 = box(plt_jak2)
box_9pupd = box(plt_9pupd)
box_ppm1d = box(plt_ppm1d)
box_u2af1 = box(plt_u2af1)
box_6pupd = box(plt_6pupd)
box_11qupd = box(plt_11qupd)
box_del20q = box(plt_del20q)

x_none = 1*1.7 -1 + runif(length(plt_none),-0.2,0.2)
x_mut = 2*1.7 -1 + runif(length(plt_mut),-0.2,0.2)
x_cna = 3*1.7 -1 + runif(length(plt_cna),-0.2,0.2)
x_jak2 = 4*1.7 -1 + runif(length(plt_jak2),-0.2,0.2)
x_9pupd = 5*1.7 -1 + runif(length(plt_9pupd),-0.2,0.2)
x_ppm1d = 6*1.7 -1 + runif(length(plt_ppm1d),-0.2,0.2)
x_u2af1 = 7*1.7 -1 + runif(length(plt_u2af1),-0.2,0.2)
x_6pupd = 8*1.7 -1 + runif(length(plt_6pupd),-0.2,0.2)
x_11qupd = 9*1.7 -1 + runif(length(plt_11qupd),-0.2,0.2)
x_del20q = 10*1.7 -1 + runif(length(plt_del20q),-0.2,0.2)

m = max(mcv_none,
plt_none,
plt_mut,
plt_cna,
plt_jak2,
plt_9pupd,
plt_ppm1d,
plt_u2af1,
plt_6pupd,
plt_11qupd,
plt_del20q
)

l = min(mcv_none,
plt_none,
plt_mut,
plt_cna,
plt_jak2,
plt_9pupd,
plt_ppm1d,
plt_u2af1,
plt_6pupd,
plt_11qupd,
plt_del20q
)

dat_list = list(
plt_none,
plt_mut,
plt_cna,
plt_jak2,
plt_9pupd,
plt_ppm1d,
plt_u2af1,
plt_6pupd,
plt_11qupd,
plt_del20q
)

x_list = list(
x_none,
x_mut,
x_cna,
x_jak2,
x_9pupd,
x_ppm1d,
x_u2af1,
x_6pupd,
x_11qupd,
x_del20q
)

box_list = list(
box_none,
box_mut,
box_cna,
box_jak2,
box_9pupd,
box_ppm1d,
box_u2af1,
box_6pupd,
box_11qupd,
box_del20q
)

labels = c(
"No alteration",
"SNV/indels",
"CNAs",
"JAK2",
"9pUPD",
"PPM1D",
"U2AF1",
"6pUPD",
"11qUPD",
"del(20q)"
)

s = 7/(m-2)
print(m)
print(l)

pdf("ED_Fig_7c_Plt_plot.pdf",width=35,height=10)

plot(
NULL,NULL,
xlim = c(-1,20),
ylim = c(-10,80),
xlab = "",
ylab = "",
axes = F
)

# median in normal
segments(0,box_list[[1]][3],length(box_list)*1.7+1,box_list[[1]][3],col="gray",lty="dashed",lwd=3)

segments(0,0,0,70)
lats = (0:7)*10
for(i in 1:length(lats)){
	segments(0,lats[i],0-0.05,lats[i])
	text(-0.1,lats[i],lats[i],adj=1)
}
lat2 = setdiff((0:12)*5,lats)
for(i in 1:length(lat2)){
	segments(0,lat2[i],0-0.025,lat2[i])
}

theta = seq(-pi,pi,length=100)
r = 1
for(l in 1:length(dat_list)){

cur_dat = dat_list[[l]]
cur_x = x_list[[l]]
for(i in 1:length(cur_dat)){
	polygon(cur_x[i]+s*cos(theta)*r,cur_dat[i]+sin(theta)*r,col=rgb(0.5,0,0.5, alpha=0.3),border="transparent")
}

b = box_list[[l]]
x0 = l*1.7-0.5; w0 = 0.15; xs = c(x0-w0,x0-w0,x0,x0-w0,x0-w0,x0+w0,x0+w0,x0,x0+w0,x0+w0)
ys = c(b[2],(b[2]+3*b[3])/4,b[3],(b[4]+3*b[3])/4,b[4])
ys = c(ys,rev(ys))
polygon(xs,ys,border="black",lwd=2,col="gray")
arrows(x0,b[4],x0,max(cur_dat[which(cur_dat <= b[5])]),angle=90,lwd=2,length=w0)
arrows(x0,b[2],x0,min(cur_dat[which(cur_dat >= b[1])]),angle=90,lwd=2,length=w0)

text(l*1.7-0.75,-5,labels[l],adj=0.5,font=2)
}

dev.off()


###########
### WBC ###
###########

is_none = (1-is_mut) * (1-is_cna)

is_14qupd = as.numeric(all_mat_cna[,"14:CNN-LOH:q"]!=0 )
is_del20q = as.numeric(all_mat_cna[,"20:Deletion:q"]!=0)

WBC_none = na.omit(all_dat[which(is_none==1),"WBC"])/1000
WBC_cna = na.omit(all_dat[which(is_cna==1),"WBC"])/1000
WBC_14qupd = na.omit(all_dat[which(is_14qupd==1),"WBC"])/1000
WBC_del20q = na.omit(all_dat[which(is_del20q==1),"WBC"])/1000

box = function(x){
	q = quantile(x)
	return(c(q[2]-1.5*(q[4]-q[2]),q[2],q[3],q[4],q[4]+1.5*(q[4]-q[2])))
}

box_none = box(WBC_none)
box_cna = box(WBC_cna)
box_14qupd = box(WBC_14qupd)
box_del20q = box(WBC_del20q)

x_none = 1*1.7 -1 + runif(length(WBC_none),-0.2,0.2)
x_cna = 2*1.7 -1 + runif(length(WBC_cna),-0.2,0.2)
x_14qupd = 3*1.7 -1 + runif(length(WBC_14qupd),-0.2,0.2)
x_del20q = 4*1.7 -1 + runif(length(WBC_del20q),-0.2,0.2)

m = max(WBC_none,
WBC_cna,
WBC_14qupd,
WBC_del20q)

l = min(WBC_none,
WBC_cna,
WBC_14qupd,
WBC_del20q)

dat_list = list(
WBC_none,
WBC_cna,
WBC_14qupd,
WBC_del20q)

dat_list

x_list = list(
x_none,
x_cna,
x_14qupd,
x_del20q)

box_list = list(
box_none,
box_cna,
box_14qupd,
box_del20q)

labels = c(
"No alteration",
"CNA",
"14qUPD",
"del(20q)"
)

s = 7/(m-2)
print(m)
print(l)

pdf("ED_Fig_7c_WBC_plot.pdf",width=12,height=10)

plot(
NULL,NULL,
xlim = c(-1,15),
ylim = c(-10,25),
xlab = "",
ylab = "",
axes = F
)

segments(0,0,0,20)
lats = (0:10)*2
for(i in 1:length(lats)){
	segments(0,lats[i],0-0.05,lats[i])
	text(-0.1,lats[i],lats[i]*1000,adj=1)
}
lat2 = setdiff(0:20,lats)
for(i in 1:length(lat2)){
	segments(0,lat2[i],0-0.025,lat2[i])
}

theta = seq(-pi,pi,length=100)
r = 0.15
for(l in 1:length(dat_list)){

cur_dat = dat_list[[l]]
cur_x = x_list[[l]]
for(i in 1:length(cur_dat)){
	polygon(cur_x[i]+s*cos(theta)*r,cur_dat[i]+sin(theta)*r,col=rgb(0.5,0,0.5, alpha=0.3),border="transparent")
}

b = box_list[[l]]
x0 = l*1.7-0.5; w0 = 0.15; xs = c(x0-w0,x0-w0,x0,x0-w0,x0-w0,x0+w0,x0+w0,x0,x0+w0,x0+w0)
ys = c(b[2],(b[2]+3*b[3])/4,b[3],(b[4]+3*b[3])/4,b[4])
ys = c(ys,rev(ys))
polygon(xs,ys,border="black",lwd=2,col="gray")
arrows(x0,b[4],x0,max(cur_dat[which(cur_dat <= b[5])]),angle=90,lwd=2,length=w0)
arrows(x0,b[2],x0,min(cur_dat[which(cur_dat >= b[1])]),angle=90,lwd=2,length=w0)

text(l*1.7-0.75,3,labels[l],adj=0.5,font=2)
}

# median in normal
segments(0,box_list[[1]][3],length(box_list)*1.7+1,box_list[[1]][3],col="gray",lty="dashed",lwd=2)
dev.off()



########################
## clone size and CBC ##
########################

pdf("ED_Fig_7d_cf_cbc_1.pdf",width=9,height=5)
par(mfrow=c(1,2))

vaf_mut = apply(all_mat_mut,1,function(x){max(as.numeric(x))})
hb_mut = all_dat[which(is_mut ==1),"Hb"]
vaf_mut = vaf_mut[which(is_mut==1)]
age_mut = all_age[which(is_mut==1)]
sex_mut = all_sex[which(is_mut==1)]

model1 = lm(hb_mut ~ vaf_mut + age_mut + sex_mut)
model = lm(hb_mut ~ vaf_mut)
coef1 = summary(model)$coef[1,1]
coef2 = summary(model)$coef[2,1]
plot(hb_mut ~ vaf_mut,col=rgb(0.5,0,0.5, alpha=0.6),pch=20,ylim=c(5,20))
abline(a=coef1,b=coef2)
abline(a=coef1,b=coef2,col="gray",lwd=2)


vaf_cna = apply(all_mat_cna,1,function(x){max(as.numeric(x))})
mcv_cna = all_dat[which(is_cna ==1),"MCV"]
vaf_cna = vaf_cna[which(is_cna==1)]
age_cna = all_age[which(is_cna==1)]
sex_cna = all_sex[which(is_cna==1)]

model1 = lm(mcv_cna ~ vaf_cna + age_cna + sex_cna)
model = lm(mcv_cna ~ vaf_cna)
coef1 = summary(model)$coef[1,1]
coef2 = summary(model)$coef[2,1]
plot(mcv_cna ~ vaf_cna,col=rgb(0.5,0,0.5, alpha=0.6),pch=20,ylim=c(60,140))
abline(a=coef1,b=coef2)
abline(a=coef1,b=coef2,col="gray",lwd=2)

dev.off()



pdf("ED_Fig_7d_cf_cbc_2.pdf",width=9,height=6.5)
par(mfrow=c(2,3))

vaf_ppm1d = as.numeric(all_mat_mut[,"PPM1D"])
is_ppm1d = as.numeric(vaf_ppm1d!=0)
hb_ppm1d = all_dat[which(is_ppm1d==1),"Hb"]
vaf_ppm1d = vaf_ppm1d[which(is_ppm1d==1)]
age_ppm1d = all_age[which(is_ppm1d==1)]
sex_ppm1d = all_sex[which(is_ppm1d==1)]

model1 = lm(hb_ppm1d ~ vaf_ppm1d + age_ppm1d + sex_ppm1d)
model = lm(hb_ppm1d ~ vaf_ppm1d)
coef1 = summary(model)$coef[1,1]
coef2 = summary(model)$coef[2,1]
plot(hb_ppm1d ~ vaf_ppm1d,col=rgb(0.5,0,0.5, alpha=0.6),ylim=c(8,18),pch=20)
abline(a=coef1,b=coef2)
abline(a=coef1,b=coef2,col="gray",lwd=2)


vaf_sf3b1 = as.numeric(all_mat_mut[,"SF3B1"])
is_sf3b1 = as.numeric(vaf_sf3b1!=0)
hb_sf3b1 = all_dat[which(is_sf3b1==1),"Hb"]
vaf_sf3b1 = vaf_sf3b1[which(is_sf3b1==1)]
age_sf3b1 = all_age[which(is_sf3b1==1)]
sex_sf3b1 = all_sex[which(is_sf3b1==1)]

model1 = lm(hb_sf3b1 ~ vaf_sf3b1 + age_sf3b1 + sex_sf3b1)
model = lm(hb_sf3b1 ~ vaf_sf3b1)
coef1 = summary(model)$coef[1,1]
coef2 = summary(model)$coef[2,1]
plot(hb_sf3b1 ~ vaf_sf3b1,col=rgb(0.5,0,0.5, alpha=0.6),ylim=c(8,18),pch=20)
abline(a=coef1,b=coef2)
abline(a=coef1,b=coef2,col="gray",lwd=2)


vaf_sf3b1 = as.numeric(all_mat_mut[,"SF3B1"])
is_sf3b1 = as.numeric(vaf_sf3b1!=0)
mcv_sf3b1 = all_dat[which(is_sf3b1==1),"MCV"]
vaf_sf3b1 = vaf_sf3b1[which(is_sf3b1==1)]
age_sf3b1 = all_age[which(is_sf3b1==1)]
sex_sf3b1 = all_sex[which(is_sf3b1==1)]

model1 = lm(mcv_sf3b1 ~ vaf_sf3b1 + age_sf3b1 + sex_sf3b1)
model = lm(mcv_sf3b1 ~ vaf_sf3b1)
coef1 = summary(model)$coef[1,1]
coef2 = summary(model)$coef[2,1]
plot(mcv_sf3b1 ~ vaf_sf3b1,col=rgb(0.5,0,0.5, alpha=0.6),pch=20)
abline(a=coef1,b=coef2)
abline(a=coef1,b=coef2,col="gray",lwd=2)



vaf_u2af1 = as.numeric(all_mat_mut[,"U2AF1"])
is_u2af1 = as.numeric(vaf_u2af1!=0)
plt_u2af1 = all_dat[which(is_u2af1==1),"PLT"]
vaf_u2af1 = vaf_u2af1[which(is_u2af1==1)]
age_u2af1 = all_age[which(is_u2af1==1)]
sex_u2af1 = all_sex[which(is_u2af1==1)]

model1 = lm(plt_u2af1 ~ vaf_u2af1 + age_u2af1 + sex_u2af1)
model = lm(plt_u2af1 ~ vaf_u2af1)
coef1 = summary(model)$coef[1,1]
coef2 = summary(model)$coef[2,1]
plot(plt_u2af1 ~ vaf_u2af1,col=rgb(0.5,0,0.5, alpha=0.6),pch=20,ylim=c(0,50))
abline(a=coef1,b=coef2,col="gray",lwd=2)



vaf_6pupd = as.numeric(all_mat_cna[,"6:CNN-LOH:p"])
is_6pupd = as.numeric(vaf_6pupd !=0)
plt_6pupd = all_dat[which(is_6pupd ==1),"PLT"]
vaf_6pupd = vaf_6pupd[which(is_6pupd==1)]
age_6pupd = all_age[which(is_6pupd==1)]
sex_6pupd = all_sex[which(is_6pupd==1)]

model1 = lm(plt_6pupd ~ vaf_6pupd + age_6pupd + sex_6pupd)
model = lm(plt_6pupd ~ vaf_6pupd)
coef1 = summary(model)$coef[1,1]
coef2 = summary(model)$coef[2,1]
plot(plt_6pupd ~ vaf_6pupd,col=rgb(0.5,0,0.5, alpha=0.6),pch=20,ylim=c(0,50))
abline(a=coef1,b=coef2)
abline(a=coef1,b=coef2,col="gray",lwd=2)


dev.off()


