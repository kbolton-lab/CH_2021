
## import data
mut_path = "MatMed_2021_mut_table.txt"
cna_path = "MatMed_2021_cna_table.csv"  
id = all_id = as.character(1:11234)


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
is_all = as.numeric((is_mut + is_cna) != 0)
is_none = (1-is_mut) * (1-is_cna)

vaf = as.numeric(apply(all_mat_mut,1,function(x){max(as.numeric(x))}))
cf = as.numeric(apply(all_mat_cna,1,function(x){max(as.numeric(x))}))
all_n_mut = n_mut
all_n_cna = n_cna

all_mat_mut_2 = all_mat_mut
all_mat_cna_2 = all_mat_cna

all_mat_mut = all_mat_mut[,which(apply(all_mat_mut != 0,2,sum) >= 20)]
all_mat_cna = all_mat_cna[,which(apply(all_mat_cna != 0,2,sum) >= 15)]



## analysis ##

dat = data.frame()

rownames(dat) = all_id
dat$WBC_ul_WBC_selected = ~~ ## WBC data
dat$Hb_g_per_dl_Hb_selected = ~~ ## Hb data
dat$Ht_percent_Hb_selected = ~~ ## Ht data
dat$PLT_10_4_per_ul_Plt_selected = ~~ ## PLT data

WBC_l = as.numeric(dat$WBC_ul_WBC_selected <= 3000)
WBC_h = as.numeric(dat$WBC_ul_WBC_selected > 10000)
names(WBC_l) = names(WBC_h) = all_id

Hb_l = as.numeric(dat$Hb_g_per_dl_Hb_selected <= 10)
Hb_h = as.numeric((dat$Hb_g_per_dl_Hb_selected >= 16.5 & all_sex == 1) | (dat$Hb_g_per_dl_Hb_selected >= 16 & all_sex == 1))
names(Hb_l) = names(Hb_h) = all_id

Ht_h = as.numeric(dat$Ht_percent_Hb_selected > 49)
names(Ht_h) = all_id

PLT_l = as.numeric(dat$PLT_10_4_per_ul_Plt_selected <= 10)
PLT_h = as.numeric(dat$PLT_10_4_per_ul_Plt_selected > 45)
names(PLT_l) = names(PLT_h) = all_id



print("All available")
sum(!is.na(dat$WBC_ul_WBC_selected) & !is.na(dat$Hb_g_per_dl_Hb_selected) & !is.na(dat$Ht_percent_Hb_selected)  & !is.na(dat$PLT_10_4_per_ul_Plt_selected))

id_1 = all_id[which(!is.na(dat$WBC_ul_WBC_selected) & !is.na(dat$Hb_g_per_dl_Hb_selected) & !is.na(dat$Ht_percent_Hb_selected)  & !is.na(dat$PLT_10_4_per_ul_Plt_selected))]
id_2 = setdiff(all_id,id_1)

print("Normal")
sum(na.omit(WBC_h+WBC_l+Hb_h+Hb_l+PLT_h+PLT_l+Ht_h)==0)
id_3 = all_id[which((WBC_h+WBC_l+Hb_h+Hb_l+PLT_h+PLT_l+Ht_h)==0)]

print("Abnormal")
sum(na.omit(WBC_h+WBC_l+Hb_h+Hb_l+PLT_h+PLT_l+Ht_h)!=0)
id_4 = all_id[which((WBC_h+WBC_l+Hb_h+Hb_l+PLT_h+PLT_l+Ht_h)!=0)]


print("Some available")
sum(!is.na(dat$WBC_ul_WBC_selected) | !is.na(dat$Hb_g_per_dl_Hb_selected) | !is.na(dat$Ht_percent_Hb_selected) | !is.na(dat$MCV_fl_Hb_selected) | !is.na(dat$MCH_pg_Hb_selected) | !is.na(dat$MCHC_percent_Hb_selected) | !is.na(dat$PLT_10_4_per_ul_Plt_selected))

sum(na.omit(WBC_l==1 | WBC_h==1 | Hb_l==1 | Hb_h==1 | Ht_h==1 | PLT_l==1 | PLT_h==1))



cytopenia_al_1 = as.numeric(WBC_l + Hb_l + PLT_l >= 1)
cytopenia_al_2 = as.numeric(WBC_l + Hb_l + PLT_l >= 2)
any_abnormality = as.numeric(WBC_l + Hb_l + PLT_l + WBC_h + Hb_h + PLT_h >= 1)

cbc = cbind(WBC_l,WBC_h,Hb_l,Hb_h,PLT_l,PLT_h,cytopenia_al_1,cytopenia_al_2,any_abnormality)


for(i in 1:ncol(cbc)){
	c = na.omit(cbc[,i])
	print(colnames(cbc)[i])
	print(sum(c))
	print(sum(c)*100/length(c))
}


genotype = cbind(all_mat_mut!=0,all_mat_cna!=0)

mut = (all_mat_mut!=0)
is_mut = as.numeric(apply(mut,1,sum) >= 1)
mut = mut[,which(apply(mut,2,sum) >= 35)]  
mut = mut[,order(apply(mut,2,sum),decreasing=T)]

cna = (all_mat_cna!=0)
is_cna = as.numeric(apply(cna,1,sum) >= 1)
cna = cna[,which(apply(cna,2,sum) >= 35)]
cna = cna[,order(apply(cna,2,sum),decreasing=T)]

print("Hb_normal")
sum(na.omit(is_mut*as.numeric(Hb_l+Hb_h==0)))
sum(na.omit(is_cna*as.numeric(Hb_l+Hb_h==0)))
sum(na.omit(is_mut*is_cna*as.numeric(Hb_l+Hb_h==0)))
sum(na.omit(as.numeric(is_mut+is_cna!=0)*as.numeric(Hb_l+Hb_h==0)))

print("Hb_low")
sum(na.omit(is_mut*Hb_l))
sum(na.omit(is_cna*Hb_l))
sum(na.omit(is_mut*is_cna*Hb_l))
sum(na.omit(as.numeric(is_mut+is_cna!=0)*Hb_l))

print("Hb_high")
sum(na.omit(is_mut*Hb_h))
sum(na.omit(is_cna*Hb_h))
sum(na.omit(is_mut*is_cna*Hb_h))
sum(na.omit(as.numeric(is_mut+is_cna!=0)*Hb_h))

print("Plt_high")
sum(na.omit(is_mut*PLT_h))
sum(na.omit(is_cna*PLT_h))
sum(na.omit(is_mut*is_cna*PLT_h))
sum(na.omit(as.numeric(is_mut+is_cna!=0)*PLT_h))

print("Plt_low")
sum(na.omit(is_mut*PLT_l))
sum(na.omit(is_cna*PLT_l))
sum(na.omit(is_mut*is_cna*PLT_l))
sum(na.omit(as.numeric(is_mut+is_cna!=0)*PLT_l))

print("Plt_normal")
sum(na.omit(is_mut*as.numeric(PLT_l+PLT_h==0)))
sum(na.omit(is_cna*as.numeric(PLT_l+PLT_h==0)))
sum(na.omit(is_mut*is_cna*as.numeric(PLT_l+PLT_h==0)))
sum(na.omit(as.numeric(is_mut+is_cna!=0)*as.numeric(PLT_l+PLT_h==0)))

print("WBC_high")
sum(na.omit(is_mut*WBC_h))
sum(na.omit(is_cna*WBC_h))
sum(na.omit(is_mut*is_cna*WBC_h))
sum(na.omit(as.numeric(is_mut+is_cna!=0)*WBC_h))

print("WBC_low")
sum(na.omit(is_mut*WBC_l))
sum(na.omit(is_cna*WBC_l))
sum(na.omit(is_mut*is_cna*WBC_l))
sum(na.omit(as.numeric(is_mut+is_cna!=0)*WBC_l))

print("WBC_normal")
sum(na.omit(is_mut*as.numeric(WBC_l+WBC_h==0)))
sum(na.omit(is_cna*as.numeric(WBC_l+WBC_h==0)))
sum(na.omit(is_mut*is_cna*as.numeric(WBC_l+WBC_h==0)))
sum(na.omit(as.numeric(is_mut+is_cna!=0)*as.numeric(WBC_l+WBC_h==0)))

print("Ht_normal")
sum(na.omit(is_mut*as.numeric(Ht_h==0)))
sum(na.omit(is_cna*as.numeric(Ht_h==0)))
sum(na.omit(is_mut*is_cna*as.numeric(Ht_h==0)))
sum(na.omit(as.numeric(is_mut+is_cna!=0)*as.numeric(Ht_h==0)))

print("Ht_high")
sum(na.omit(is_mut*as.numeric(Ht_h==1)))
sum(na.omit(is_cna*as.numeric(Ht_h==1)))
sum(na.omit(is_mut*is_cna*as.numeric(Ht_h==1)))
sum(na.omit(as.numeric(is_mut+is_cna!=0)*as.numeric(Ht_h==1)))


#genotype = cbind(is_mut,mut,is_cna,cna)
genotype = cbind(mut,cna)
#genotype = genotype[,which(apply(na.omit(genotype),2,sum)/10852 > 0.003)]
head(genotype)
# genotype = genotype[,-grep("Unknown",colnames(genotype))]

e = pv = matrix(nrow=ncol(genotype),ncol=ncol(cbc))
rownames(e) = rownames(pv) = colnames(genotype)
colnames(e) = colnames(pv) = colnames(cbc)

for(i in 1:ncol(genotype)){
	g = genotype[,i]
	for (j in 1:ncol(cbc)){
		c = cbc[,j]
		cur_data = na.omit(cbind(g,c))
		g1 = as.numeric(cur_data[,1])
		c1 = cur_data[,2]
		pp = sum(g1*c1)
		pn = sum(g1) - pp
		np = sum(c1) - pp
		nn = nrow(cur_data) - pp - pn - np
		pp = pp + 0.1
		pn = pn + 0.1
		np = np + 0.1
		nn = nn + 0.1
		res = fisher.test(matrix(c(pp,pn,np,nn),ncol=2))
		pv[i,j] = res$p.value
		e[i,j] = res$estimate
	}
}

fdr = matrix(p.adjust(as.vector(pv),"BH"),nrow=ncol(genotype))
colnames(fdr) = colnames(pv)
rownames(fdr) = rownames(pv)


e = log10(e + 0.00001)

cbc_sum = apply(na.omit(cbc),2,sum)


x_wid = nrow(fdr)
y_wid = ncol(fdr)
wid = max(x_wid,y_wid)

# plot
pdf("clinical_genotype_fisher.pdf",width=10,height=10)
wid = 35
plot(NULL, NULL,
xlim=c(-5,wid+5),
ylim=c(-5,wid+5),
axes=FALSE,
xlab="", ylab="",
main=NULL
)

for (i in 1:x_wid){
	
	if(i <= 10){
		x_add = 0
	}else{
		x_add = wid/30
	}

	for (j in 1:y_wid){

		or = e[i,j]
		q = fdr[i,j]
		
		size = 0
		if(or > 0){
			size = (q > 0.1) * 0.1 + ( q < 0.1 & q >= 0.01 ) * 0.3 + ( q < 0.01 & q >= 0.001 ) * 0.6 + ( q < 0.001 ) * 0.99
		}

		color_id = round( or * 50 ) + 51
		if(color_id > 101) {color_id = 101} else if(color_id < 0) {color_id = 1}
		bwr = colorRampPalette(colors=c("blue", "white", "red"))(101)
		color = bwr[color_id]
		
		rect(i-0.5-size/2+x_add,(y_wid-j+1)-0.5-size/2,i-0.5+size/2+x_add,(y_wid-j+1)-0.5+size/2,col=color,border=F)
		rect(i-1+x_add,(y_wid-j+1)-1,i+x_add,(y_wid-j+1),density=0,border="gray",lwd=3.5)
	}
}

n_mut = 11
font_label = c(rep(4,n_mut-1),rep(2,x_wid+1-n_mut))

x_lab = rownames(e)


for (i in 1:x_wid){
	if(i >= 11){
		text(i-1+0.5+wid/30,-0.6,cex=0.7,x_lab[i],adj=1,srt=50,font=font_label[i])
	}else{
		text(i-1+0.5,-0.6,cex=0.7,x_lab[i],adj=1,srt=50,font=font_label[i])
	}
}

y_lab = colnames(e)
#y_lab = c("WBC low","WBC high","Hb low","Hb high","Plt low","Plt high","Cytopenia (All)","Cytopenia (Multi)","Any abnormality")
for (j in 1:y_wid){
	text(-0.4,(y_wid-j+1)-0.5,cex=0.7,y_lab[j],adj=1,srt=0,font=2)
}


# color scale
for (i in 51:101){
	rect(x_wid+3,(i/101)*4+5.3-2,x_wid+3.5,((i+1)/101)*4+5.3-2,col=bwr[i],border=F)
}

text(x_wid+3.25,2+5.3+0.7,cex=0.5,"Odds ratio",adj=0.5)
text(x_wid+4.3,(101/101)*4+5.3-2,cex=0.45,"> 10",adj=0)
text(x_wid+4.3,(51/101)*4+5.3-2,cex=0.45,"< 1",adj=0)
#text(x_wid+4.3,(1/101)*2+5.3,cex=0.45,"< 1",adj=0)

# size scale (positive)
x_loc = x_wid+3.25; leg_col = "gray"
size = 0.1; y_loc = (1/101)*4
rect(x_loc-size/2,y_loc-size/2,x_loc+size/2,y_loc+size/2,col=leg_col,border=F)
text(x_loc+1.5,y_loc,">0.1",cex=0.45,adj=0) 
size = 0.3; y_loc = (20/101)*4
rect(x_loc-size/2,y_loc-size/2,x_loc+size/2,y_loc+size/2,col=leg_col,border=F)
text(x_loc+1.5,y_loc,"<0.1",cex=0.45,adj=0)
size = 0.5; y_loc = (40/101)*4
rect(x_loc-size/2,y_loc-size/2,x_loc+size/2,y_loc+size/2,col=leg_col,border=F)
text(x_loc+1.5,y_loc+0.1,"<0.01",cex=0.45,adj=0)
size = 0.9; y_loc = (65/101)*4
rect(x_loc-size/2,y_loc-size/2,x_loc+size/2,y_loc+size/2,col=leg_col,border=F)
text(x_loc+1.5,y_loc,"<0.001",cex=0.45,adj=0)

text(x_loc-0.06-0.4,4.2,"q",font=3,cex=0.5,adj=1)
text(x_loc+0.06-0.4,4.2,"value",cex=0.5,adj=0)

dev.off()






genotype = cbind(all_mat_mut!=0,all_mat_cna!=0)

TET2 = genotype[,"4:CNN-LOH:q"]
SF3B1 = genotype[,"SF3B1"]
JAK2 = genotype[,"JAK2"]

Hb = dat$Hb
WBC = dat$WBC_ul_WBC_selected
PLT = dat$PLT
MCV =  dat$MCV

b1 = boxplot(Hb[which(Hb>=0 & Hb <30)] ~ TET2[which(Hb>=0 & Hb <30)])
dev.off()
b2 = boxplot(MCV[which(MCV>=75 & MCV <115)] ~ SF3B1[which(MCV>=75 & MCV <115)])
dev.off()
b3 = boxplot(PLT ~ JAK2,ylim=c(0,1000))
dev.off()


mi = min(na.omit(Hb))
ma = max(na.omit(Hb))
mim = mi - (ma - mi)/20
map = ma + (ma - mi)/20
dw = (ma - mi)/20


pdf("clinical_genotype_cbc.pdf",width=30,height=10)
par(mfrow=c(1,3))

plot(
NULL,NULL,
xlim=c(-0.1,2.5),ylim=c(mim-dw,map+2.5*dw),
axes=FALSE,
xlab="",ylab=""
)

# y axis
ax = 0.25
lat = c(5,10,15,20)
lat2 = c((8:9)/2,(11:19)/2,(21:29)/2,(31:39)/2)
for(i in 1:length(lat)){
	segments(ax,lat[i],ax-0.05,lat[i],lwd=3)
	text(ax-0.1,lat[i],lat[i],adj=1,cex=2.5,font=2)
}
for(i in 1:length(lat2)){
	segments(ax,lat2[i],ax-0.02,lat2[i],lwd=3)
}
segments(ax,mim,ax,map,lwd=3)
text(ax-0.37,(map+mim)/2,"Hemoglobin (g/dL)",srt=90,font=2,cex=4)

le = 0.07 #0.05
md = 0.20
re = 0.33 #0.35

b = b1
x0 = 0.9
x1 = 1.9

xs = c(x0+le,x0+le,x0+md,x0+le,x0+le,x0+re,x0+re,x0+md,x0+re,x0+re)
ys = c(b$stats[2,1],(b$stats[2,1]+3*b$stats[3,1])/4,b$stats[3,1],(b$stats[4,1]+3*b$stats[3,1])/4,b$stats[4,1])
ys = c(ys,rev(ys))
polygon(xs,ys,col="gray",border="black",lwd=2)
arrows(x0+md,b$stats[4,1],x0+md,b$stats[5,1],angle=90,lwd=2)
arrows(x0+md,b$stats[2,1],x0+md,b$stats[1,1],angle=90,lwd=2)

xs = c(x1+le,x1+le,x1+md,x1+le,x1+le,x1+re,x1+re,x1+md,x1+re,x1+re)
ys = c(b$stats[2,2],(b$stats[2,2]+3*b$stats[3,2])/4,b$stats[3,2],(b$stats[4,2]+3*b$stats[3,2])/4,b$stats[4,2])
ys = c(ys,rev(ys))
polygon(xs,ys,col="gray",border="black",lwd=2)
arrows(x1+md,b$stats[4,2],x1+md,b$stats[5,2],angle=90,lwd=2)
arrows(x1+md,b$stats[2,2],x1+md,b$stats[1,2],angle=90,lwd=2)

Hb_TET2_n = na.omit(Hb[which(TET2==0)])
Hb_TET2_p = na.omit(Hb[which(TET2==1)])

points(runif(length(Hb_TET2_n),x0-re,x0-le),Hb_TET2_n,cex=2,pch=21,col=NULL,bg=rgb(0.5,0,0.5, alpha=0.3))
points(runif(length(Hb_TET2_p),x1-re,x1-le),Hb_TET2_p,cex=2,pch=21,col=NULL,bg=rgb(0.5,0,0.5, alpha=0.3))

# genotype
text(x0-0.1,mim-dw,"TET2",font=4,cex=4)
text(x0+0.25,mim-dw,"(-)",font=2,cex=4)
text(x1-0.1,mim-dw,"TET2",font=4,cex=4)
text(x1+0.25,mim-dw,"(+)",font=2,cex=4)

# significance
#segments(x0,map+1.5*dw,x1,map+1.5*dw)
#segments(x0,map+1.5*dw,x0,map+1*dw)
#segments(x1,map+1.5*dw,x1,map+1*dw)
#text((x0+x1)/2,map+2*dw,"*",cex=4.5)


# MCV on SF3B1 mutations
#mi = min(na.omit(MCV))
#ma = max(na.omit(MCV))
mi = 75
ma = 115
mim = mi - (ma - mi)/20
map = ma + (ma - mi)/20
dw = (ma - mi)/20 

plot(
NULL,NULL,
xlim=c(-0.1,2.5),ylim=c(mim-dw,map+2.5*dw),
axes=FALSE,
xlab="",ylab=""
)

# y axis
ax = 0.25
lat = c(80,90,100,110)
lat2 = c(70+(3:4)*2,80+(1:4)*2,90+(1:4)*2,100+(1:4)*2,110+(1:4)*2)
for(i in 1:length(lat)){
	segments(ax,lat[i],ax-0.05,lat[i],lwd=3)
	text(ax-0.1,lat[i],lat[i],adj=1,cex=2.5,font=2)
}
for(i in 1:length(lat2)){
	segments(ax,lat2[i],ax-0.02,lat2[i],lwd=3)
}
segments(ax,mim,ax,map,lwd=3)
text(ax-0.37,(map+mim)/2,"MCV (fL)",srt=90,cex=4,font=2)

le = 0.07
md = 0.20
re = 0.33

b = b2
x0 = 0.9
x1 = 1.9

xs = c(x0+le,x0+le,x0+md,x0+le,x0+le,x0+re,x0+re,x0+md,x0+re,x0+re)
ys = c(b$stats[2,1],(b$stats[2,1]+3*b$stats[3,1])/4,b$stats[3,1],(b$stats[4,1]+3*b$stats[3,1])/4,b$stats[4,1])
ys = c(ys,rev(ys))
polygon(xs,ys,col="gray",border="black",lwd=2)
arrows(x0+md,b$stats[4,1],x0+md,b$stats[5,1],angle=90,lwd=2)
arrows(x0+md,b$stats[2,1],x0+md,b$stats[1,1],angle=90,lwd=2)

xs = c(x1+le,x1+le,x1+md,x1+le,x1+le,x1+re,x1+re,x1+md,x1+re,x1+re)
ys = c(b$stats[2,2],(b$stats[2,2]+3*b$stats[3,2])/4,b$stats[3,2],(b$stats[4,2]+3*b$stats[3,2])/4,b$stats[4,2])
ys = c(ys,rev(ys))
polygon(xs,ys,col="gray",border="black",lwd=2)
arrows(x1+md,b$stats[4,2],x1+md,b$stats[5,2],angle=90,lwd=2)
arrows(x1+md,b$stats[2,2],x1+md,b$stats[1,2],angle=90,lwd=2)

SF3B1_2 = SF3B1[which(MCV>=75 & MCV < 115)]
MCV_2 = MCV[which(MCV>=75 & MCV < 115)]

MCV_SF3B1_n = na.omit(MCV_2[which(SF3B1_2==0)])
MCV_SF3B1_p = na.omit(MCV_2[which(SF3B1_2==1)])

points(runif(length(MCV_SF3B1_n),x0-re,x0-le),MCV_SF3B1_n,cex=2,pch=21,col=NULL,bg=rgb(0.5,0,0.5, alpha=0.3))
points(runif(length(MCV_SF3B1_p),x1-re,x1-le),MCV_SF3B1_p,cex=2,pch=21,col=NULL,bg=rgb(0.5,0,0.5, alpha=0.3))

# genotype
text(x0-0.1,mim-dw,"SF3B1",font=4,cex=4)
text(x0+0.25,mim-dw,"(-)",font=2,cex=4)
text(x1-0.1,mim-dw,"SF3B1",font=4,cex=4)
text(x1+0.25,mim-dw,"(+)",font=2,cex=4)

# significance
#segments(x0,map+1.5*dw,x1,map+1.5*dw)
#segments(x0,map+1.5*dw,x0,map+1*dw)
#segments(x1,map+1.5*dw,x1,map+1*dw)
#text((x0+x1)/2,map+2*dw,"**",cex=4.5)


# PLT on JAK2
mi = min(na.omit(PLT))
#ma = max(na.omit(PLT))

ma = 150
mim = mi - (ma - mi)/20
map = ma + (ma - mi)/20
dw = (ma - mi)/20 

plot(
NULL,NULL,
xlim=c(-0.1,2.5),ylim=c(mim-dw,map+2.5*dw),
axes=FALSE,
xlab="",ylab=""
)

# y axis
ax = 0.25
lat = c(0,20,40,60,80,100,120,140)
lat2 = (1:28)*5
for(i in 1:length(lat)){
	segments(ax,lat[i],ax-0.05,lat[i],lwd=3)
	text(ax-0.1,lat[i],lat[i],adj=1,cex=2.5,font=2)
}
for(i in 1:length(lat2)){
	segments(ax,lat2[i],ax-0.02,lat2[i],lwd=3)
}
segments(ax,mim,ax,map,lwd=3)
text(ax-0.37,(map+mim)/2,"Platelet (10^4/uL)",srt=90,font=2,cex=4)

le = 0.07
md = 0.20
re = 0.33

b = b3
x0 = 0.9
x1 = 1.9

xs = c(x0+le,x0+le,x0+md,x0+le,x0+le,x0+re,x0+re,x0+md,x0+re,x0+re)
ys = c(b$stats[2,1],(b$stats[2,1]+3*b$stats[3,1])/4,b$stats[3,1],(b$stats[4,1]+3*b$stats[3,1])/4,b$stats[4,1])
ys = c(ys,rev(ys))
polygon(xs,ys,col="gray",border="black",lwd=2)
arrows(x0+md,b$stats[4,1],x0+md,b$stats[5,1],angle=90,lwd=2)
arrows(x0+md,b$stats[2,1],x0+md,b$stats[1,1],angle=90,lwd=2)

xs = c(x1+le,x1+le,x1+md,x1+le,x1+le,x1+re,x1+re,x1+md,x1+re,x1+re)
ys = c(b$stats[2,2],(b$stats[2,2]+3*b$stats[3,2])/4,b$stats[3,2],(b$stats[4,2]+3*b$stats[3,2])/4,b$stats[4,2])
ys = c(ys,rev(ys))
polygon(xs,ys,col="gray",border="black",lwd=2)
arrows(x1+md,b$stats[4,2],x1+md,b$stats[5,2],angle=90,lwd=2)
arrows(x1+md,b$stats[2,2],x1+md,b$stats[1,2],angle=90,lwd=2)

# genotype
text(x0-0.1,mim-dw,"JAK2",font=4,cex=4)
text(x0+0.25,mim-dw,"(-)",font=2,cex=4)
text(x1-0.1,mim-dw,"JAK2",font=4,cex=4)
text(x1+0.25,mim-dw,"(+)",font=2,cex=4)

# significance
#segments(x0,map+1.5*dw,x1,map+1.5*dw)
#segments(x0,map+1.5*dw,x0,map+1*dw)
#segments(x1,map+1.5*dw,x1,map+1*dw)
#text((x0+x1)/2,map+2*dw,"***",cex=4.5)

JAK2_2 = JAK2[which(PLT>=0 & PLT <150)]
PLT_2 = PLT[which(PLT>=0 & PLT <150)]

PLT_JAK2_n = na.omit(PLT_2[which(JAK2_2==0)])
PLT_JAK2_p = na.omit(PLT_2[which(JAK2_2==1)])

points(runif(length(PLT_JAK2_n),x0-re,x0-le),PLT_JAK2_n,cex=2,pch=21,col=NULL,bg=rgb(0.5,0,0.5, alpha=0.3))
points(runif(length(PLT_JAK2_p),x1-re,x1-le),PLT_JAK2_p,cex=2,pch=21,col=NULL,bg=rgb(0.5,0,0.5, alpha=0.3))

dev.off()




## TET2 4qUPD ##

genotype = cbind(all_mat_mut_2!=0,all_mat_cna_2!=0)

mut_path = "MatMed_2021_mut_table.txt"
call = read.table(mut_path,header=T,stringsAsFactor=F,sep="\t",quote="",colClass="character")
call = call[which(is.element(call$id,all_id)),]
call = call[which(call$Gene.refGene == "TET2"),]

fac_none = as.numeric(apply(genotype,1,sum)==0)
fac_tet2 = as.numeric(genotype[,"TET2"]!=0)
fac_upd = as.numeric(genotype[,"4:CNN-LOH:q"] != 0)
fac_del = as.numeric(genotype[,"4:Deletion:q"] != 0)
fac_loh = as.numeric(genotype[,"4:CNN-LOH:q"] + genotype[,"4:Deletion:q"] + genotype[,"4:Duplication:q"] + genotype[,"4:Unknown:q"] != 0)
cbc = dat[which(fac_none+fac_tet2+fac_loh!=0),]


all_cases = all_id[which(fac_none+fac_tet2+fac_loh!=0)]
upd_cases = all_id[which(fac_loh==1)]
upd_cases_2 = all_id[which(fac_upd==1)]
del_cases = all_id[which(fac_del==1)]

upd_cases
upd_cases_2

Hbs = dat$Hb
names(Hbs) = all_id
Hbs = Hbs[all_cases]
names(Hbs) = all_cases
Hbs
if(F){
	library(car)
		nna = which(!is.na(cbc[,"Hb_g_per_dl_Hb_selected"]))
		v = na.omit(cbc[,"Hb_g_per_dl_Hb_selected"])
		v1 = powerTransform(v)
		v2 = na.omit(bcPower(v, v1$lambda))
		print(shapiro.test(sample(v2,1000,replace=F)))
		cbc[nna,"Hb_g_per_dl_Hb_selected"] = v2
}


res = as.data.frame(matrix(0,nrow=length(all_cases),ncol=9))
rownames(res) = all_cases
colnames(res) = c("UPD","UPD2","del","VAF","read","call","plt","sex","age")


for(i in 1:nrow(res)){
	cur_id = all_cases[i]
	res[i,"UPD"] = as.numeric(is.element(cur_id,upd_cases))
	res[i,"UPD2"] = as.numeric(is.element(cur_id,upd_cases_2))
	res[i,"del"] = as.numeric(is.element(cur_id,del_cases))
	res[i,"sex"] = as.numeric(all_sex[cur_id]==1)
	res[i,"age"] = as.numeric(all_age[cur_id])
	res[i,"plt"] = as.numeric(Hbs[cur_id])
	
	if(sum(call$id == cur_id) == 0){
		next
	}

	res[i,"VAF"] = max(call$misRate[which(call$id == cur_id)])
	#res[i,"VAF"] = sum(call$misRate[which(call$id == cur_id)])
	res[i,"read"] = max(call$variantNum[which(call$id == cur_id)]) 

	if(is.element(cur_id ,call$id)){
		res[i,"call"] = sum(call$id == cur_id)
		#res[i,"call"] = sum(call$misRate[which(call$id == cur_id)] >= 0.05)
	}
}

res

mm = min(na.omit(as.numeric(res[which(res[,"UPD"] == 1),"plt"])))
mm
rownames(res)[which(res[,"plt"] == mm)]

print("none vs AI")
res_tet2 = as.data.frame(res[which(res[,"call"] == 0 | res[,"UPD"] == 1),])
summary(lm(plt ~ UPD + age + sex,data=res_tet2))

print("none vs UPD")
res_tet2 = as.data.frame(res[which(res[,"call"] == 0 | res[,"UPD2"] == 1),])
summary(lm(plt ~ UPD2 + age + sex,data=res_tet2))

print("none vs 1")
res_tet2 = as.data.frame(res[which(res[,"call"] <= 1 & res[,"UPD"] == 0),])
res_tet2$call = as.numeric(res_tet2$call >= 1)
summary(lm(plt ~ call + age + sex,data=res_tet2))

print("none vs 2")
res_tet2 = as.data.frame(res[which((res[,"call"] >= 2 | res[,"call"] >= 0) & res[,"UPD"] == 0),])
res_tet2$call = as.numeric(res_tet2$call >= 2)
summary(lm(plt ~ call + age + sex,data=res_tet2))

print("1 vs 2")
res_tet2 = as.data.frame(res[which((res[,"call"] >= 1) & res[,"UPD"] == 0),])
res_tet2$call = as.numeric(res_tet2$call >= 2)
summary(lm(plt ~ call + age + sex,data=res_tet2))

print("1 vs AI")
res_tet2 = as.data.frame(res[which(res[,"call"] == 1 | res[,"UPD"] == 1),])
summary(lm(plt ~ UPD + age + sex,data=res_tet2))

print("1 vs UPD")
res_tet2 = as.data.frame(res[which(res[,"call"] == 1 | res[,"UPD2"] == 1),])
summary(lm(plt ~ UPD2 + age + sex,data=res_tet2))

print("2 vs AI")
res_tet2 = as.data.frame(res[which(res[,"call"] >= 2 | res[,"UPD"] == 1),])
summary(lm(plt ~ UPD + age + sex,data=res_tet2))

print("2 vs UPD")
res_tet2 = as.data.frame(res[which(res[,"call"] >= 2 | res[,"UPD2"] == 1),])
summary(lm(plt ~ UPD2 + age + sex,data=res_tet2))

print("none vs 3")
res_tet2 = as.data.frame(res[which(res[,"call"] >= 3 | res[,"call"] == 0 | res[,"UPD"] == 0),])
res_tet2$call = as.numeric(res_tet2$call >= 3)
summary(lm(plt ~ call + age + sex,data=res_tet2))

res_tet2 = as.data.frame(res[which(res[,"UPD"] == 1),])
res_tet2$AI = 1
res_tet3 = as.data.frame(res[which(res[,"UPD2"] == 1),])
res_tet3$AI = 0
res_tet4 = rbind(res_tet2,res_tet3)
summary(lm(plt ~ AI + age + sex,data=res_tet4))


print("n of mut")
res_tet2 = as.data.frame(res[which(res[,"call"] >= 0 | res[,"UPD"] == 0),])
summary(lm(plt ~ call + age + sex,data=res_tet2))

#q()

hb_none = na.omit(res[which(res[,"call"] == 0 & res[,"UPD"] == 0),"plt"])
hb_mut_1 = na.omit(res[which(res[,"call"] == 1 & res[,"UPD"] == 0),"plt"])
hb_mut_2 = na.omit(res[which(res[,"call"] >= 2 & res[,"UPD"] == 0),"plt"])
hb_upd = na.omit(res[which(res[,"UPD2"] == 1),"plt"])
hb_ai = na.omit(res[which(res[,"UPD"] == 1),"plt"])

box = function(x){
	q = quantile(x)
	return(c(q[2]-1.5*(q[4]-q[2]),q[2],q[3],q[4],q[4]+1.5*(q[4]-q[2])))
}

box_none = box(hb_none)
box_mut_1 = box(hb_mut_1)
box_mut_2 = box(hb_mut_2)
box_upd = box(hb_upd)
box_ai = box(hb_ai)

x_none = 1*1.7 -1 + runif(length(hb_none),-0.2,0.2)
x_mut_1 = 2*1.7 -1 + runif(length(hb_mut_1),-0.2,0.2)
x_mut_2 = 3*1.7 -1 + runif(length(hb_mut_2),-0.2,0.2)
x_upd = 4*1.7 -1 + runif(length(hb_upd),-0.2,0.2)
x_ai = 5*1.7 -1 + runif(length(hb_ai),-0.2,0.2)

m = max(hb_none,hb_mut_1,hb_mut_2,hb_ai)
l = min(hb_none,hb_mut_1,hb_mut_2,hb_ai)
s = 7/(m-2)
print(m)
print(l)

pdf("TET2_4pUPD_Hb.pdf",width=12,height=10)

plot(
NULL,NULL,
xlim = c(-1,10),
ylim = c(2,m*(6/5)),
xlab = "",
ylab = "",
axes = F
)

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

# text(-0.2,(min(lats)+max(l,"VAF of V617F (%)",srt=90,adj=0.5)
theta = seq(-pi,pi,length=100)
r = 0.1
for(i in 1:length(x_none)){
	polygon(x_none[i]+s*cos(theta)*r,hb_none[i]+sin(theta)*r,col=rgb(0.5,0,0.5, alpha=0.3),border="transparent")
}
for(i in 1:length(x_mut_1)){
	polygon(x_mut_1[i]+s*cos(theta)*r,hb_mut_1[i]+sin(theta)*r,col=rgb(0.5,0,0.5, alpha=0.3),border="transparent")
}
for(i in 1:length(x_mut_2)){
	polygon(x_mut_2[i]+s*cos(theta)*r,hb_mut_2[i]+sin(theta)*r,col=rgb(0.5,0,0.5, alpha=0.3),border="transparent")
}
for(i in 1:length(x_upd)){
	polygon(x_upd[i]+s*cos(theta)*r,hb_upd[i]+sin(theta)*r,col=rgb(0.5,0,0.5, alpha=0.3),border="transparent")
}
for(i in 1:length(x_ai)){
	polygon(x_ai[i]+s*cos(theta)*r,hb_ai[i]+sin(theta)*r,col=rgb(0.5,0,0.5, alpha=0.3),border="transparent")
}


b = box_none
x0 = 1*1.7-0.5; w0 = 0.15; xs = c(x0-w0,x0-w0,x0,x0-w0,x0-w0,x0+w0,x0+w0,x0,x0+w0,x0+w0)
ys = c(b[2],(b[2]+3*b[3])/4,b[3],(b[4]+3*b[3])/4,b[4])
ys = c(ys,rev(ys))
polygon(xs,ys,border="black",lwd=2,col="gray")
arrows(x0,b[4],x0,max(hb_none[which(hb_none <= b[5])]),angle=90,lwd=2,length=w0)
arrows(x0,b[2],x0,min(hb_none[which(hb_none >= b[1])]),angle=90,lwd=2,length=w0)
# arrows(x0,b[4],x0,b[5],angle=90,lwd=2,length=w0)
# arrows(x0,b[2],x0,b[1],angle=90,lwd=2,length=w0)

b = box_mut_1
x0 = 2*1.7-0.5; w0 = 0.15; xs = c(x0-w0,x0-w0,x0,x0-w0,x0-w0,x0+w0,x0+w0,x0,x0+w0,x0+w0)
ys = c(b[2],(b[2]+3*b[3])/4,b[3],(b[4]+3*b[3])/4,b[4])
ys = c(ys,rev(ys))
polygon(xs,ys,border="black",lwd=2,col="gray")
arrows(x0,b[4],x0,max(hb_mut_1[which(hb_mut_1 <= b[5])]),angle=90,lwd=2,length=w0)
arrows(x0,b[2],x0,min(hb_mut_1[which(hb_mut_1 >= b[1])]),angle=90,lwd=2,length=w0)
# arrows(x0,b[4],x0,b[5],angle=90,lwd=2,length=w0)
# arrows(x0,b[2],x0,b[1],angle=90,lwd=2,length=w0)

b = box_mut_2
x0 = 3*1.7-0.5; w0 = 0.15; xs = c(x0-w0,x0-w0,x0,x0-w0,x0-w0,x0+w0,x0+w0,x0,x0+w0,x0+w0)
ys = c(b[2],(b[2]+3*b[3])/4,b[3],(b[4]+3*b[3])/4,b[4])
ys = c(ys,rev(ys))
polygon(xs,ys,border="black",lwd=2,col="gray")
arrows(x0,b[4],x0,max(hb_mut_2[which(hb_mut_2 <= b[5])]),angle=90,lwd=2,length=w0)
arrows(x0,b[2],x0,min(hb_mut_2[which(hb_mut_2 >= b[1])]),angle=90,lwd=2,length=w0)
# arrows(x0,b[4],x0,b[5],angle=90,lwd=2,length=w0)
# arrows(x0,b[2],x0,b[1],angle=90,lwd=2,length=w0)

b = box_upd
x0 = 4*1.7-0.5; w0 = 0.15; xs = c(x0-w0,x0-w0,x0,x0-w0,x0-w0,x0+w0,x0+w0,x0,x0+w0,x0+w0)
ys = c(b[2],(b[2]+3*b[3])/4,b[3],(b[4]+3*b[3])/4,b[4])
ys = c(ys,rev(ys))
polygon(xs,ys,border="black",lwd=2,col="gray")
arrows(x0,b[4],x0,max(hb_upd[which(hb_upd <= b[5])]),angle=90,lwd=2,length=w0)
arrows(x0,b[2],x0,min(hb_upd[which(hb_upd >= b[1])]),angle=90,lwd=2,length=w0)
#arrows(x0,b[4],x0,b[5],angle=90,lwd=2,length=w0)
#arrows(x0,b[2],x0,b[1],angle=90,lwd=2,length=w0)

b = box_ai
x0 = 5*1.7-0.5; w0 = 0.15; xs = c(x0-w0,x0-w0,x0,x0-w0,x0-w0,x0+w0,x0+w0,x0,x0+w0,x0+w0)
ys = c(b[2],(b[2]+3*b[3])/4,b[3],(b[4]+3*b[3])/4,b[4])
ys = c(ys,rev(ys))
polygon(xs,ys,border="black",lwd=2,col="gray")
arrows(x0,b[4],x0,max(hb_ai[which(hb_ai <= b[5])]),angle=90,lwd=2,length=w0)
arrows(x0,b[2],x0,min(hb_ai[which(hb_ai >= b[1])]),angle=90,lwd=2,length=w0)
#arrows(x0,b[4],x0,b[5],angle=90,lwd=2,length=w0)
#arrows(x0,b[2],x0,b[1],angle=90,lwd=2,length=w0)

text(1*1.7-0.75,3,"No alteration",adj=0.5,font=2)
text(2*1.7-0.75,3,"One mutation",adj=0.5,font=2)
text(3*1.7-0.75,3,">=2 mutations",adj=0.5,font=2)
text(4*1.7-0.75,3,"4qUPD",adj=0.5,font=2)
text(5*1.7-0.75,3,"4qAI",adj=0.5,font=2)

dev.off()






