
mut_path = "MatMed_2021_mut_table.txt"
cna_path = "MatMed_2021_cna_table.csv"  
id = as.character(1:11234)

## get centromere position ##
arminfo = read.table("cytoband/arm_chr_length.txt")
cent = arminfo[1:22,4]


getCnaLabel2 = function(s){
	splited = unlist(strsplit(s,"_"))
	type=NA
	
	if(length(splited)==1){
		return(s)
	}
	
	if(splited[1]=="cna"){
		return(splited[2])
	}
	prf = ifelse(splited[2]=="CNN-LOH","",ifelse(splited[2]=="Duplication","+",ifelse(splited[2]=="Deletion","del(","Unknown")))
	suf = ifelse(splited[2]=="CNN-LOH","UPD",ifelse(splited[2]=="Duplication","",ifelse(splited[2]=="Deletion",")","Unknown")))
	return(paste0(prf,splited[1],splited[3],suf))
}


mut = read.table(mut_path,header=T,stringsAsFactor=F,sep="\t",quote="",colClass="character")
mut = mut[which(!is.na(mut$Chr)),]

mat_mut = matrix(0,ncol=length(unique(mut$Gene.refGene)),nrow=length(id))
colnames(mat_mut) = unique(mut$Gene.refGene); rownames(mat_mut) = id

for(i in 1:nrow(mut)){
	cur_id = mut[i,1]
	cur_gene = mut$Gene.refGene[i]
	cur_vaf = as.numeric(mut$misRate[i])
	if(! is.element(cur_id,id)){
		next
	}
	if(cur_vaf > mat_mut[cur_id,cur_gene]){
		mat_mut[cur_id,cur_gene] = cur_vaf
	}
}

cna = read.csv(cna_path,header=T,stringsAsFactor=F,quote="",colClass="character")
cna = cna[which(!is.na(cna$chr)),]
label_cna = paste0(sort(rep(1:22,4*2)),"_",sort(rep(c("CNN-LOH","Duplication","Deletion","unknown"),2)),"_",c("p","q"))
mat_cna =  matrix(0,ncol=length(label_cna),nrow=length(id))
colnames(mat_cna) = label_cna; rownames(mat_cna) = id

num_14q_1 = which(as.numeric(cna$chr) == 14 & cna$COPY_CHANGE == "loss" & cna$start <= 25)
num_14q_2 = which(as.numeric(cna$chr) == 14 & cna$COPY_CHANGE == "loss" & cna$start >= 25)
num_14q_1
num_14q_2

for(i in 1:nrow(cna)){
	cur_id = cna$id[i]
	if(! is.element(cur_id,id)){
		next
	}

	if(is.element(i,num_14q_2)){
		next
	}
	
	cur_type = "unknown"
	if(cna$COPY_CHANGE[i] == "loss"){
		cur_type = "Deletion"
	}else if(cna$COPY_CHANGE[i] == "neutral"){
		cur_type = "CNN-LOH"
	}else if(cna$COPY_CHANGE[i] == "gain"){
		cur_type = "Duplication"
	}
	cur_cna = paste0(cna$chr[i],"_",cur_type)

	cur_cent = cent[as.numeric(cna$chr[i])]
	if(as.numeric(cna$start[i])*(10^6) <= cur_cent & as.numeric(cna$end[i])*(10^6) <= cur_cent){
		cur_cna = paste0(cur_cna,"_p")
	}else if(as.numeric(cna$start[i])*(10^6) >= cur_cent & as.numeric(cna$end[i])*(10^6) >= cur_cent){
		cur_cna = paste0(cur_cna,"_q")
	}else if(as.numeric(cna$start[i])*(10^6) >= cur_cent & as.numeric(cna$end[i])*(10^6) >= cur_cent){
		cur_cna = c(paste0(cur_cna,"_p"),paste0(cur_cna,"_q"))
	}
	if(! is.element(cur_cna,label_cna)){
		next
	}
	if(cur_type == "unknown"){
		cur_cf = 0.00001
	}
	cur_cf = cna$CELL_FRAC[i]
	if(mat_cna[cur_id,cur_cna] < cur_cf & !is.na(cur_cf)){
		mat_cna[cur_id,cur_cna] = cur_cf
	}
}

mut_mat = (mat_mut != 0)
cna_mat = (mat_cna != 0)

cna_2pLOH = cna_mat[,"2_CNN-LOH_p"] | cna_mat[,"2_Deletion_p"] | cna_mat[,"2_unknown_p"] 
cna_4qLOH = cna_mat[,"4_CNN-LOH_q"] | cna_mat[,"4_Deletion_q"] | cna_mat[,"4_unknown_q"]
cna_17pLOH = cna_mat[,"17_CNN-LOH_p"] | cna_mat[,"17_Deletion_p"] | cna_mat[,"17_unknown_p"]
ex = c("2_CNN-LOH_p","2_Deletion_p","4_CNN-LOH_q","4_Deletion_q","17_CNN-LOH_p","17_Deletion_p","2_unknown_p","4_unknown_q","17_unknown_p")
cna_mat = cbind(cna_2pLOH,cna_4qLOH,cna_17pLOH,cna_mat[,-which(is.element(colnames(cna_mat),ex))])

colnames(cna_mat)

cna_result_filt = cna_mat[,which(apply(cna_mat,2,sum)>=15)]
mut_result_filt = mut_mat[,which(apply(mut_mat,2,sum)>=20)]

colnames(cna_result_filt)

colnames(mut_result_filt)[which(colnames(mut_result_filt)=="U2AF1;U2AF1L5")] = "U2AF1"
mut_result_filt = mut_result_filt[,order(apply(mut_result_filt,2,sum),decreasing=F)]

ordered = c(
"cna_2pLOH",
"cna_4qLOH",
"cna_17pLOH",
"1_CNN-LOH_p",
"1_CNN-LOH_q",
"6_CNN-LOH_p",
"9_CNN-LOH_p",
"9_CNN-LOH_q",
"11_CNN-LOH_q",
"14_CNN-LOH_q",
"16_CNN-LOH_p",
"17_CNN-LOH_q",
"5_Deletion_q",
"6_Deletion_q",
"11_Deletion_q",
"13_Deletion_q",
"14_Deletion_q",
"20_Deletion_q",
"15_Duplication_q",
"21_Duplication_q",
"22_Duplication_q")

order_label = c(
"2pLOH","4qLOH","17pLOH",
"1pUPD","1qUPD",
"6pUPD","9pUPD",
"9qUPD",
"11qUPD",
"14qUPD","16pUPD",
"17qUPD",
"del(5q)","del(6q)",
"del(11q)","del(13q)",
"del(14q)","del(20q)",
"+15q",
"+21q","+22q"
)


# odd_ratio
odds_ratio <- function(a, b, c, d, correct=FALSE){

  if (correct | sign(a)*sign(b)*sign(c)*sign(d) == 0) {
    a <- a+0.5
    b <- b+0.5
    c <- c+0.5
    d <- d+0.5
  }
  or <- (a/b)*(d/c)
  SE <- sqrt(1/a + 1/b + 1/c + 1/d)
  q <- qnorm(0.975)
  ci_95l <- or * exp(-q * SE)
  ci_95u <- or * exp(q * SE)
  x<- pnorm(log(or)/SE)
  if(x<=0.5){
    p<- x*2
  }else{
    p<- (1-x)*2
  }
  c(or=or, p=p, ci_95l=ci_95l, ci_95u=ci_95u)
}



# figure 1
# make input
mut_result_filt_2 = cna_result_filt # cna_mat[,which(apply(cna_mat,2,sum)>=10)]
cna_result_filt_2 = mut_result_filt # mut_mat[,which(apply(mut_mat,2,sum)>=10)]

rectanCor_or = rectanCor_p = matrix(NA,ncol=ncol(cna_result_filt_2),nrow=ncol(mut_result_filt_2))
colnames(rectanCor_or) = colnames(rectanCor_p) = colnames(cna_result_filt_2)
rownames(rectanCor_or) = rownames(rectanCor_p) = colnames(mut_result_filt_2)

for (i in 1:ncol(cna_result_filt_2)){
	for (j in 1:ncol(mut_result_filt_2)){
		pp = sum((cna_result_filt_2[,i]==1)*(mut_result_filt_2[,j]==1))
		pn = sum((cna_result_filt_2[,i]==1)*(mut_result_filt_2[,j]!=1))
		np = sum((cna_result_filt_2[,i]!=1)*(mut_result_filt_2[,j]==1))
		nn = sum((cna_result_filt_2[,i]!=1)*(mut_result_filt_2[,j]!=1))
		rectanCor_or[j,i] = odds_ratio(pp,pn,np,nn)["or"]
		rectanCor_p[j,i]  = odds_ratio(pp,pn,np,nn)["p"]
	}
}

rectanCor_p


res = matrix(NA,ncol=5,nrow=ncol(cna_result_filt_2)*ncol(mut_result_filt_2))
k = 1
for (i in 1:ncol(cna_result_filt_2)){
	for (j in 1:ncol(mut_result_filt_2)){
		pp = sum((cna_result_filt_2[,i]==1)*(mut_result_filt_2[,j]==1))
    pn = sum((cna_result_filt_2[,i]==1)*(mut_result_filt_2[,j]!=1))
    np = sum((cna_result_filt_2[,i]!=1)*(mut_result_filt_2[,j]==1))
    nn = sum((cna_result_filt_2[,i]!=1)*(mut_result_filt_2[,j]!=1))
		res[k,] = c(colnames(mut_result_filt_2)[j],colnames(cna_result_filt_2)[i],odds_ratio(pp,pn,np,nn)[c("or","p")],pp)
		k = k + 1
	}
}

P = as.numeric(res[,4])
Q = p.adjust(P,"BH")
res[,4] = Q

res_sel = res[which(as.numeric(res[,4])<0.1),]
mut_lab = unique(res_sel[,1])
cna_lab = unique(res_sel[,1])

x_wid = nrow(rectanCor_or)
y_wid = ncol(rectanCor_or)
wid = max(x_wid,y_wid)

pdf("rectan_cor_1.pdf",width=9.5,height=10)

plot(NULL, NULL,
xlim=c(-3,wid+3+10),ylim=c(-3,wid+3+10),
axes=FALSE,xlab="", ylab="",main=NULL)

for (i in 1:x_wid){
	for (j in 1:y_wid){
		j2 = y_wid - j
		
		n_mut = rownames(rectanCor_p)[i]
		n_cna = colnames(rectanCor_p)[j]
		c = intersect(which(res[,1] == n_mut),which(res[,2] == n_cna))
		
		or = as.numeric(res[c,3])
		q = as.numeric(res[c,4])
		n_pp = as.numeric(res[c,5])
		
		size = ( q >= 0.1 ) * 0.1 +
			( q < 0.1 & q >= 0.01 ) * 0.3 + 
			( q < 0.01 & q >= 0.001 ) * 0.6 + 
			( q < 0.001 ) * 1

		color_id = round( log10(or) * 25 ) + 51
		if(color_id > 101){color_id = 101}else if(color_id < 0) {color_id = 1}
		bwr = colorRampPalette(colors=c("blue", "white", "red"))(101)
		color = bwr[color_id]
	
		rect(i-0.5-size/2,(y_wid-j2+1)-0.5-size/2,i-0.5+size/2,(y_wid-j2+1)-0.5+size/2,col=color,border=F)
		rect(i-1,(y_wid-j2+1)-1,i,(y_wid-j2+1),density=0,border="gray",lwd=3)
		
		if(n_pp >= 5){
			segments(i-0.5,(y_wid-j2+1)-0.5-0.15,i-0.5,(y_wid-j2+1)-0.5+0.15,lwd=0.5)
			segments(i-0.5-0.1,(y_wid-j2+1)-0.5+0.1,i-0.5+0.1,(y_wid-j2+1)-0.5-0.1,lwd=0.5)
			segments(i-0.5-0.1,(y_wid-j2+1)-0.5-0.1,i-0.5+0.1,(y_wid-j2+1)-0.5+0.1,lwd=0.5)
		}
	}
}

font_label = c(rep(3,13),rep(1,x_wid-13))

for (i in 1:x_wid){
	text(i-1+0.6,1-0.5,cex=0.65,getCnaLabel2(rownames(rectanCor_p)[i]),adj=1,srt=45,font=2)
}

for (j in 1:y_wid){
	j2 = y_wid - j
	text(-0.4,(y_wid-j2+1)-0.5,cex=0.65,colnames(rectanCor_p)[j],adj=1,font=4)
}

dev.off()

pdf("rectan_cor_1_legend.pdf",width=9.5,height=10)

plot(NULL, NULL,
xlim=c(-3,wid+3+10),ylim=c(-3,wid+3+10),
axes=FALSE,xlab="", ylab="",main=NULL)

legend_x = 15

# size scale (positive)
size = 0.3
x_loc = legend_x-4.5; y_loc = y_wid*(7/6)
rect(x_loc-size/2,y_loc-size/2,x_loc+size/2,y_loc+size/2,col="gray",border=F)
text(x_loc+1,y_loc,"< 0.1",cex=0.7,font=2,adj=0)
size = 0.6
y_loc = y_loc-1.5
rect(x_loc-size/2,y_loc-size/2,x_loc+size/2,y_loc+size/2,col="gray",border=F)
text(x_loc+1,y_loc,"< 0.01",cex=0.7,font=2,adj=0)
size = 0.9
y_loc = y_loc-2
rect(x_loc-size/2,y_loc-size/2,x_loc+size/2,y_loc+size/2,col="gray",border=F)
text(x_loc+1,y_loc,"< 0.001",cex=0.7,font=2,adj=0)

text(legend_x-4.5,y_loc+5,"q",cex=0.7,font=4,adj=0.5)
text(legend_x-3,y_loc+5,"value",cex=0.7,font=2,adj=0.5)


# odds ratio
text(legend_x+4.5,y_loc+5,cex=0.7,"Odds ratio",font=2,adj=0.5)
text(legend_x+5,y_loc+0,cex=0.7,"< 0.01",font=2,adj=0)
text(legend_x+5,y_loc+1.5,cex=0.7,"1",font=2,adj=0)
text(legend_x+5,y_loc+3,cex=0.7,"> 100",font=2,adj=0)

# color scale
for (i in 1:101){
	rect(legend_x+3-0.3,y_loc+(i-1)/(101/3),legend_x+3+0.3,y_loc+i/(101/3),col=rev(bwr)[i],border=F)
}

dev.off()


## Correlation of all combinations ##
## Draw ED Fig.2d ##
# make input
cna_result_filt = cna_mat[,which(apply(cna_mat,2,sum)>=15)]
mut_result_filt = mut_mat[,which(apply(mut_mat,2,sum)>=20)]

colnames(mut_result_filt)[which(colnames(mut_result_filt)=="U2AF1;U2AF1L5")] = "U2AF1"
mut_result_filt = mut_result_filt[,order(apply(mut_result_filt,2,sum),decreasing=F)]


ordered = c(
"cna_2pLOH",
"cna_4qLOH",
"cna_17pLOH",
"1_CNN-LOH_p",
"1_CNN-LOH_q",
"6_CNN-LOH_p",
"9_CNN-LOH_p",
"9_CNN-LOH_q",
"11_CNN-LOH_q",
"14_CNN-LOH_q",
"16_CNN-LOH_p",
"17_CNN-LOH_q",
"5_Deletion_q",
"6_Deletion_q",
"11_Deletion_q",
"13_Deletion_q",
"14_Deletion_q",
"20_Deletion_q",
"15_Duplication_q",
"21_Duplication_q",
"22_Duplication_q"
)

order_label = c(
"2pLOH","4qLOH","17pLOH",
"1pUPD","1qUPD",
"6pUPD","9pUPD",
"9qUPD",
"11qUPD",
"14qUPD","16pUPD",
"17qUPD",
"del(5q)","del(6q)",
"del(11q)","del(13q)",
"del(14q)","del(20q)",
"+15q",
"+21q","+22q"
)

cna_result_filt = cna_result_filt[,rev(ordered)]

mut_result_filt_tmp = cbind(cna_result_filt,mut_result_filt)
cna_result_filt_tmp = cbind(cna_result_filt,mut_result_filt)

mut_result_filt = mut_result_filt_tmp
cna_result_filt = cna_result_filt_tmp

rectanCor_or = rectanCor_p = matrix(NA,ncol=ncol(cna_result_filt),nrow=ncol(mut_result_filt))
colnames(rectanCor_or) = colnames(rectanCor_p) = colnames(cna_result_filt)
rownames(rectanCor_or) = rownames(rectanCor_p) = colnames(mut_result_filt)

for (i in 1:ncol(cna_result_filt)){
	for (j in 1:ncol(mut_result_filt)){
		pp = sum((cna_result_filt[,i]==1)*(mut_result_filt[,j]==1))
		pn = sum((cna_result_filt[,i]==1)*(mut_result_filt[,j]!=1))
		np = sum((cna_result_filt[,i]!=1)*(mut_result_filt[,j]==1))
		nn = sum((cna_result_filt[,i]!=1)*(mut_result_filt[,j]!=1))
		rectanCor_or[j,i] = odds_ratio(pp,pn,np,nn)["or"]
		rectanCor_p[j,i]  = odds_ratio(pp,pn,np,nn)["p"]
	}
}


res = matrix(NA,ncol=5,nrow=ncol(cna_result_filt)*ncol(mut_result_filt))
k = 1
for (i in 1:ncol(cna_result_filt)){
	for (j in 1:ncol(mut_result_filt)){
		pp = sum((cna_result_filt[,i]==1)*(mut_result_filt[,j]==1))
    pn = sum((cna_result_filt[,i]==1)*(mut_result_filt[,j]!=1))
    np = sum((cna_result_filt[,i]!=1)*(mut_result_filt[,j]==1))
    nn = sum((cna_result_filt[,i]!=1)*(mut_result_filt[,j]!=1))
		res[k,] = c(colnames(mut_result_filt)[j],colnames(cna_result_filt)[i],odds_ratio(pp,pn,np,nn)[c("or","p")],pp)
		k = k + 1
	}
}

res = res[which(res[,1]!=res[,2]),]
a = rep(NA,nrow(res))
for(i in 1:nrow(res)){b=sort(c(res[i,1],res[i,2]));a[i]=paste0(b[1],",",b[2])}
res = res[which(!duplicated(a)),]

P = as.numeric(res[,4])
Q = p.adjust(P,"BH")
res[,4] = Q
res_sel = res[which(as.numeric(res[,4])<0.1),]
mut_lab = unique(res_sel[,1])
cna_lab = unique(res_sel[,1])

x_wid = nrow(rectanCor_or)
y_wid = ncol(rectanCor_or)
wid = max(x_wid,y_wid)


pdf("rectan_cor_2.pdf",width=10,height=10)
plot(NULL, NULL,xlim=c(-3,wid+3+10),ylim=c(-3,wid+3+10),
axes=FALSE,xlab="", ylab="",main=NULL)

for(i in 1:nrow(res)){
	x = which(rev(colnames(mut_result_filt_tmp)) == res[i,1])
	y = which(rev(colnames(mut_result_filt_tmp)) == res[i,2])
	
	or = as.numeric(res[i,3])
	q = as.numeric(res[i,4])
	n_pp = as.numeric(res[i,5]) 

	size = (q >= 0.1) * 0.1 + ( q < 0.1 & q >= 0.01 ) * 0.3 + ( q < 0.01 & q >= 0.001 ) * 0.6 + ( q < 0.001 ) * 1

	color_id = round( log10(or) * 25 ) + 51
	if(color_id > 101){color_id = 101}else if(color_id < 0) {color_id = 1}
	bwr = colorRampPalette(colors=c("blue", "white", "red"))(101)
	color = bwr[color_id]  

	rect(x-0.5-size/2,(y_wid-y+1)-0.5-size/2,x-0.5+size/2,(y_wid-y+1)-0.5+size/2,col=color,border=F)
	rect(x-1,(y_wid-y+1)-1,x,(y_wid-y+1),density=0,border="gray",lwd=3) 
	
	if(n_pp >= 5){
		segments(x-0.5,(y_wid-y+1)-0.5-0.15,x-0.5,(y_wid-y+1)-0.5+0.15,lwd=0.5)
		segments(x-0.5-0.1,(y_wid-y+1)-0.5+0.1,x-0.5+0.1,(y_wid-y+1)-0.5-0.1,lwd=0.5)
		segments(x-0.5-0.1,(y_wid-y+1)-0.5-0.1,x-0.5+0.1,(y_wid-y+1)-0.5+0.1,lwd=0.5)
	}  
}

labels = c(
rev(colnames(mut_result_filt_tmp))[1:12],
"2pLOH","4qLOH","17pLOH",
"1pUPD","1qUPD","6pUPD","9pUPD","9qUPD",
"11qUPD","14qUPD","16pUPD","17qUPD",
"del(5q)","del(6q)","del(11q)","del(13q)",
"del(14q)","del(20q)","+15q","+21q","+22q"
)
font_label = c(rep(4,14),rep(2,length(labels)-12))


for (i in 1:(ncol(mut_result_filt_tmp)-1)){
	text(i-1+0.4,-0.7,cex=0.45,getCnaLabel2(rev(colnames(mut_result_filt_tmp))[i]),adj=1,srt=45,font=font_label[i])
}

for (i in 2:ncol(mut_result_filt_tmp)){
	text(-0.4,(ncol(mut_result_filt_tmp)-i+1)-0.5,cex=0.45,getCnaLabel2(rev(colnames(mut_result_filt_tmp))[i]),adj=1,srt=0,font=font_label[i])
}

dev.off()


pdf("rectan_cor_2_legend.pdf",width=10,height=10)
plot(NULL, NULL,xlim=c(-3,wid+3+10),ylim=c(-3,wid+3+10),
axes=FALSE,xlab="", ylab="",main=NULL)

legend_x = 15

# size scale (positive)
size = 0.3
x_loc = legend_x-4.5; y_loc = y_wid*(7/6)
rect(x_loc-size/2,y_loc-size/2,x_loc+size/2,y_loc+size/2,col="gray",border=F)
text(x_loc+1,y_loc,"< 0.1",cex=0.7,font=2,adj=0)
size = 0.6
y_loc = y_loc-1.5
rect(x_loc-size/2,y_loc-size/2,x_loc+size/2,y_loc+size/2,col="gray",border=F)
text(x_loc+1,y_loc,"< 0.01",cex=0.7,font=2,adj=0)
size = 0.9
y_loc = y_loc-2
rect(x_loc-size/2,y_loc-size/2,x_loc+size/2,y_loc+size/2,col="gray",border=F)
text(x_loc+1,y_loc,"< 0.001",cex=0.7,font=2,adj=0)

text(legend_x-4.5,y_loc+5,"q",cex=0.7,font=4,adj=0.5)
text(legend_x-3,y_loc+5,"value",cex=0.7,font=2,adj=0.5)


# odds ratio
text(legend_x+4.5,y_loc+5,cex=0.7,"Odds ratio",font=2,adj=0.5)
text(legend_x+5,y_loc+0,cex=0.7,"< 0.01",font=2,adj=0)
text(legend_x+5,y_loc+1.5,cex=0.7,"1",font=2,adj=0)
text(legend_x+5,y_loc+3,cex=0.7,"> 100",font=2,adj=0)

# color scale
for (i in 1:101){
	rect(legend_x+3-0.3,y_loc+(i-1)/(101/3),legend_x+3+0.3,y_loc+i/(101/3),col=rev(bwr)[i],border=F)
}

dev.off()


## End of analysis ##
q()





