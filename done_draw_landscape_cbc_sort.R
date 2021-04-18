
## ED Fig.7a ##

## import genome data ##
mut_path = "MatMed_2021_mut_table.txt"
cna_path = "MatMed_2021_cna_table.csv"
id = 1:11234

## import age information ##
age = ~~

## get info of chr length ##
arminfo = read.table(file="cytoBand/arm_chr_length.txt")
cent = arminfo[1:22,4]

mut = read.table(mut_path,header=T,stringsAsFactor=F,quote="",sep="\t",comment.char="")
mut = mut[which(!is.na(mut$Chr)),]

mat_mut = matrix(0,ncol=length(unique(mut$Gene.refGene)),nrow=length(id))
colnames(mat_mut) = unique(mut$Gene.refGene); rownames(mat_mut) = id

mut_count = matrix(0,ncol=7,nrow=length(id))
rownames(mut_count) = id

for(i in 1:nrow(mut)){
	cur_id = mut$id[i]
	cur_gene = mut$Gene.refGene[i]
	if(! is.element(cur_id,id)){
		next
	}
	
	if(mut$Merge_Func[i] == "stopgain"){
		mut_count[cur_id,3] = mut_count[cur_id,3] + 1
	}else if(mut$Merge_Func[i] == "frameshift deletion" | mut$Merge_Func[i] == "frameshift insertion"){
		mut_count[cur_id,4] = mut_count[cur_id,4] + 1
	}else if(mut$Merge_Func[i] == "splicing"){
		mut_count[cur_id,5] = mut_count[cur_id,5] + 1
	}else if(mut$Merge_Func[i] == "nonsynonymous SNV" | mut$Merge_Func[i] == "stoploss"){
		mut_count[cur_id,7] = mut_count[cur_id,7] + 1
	}else if(mut$Merge_Func[i] == "nonframeshift deletion" | mut$Merge_Func[i] == "nonframeshift insertion"){
		mut_count[cur_id,6] = mut_count[cur_id,6] + 1
	}
	
	if(mat_mut[cur_id,cur_gene] != 0){
		mat_mut[cur_id,cur_gene] = 2
	}else if(mut$Merge_Func[i] == "stopgain"){
		mat_mut[cur_id,cur_gene] = 3
	}else if(mut$Merge_Func[i] == "frameshift deletion" | mut$Merge_Func[i] == "frameshift insertion"){
		mat_mut[cur_id,cur_gene] = 4
	}else if(mut$Merge_Func[i] == "splicing"){
		mat_mut[cur_id,cur_gene] = 5
	}else if(mut$Merge_Func[i] == "nonsynonymous SNV" | mut$Merge_Func[i] == "stoploss"){
		mat_mut[cur_id,cur_gene] = 7
	}else if(mut$Merge_Func[i] == "nonframeshift deletion" | mut$Merge_Func[i] == "nonframeshift insertion"){
		mat_mut[cur_id,cur_gene] = 6
	}else if(mut$Merge_Func[i]== "focal deletion"){
		mat_mut[cur_id,cur_gene] = 1
	}
}

cna = read.csv(cna_path,header=T,stringsAsFactor=F,quote="",colClass="character")
cna = cna[which(!is.na(cna$chr)),]
label_cna = paste0(sort(rep(1:22,4*2)),":",sort(rep(c("CNN-LOH","Duplication","Deletion","Unknown"),2)),":",c("p","q"))
mat_cna =  matrix(0,ncol=length(label_cna),nrow=length(id))
colnames(mat_cna) = label_cna; rownames(mat_cna) = id

for(i in 1:nrow(cna)){
	cur_id = cna$id[i]
	if(! is.element(cur_id,id)){
		next
	}

	cur_type = "Unknown"
	if(cna$COPY_CHANGE[i] == "loss"){
		cur_type = "Deletion"
	}else if(cna$COPY_CHANGE[i] == "neutral"){
		cur_type = "CNN-LOH"
	}else if(cna$COPY_CHANGE[i] == "gain"){
		cur_type = "Duplication"
	}
	
	cur_cna = paste0(cna$chr[i],":",cur_type)
	cur_cent = cent[as.numeric(cna$chr[i])]
	if(as.numeric(cna$start[i])*(10^6) <= cur_cent & as.numeric(cna$end[i])*(10^6) <= cur_cent){
		cur_cna = paste0(cur_cna,":p")
	}else if(as.numeric(cna$start[i])*(10^6) >= cur_cent & as.numeric(cna$end[i])*(10^6) >= cur_cent){
		cur_cna = paste0(cur_cna,":q")
	}else if(as.numeric(cna$start[i])*(10^6) <= cur_cent & as.numeric(cna$end[i])*(10^6) >= cur_cent){
		cur_cna = c(paste0(cur_cna,":p"),paste0(cur_cna,":q"))
	}
	if(! is.element(cur_cna,label_cna)){
		next
	}
	mat_cna[cur_id,cur_cna] = mat_cna[cur_id,cur_cna] + 1
}

cna_count = matrix(NA,ncol=4,nrow=nrow(mat_cna))
rownames(cna_count) = rownames(mat_cna)
upd = sort(c((1:(ncol(mat_cna)/8))*8-6,(1:(ncol(mat_cna)/8))*8-7))
cna_count[,4] = apply(mat_cna[,upd],1,sum) # UPD
del = sort(c((1:(ncol(mat_cna)/8))*8-4,(1:(ncol(mat_cna)/8))*8-5))
cna_count[,3] = apply(mat_cna[,del],1,sum) # Deletion
dup = sort(c((1:(ncol(mat_cna)/8))*8-2,(1:(ncol(mat_cna)/8))*8-3))
cna_count[,2] = apply(mat_cna[,dup],1,sum) # Duplication
unk = sort(c((1:(ncol(mat_cna)/8))*8-0,(1:(ncol(mat_cna)/8))*8-1))
cna_count[,1] = apply(mat_cna[,unk],1,sum) # Duplication 

mat_cna = mat_cna[,which(apply(mat_cna,2,sum)!=0)]

mut_label = names(sort(apply(mat_mut,2,sum),decreasing=T)) 
cna_label = names(sort(apply(mat_cna,2,sum),decreasing=T))

mat_mut = mat_mut[,mut_label]
mat_cna = mat_cna[,cna_label]

cur_order = id

mat_cna_2 = mat_cna
mat_mut_2 = mat_mut

for(i in ncol(mat_cna_2):1){
	cur_order = cur_order[order(mat_cna_2[,i],decreasing=T)]
	mat_cna_2 = mat_cna_2[cur_order,]
	mat_mut_2 = mat_mut_2[cur_order,]
}

for(i in ncol(mat_mut_2):1){
	cur_order = cur_order[order(mat_mut_2[,i],decreasing=T)]
	mat_cna_2 = mat_cna_2[cur_order,]
	mat_mut_2 = mat_mut_2[cur_order,]
}


## import CBC information ##
id_normal = ~~ ## id of cbc normal cases
id_abnormal = ~~ ## id of cbc abnormal cases
id_unknown = setdiff(id,c(id_normal,id_abnormal))

order_normal = cur_order[which(is.element(cur_order,id_normal))]
order_abnormal = cur_order[which(is.element(cur_order,id_abnormal))]
order_unknown = cur_order[which(is.element(cur_order,id_unknown))]

cur_order = c(order_normal,order_abnormal,order_unknown)



num_cutoff_cna = 20
mat_cna_3 = mat_cna[cur_order,which(apply(mat_cna[cur_order,],2,sum) >= num_cutoff_cna)]
other_cna = as.numeric(apply(mat_cna[cur_order,which(apply(mat_cna[cur_order,],2,sum) < num_cutoff_cna)],1,sum)!=0) 

mat = cbind(mat_mut[cur_order,],mat_cna_3,other_cna)
mat = mat[which(apply(mat!=0,1,sum)!=0),]

a_sel = length(id_normal)
a_exc = length(id_abnormal)
a_unk = length(id_unknown)

n_sel = sum(apply(mat[intersect(rownames(mat),order_normal),],1,sum)!=0)
n_exc = sum(apply(mat[intersect(rownames(mat),order_abnormal),],1,sum)!=0)
n_unk = sum(apply(mat[intersect(rownames(mat),order_unknown),],1,sum)!=0)

mut_count = mut_count[cur_order,]
cna_count = cna_count[cur_order,]

all_count = cbind(mut_count,cna_count)
all_remain = which(apply(all_count!=0,1,sum)!=0)

mut_count = mut_count[all_remain,]
cna_count = cna_count[all_remain,]

getCnaLabel = function(s){
	splited = unlist(strsplit(s,":"))
	type=NA
	prf = ifelse(splited[2]=="CNN-LOH","",ifelse(splited[2]=="Duplication","+",ifelse(splited[2]=="Deletion","del(","")))
	suf = ifelse(splited[2]=="CNN-LOH","UPD",ifelse(splited[2]=="Duplication","",ifelse(splited[2]=="Deletion",")","AI")))
	return(paste0(prf,splited[1],splited[3],suf))
}


label_y = c(
mut_label,
sapply(cna_label,getCnaLabel),
"others"
)
label_y[which(label_y == "U2AF1;U2AF1L5")] = "U2AF1"

col1 = "gray80"
col2 = "gray85"
label_col = rep(c(col1,col2),100)


label_font = c(
rep(4,length(mut_label)),
rep(2,length(label_y)-length(mut_label))
)

prop = 10
d = 1.1
x_n = nrow(mat); x_wid = nrow(mat)/prop
y_n = ncol(mat); y_wid = ncol(mat)

tag = "landscape_cbc_sorted.pdf"
pdf(tag)
plot(NULL, NULL,xlim=c(-15,x_wid+80+200),ylim=c(-y_wid-1,7),
axes=FALSE,xlab="", ylab="",main=NULL)


#### editing ###
get_col = function(n){
	if(n==-1){ # 1: Duplication
		col = "firebrick3"
	}else if(n==-2){ # 2: Deletion
		col = "dodgerblue3"
	 }else if(n==-5){
		col = "tan"
	}else if(n==-3){ # 3: CNN-LOH
		col = "chartreuse3"
	}else if(n==7){ # 4: missense
		col = "deepskyblue1"
	}else if(n==5){ # 5: splice site
		col = "purple"
	}else if(n==4){ # 6: frameshift indel
		col = "salmon"
	}else if(n==3){ # 7: stopgain
		col = "red"
	}else if(n==6){ # 8: inframe indel
		col = "yellow"
	}else if(n==2){
		col = "brown"
	}else if(n==1){
		col = "dodgerblue3"
	}else if(n==-4){ # -4: other CNAs
		col = "cornsilk4"
	}else{ # 0: none
		col = "gray"
	}	
	return(col)
}


## main, landscape ##
sp = 100
for (i in 1:x_n){
	if(i > n_sel + n_exc){
		l_sp = 2 * sp
	}else if(i > n_sel){
		l_sp = sp
	}else{
		l_sp = 0
	}
	
	# mutations
	for (j in 1:ncol(mat_mut)){
		cur_var = mat[i,j]
		col = get_col(cur_var)
		
		if(cur_var!=0){
			rect((i-1)/prop + l_sp,(-j+1)/d,i/prop + l_sp,-j/d,col=col,border=F)
		}else{
			rect((i-1)/prop + l_sp,(-j+1)/d,i/prop + l_sp,-j/d,col=label_col[j],border=F)
		}
	}
	
	# CNAs
	for (j in (ncol(mat_mut)+1):y_n){
		cur_var = mat[i,j]
		if(j != y_n){
			cur_type = unlist(strsplit(colnames(mat)[j],":"))[2]
			#n_col = ifelse(cur_type=="Duplication",-1,ifelse(cur_type=="Deletion",-2,if_else(cur_type=="CNN-LOH",-3,0)))
    	n_col = ifelse(cur_type=="Duplication",-1,ifelse(cur_type=="Deletion",-2,if_else(cur_type=="CNN-LOH",-3,ifelse(cur_type=="Unknown",-5,0))))
			col = get_col(n_col)
		}else{
			cur_type = "other"
			n_col = -4
			col = get_col(n_col)
		}
		
		if(cur_var!=0){
			rect((i-1)/prop + l_sp,(-j+1)/d,i/prop + l_sp,-j/d,col=col,border=F)
			mat[i,j] = n_col
		}else{
			rect((i-1)/prop + l_sp,(-j+1)/d,i/prop + l_sp,-j/d,col=label_col[j],border=F)
		}
	}
}

label_y[y_n] = "Others"

## name of alterations ##
for (i in 1:y_n){
	text(-2.5,(-i+0.5)/d,cex=0.35,label_y[i],adj=1,srt=0,font=label_font[i])
	#segments(0,(-i)/d,x_n/prop,(-i)/d,col="white",lwd=0.01)
}

## upper side, age ##
gp = 0.35
for (i in 1:x_n){
	if(i > n_sel + n_exc){
		l_sp = 2 * sp
	}else if(i > n_sel){
		l_sp = sp
	}else{
		l_sp = 0
	}

	cur_age = age[rownames(mat_mut)[i]]
	cur_age_2 = max(100 - round((cur_age-60)*(10/4)),0)
	#print(cur_age_2)
	cur_grey = paste0("grey",cur_age_2)
	
	rect((i-1)/prop + l_sp,1/d + gp - 0.1/d,i/prop + l_sp, gp + 0.1/d,col=cur_grey,border=F)
}
text(-2.5,0.55/d+gp,"Age",adj=1,cex=0.4,font=2)


## upper side, CBC category ##
gp = 1.62
for(i in 1:x_n){
	if(i > n_sel + n_exc){
		l_col = "burlywood3"
		l_sp = 2 * sp
	}else if(i > n_sel){
		l_col = "coral1"
		l_sp = sp
	}else{
		l_col = "cadetblue4"
		l_sp = 0
	}	
	rect((i-1)/prop + l_sp,1/d + gp - 0.1/d,i/prop + l_sp, gp + 0.1/d,col=l_col,border=F)
}
text(-2.5,0.55/d+gp,"Category",adj=1,cex=0.4,font=2)


## upper side, number of alterations ##
gp = 3
vat = 1.7
types = c(1:7,-5,(-1):(-3))
for(i in 1:x_n){
	if(i > n_sel + n_exc){
		l_sp = 2 * sp
	}else if(i > n_sel){
		l_sp = sp
	}else{
		l_sp = 0
	}

	counts = c(mut_count[i,],cna_count[i,])
	c_sum = c(0,cumsum(counts))
	for(k in 1:(length(c_sum)-1)){
		rect((i-1)/prop + l_sp,c_sum[k]/(vat*d)+gp,i/prop + l_sp,c_sum[k+1]/(vat*d)+gp,col=get_col(types[k]),border=F)
	}
}

text(-7,3.5/(vat*d)+gp,"Count",adj=1,cex=0.4,font=2)
segments(-2,0/(vat*d)+gp,-2,8/(vat*d)+gp,lwd=1)
for(i in 0:8){segments(-3,i/(vat*d)+gp,-2,i/(vat*d)+gp,lwd=1)}
segments(-2+n_sel/prop+sp,0/(vat*d)+gp,-2+n_sel/prop+sp,8/(vat*d)+gp,lwd=1) 
for(i in 0:8){segments(-3+n_sel/prop+sp,i/(vat*d)+gp,-2+n_sel/prop+sp,i/(vat*d)+gp,lwd=1)}
segments(-2+(n_sel+n_exc)/prop+2*sp,0/(vat*d)+gp,-2+(n_sel+n_exc)/prop+2*sp,8/(vat*d)+gp,lwd=1) 
for(i in 0:8){segments(-3+(n_sel+n_exc)/prop+2*sp,i/(vat*d)+gp,-2+(n_sel+n_exc)/prop+2*sp,i/(vat*d)+gp,lwd=1)}


## right side, frequency graph ##
sp = 100
gp = 0.1
types = c(8:1,(-1):(-4),-5)
wid_frq = 4000
for(i in 1:(y_n-1)){
	counts = rep(NA,12)
	for(j in 1:length(types)){
		counts[j] = sum(mat[1:n_sel,i] == types[j])
	}
	pct = c(0,cumsum(counts/a_sel))
	for(k in 1:(length(pct)-1)){
		#rect((x_n+pct[k]*3000)/prop+3*sp,(-i+1-gp)/d,(x_n+pct[k+1]*3000)/prop+3*sp,(-i+gp)/d,col=get_col(types[k]),border=F)
		rect((n_sel+pct[k]*wid_frq)/prop+0*sp+5,(-i+1-gp)/d,(n_sel+pct[k+1]*wid_frq)/prop+0*sp+5,(-i+gp)/d,col=get_col(types[k]),border=F)
	}
}

text((n_sel+0.07*wid_frq)/prop+0*sp+5+5,4/d,"Frequency (%)",cex=0.4,font=2)
segments((n_sel+0*wid_frq)/prop+0*sp+5,0.5/d,(n_sel+0.15*wid_frq)/prop+0*sp+5,0.5/d,lwd=1)
lat = c(0,0.05,0.1,0.15)
lat2 = 0.01*c(1:4,6:9,11:14)
for(i in 1:length(lat)){
	x = (n_sel+lat[i]*wid_frq)/prop+0*sp+5
	segments(x,0.5/d,x,0.80/d,lwd=1)
}

x = (n_sel+lat[1]*wid_frq)/prop+0*sp+5
text(x,2/d,paste0((lat[1]*100)),adj=0,cex=0.4,font=2)
x = (n_sel+lat[2]*wid_frq)/prop+0*sp+5
text(x,2/d,paste0((lat[2]*100)),adj=0.5,cex=0.4,font=2)
x = (n_sel+lat[3]*wid_frq)/prop+0*sp+5
text(x,2/d,paste0((lat[3]*100)),adj=0.5,cex=0.4,font=2)
x = (n_sel+lat[4]*wid_frq)/prop+0*sp+5
text(x,2/d,paste0((lat[4]*100)),adj=0.5,cex=0.4,font=2)

for(i in 1:length(lat2)){
	x = (n_sel+lat2[i]*wid_frq)/prop+0*sp+5
	segments(x,0.5/d,x,0.65/d,lwd=1)
}

## right side, frequency graph ##
sp = 100
gp = 0.1
types = c(8:1,(-1):(-4),-5)
wid_frq = 4000
for(i in 1:(y_n-1)){
	counts = rep(NA,12)
	for(j in 1:length(types)){
		counts[j] = sum(mat[(n_sel+1):(n_sel+n_exc),i] == types[j])
	}
	pct = c(0,cumsum(counts/a_exc))
	for(k in 1:(length(pct)-1)){
		#rect((x_n+pct[k]*3000)/prop+3*sp,(-i+1-gp)/d,(x_n+pct[k+1]*3000)/prop+3*sp,(-i+gp)/d,col=get_col(types[k]),border=F)
		rect((n_sel+n_exc+pct[k]*wid_frq)/prop+1*sp+5,(-i+1-gp)/d,(n_sel+n_exc+pct[k+1]*wid_frq)/prop+1*sp+5,(-i+gp)/d,col=get_col(types[k]),border=F)
	}
}

text((n_sel+n_exc+0.07*wid_frq)/prop+1*sp+5+5,4/d,"Frequency (%)",cex=0.4,font=2)
segments((n_sel+n_exc+0*wid_frq)/prop+1*sp+5,0.5/d,(n_sel+n_exc+0.15*wid_frq)/prop+1*sp+5,0.5/d,lwd=1)
lat = c(0,0.05,0.1,0.15)
lat2 = 0.01*c(1:4,6:9,11:14)
for(i in 1:length(lat)){
	x = (n_sel+n_exc+lat[i]*wid_frq)/prop+1*sp+5
	segments(x,0.5/d,x,0.80/d,lwd=1)
}

x = (n_sel+n_exc+lat[1]*wid_frq)/prop+1*sp+5
text(x,2/d,paste0((lat[1]*100)),adj=0,cex=0.4,font=2)
x = (n_sel+n_exc+lat[2]*wid_frq)/prop+1*sp+5
text(x,2/d,paste0((lat[2]*100)),adj=0.5,cex=0.4,font=2)
x = (n_sel+n_exc+lat[3]*wid_frq)/prop+1*sp+5
text(x,2/d,paste0((lat[3]*100)),adj=0.5,cex=0.4,font=2)
x = (n_sel+n_exc+lat[4]*wid_frq)/prop+1*sp+5
text(x,2/d,paste0((lat[4]*100)),adj=0.5,cex=0.4,font=2)

for(i in 1:length(lat2)){
	x = (n_sel+n_exc+lat2[i]*wid_frq)/prop+1*sp+5
	segments(x,0.5/d,x,0.65/d,lwd=1)
}

## right side, frequency graph ##
sp = 100
gp = 0.1
types = c(8:1,(-1):(-4),-5)
wid_frq = 4000
for(i in 1:(y_n-1)){
	counts = rep(NA,12)
	for(j in 1:length(types)){
		counts[j] = sum(mat[(n_sel+n_exc+1):nrow(mat),i] == types[j])
	}
	pct = c(0,cumsum(counts/a_unk))
	for(k in 1:(length(pct)-1)){
		#rect((x_n+pct[k]*3000)/prop+3*sp,(-i+1-gp)/d,(x_n+pct[k+1]*3000)/prop+3*sp,(-i+gp)/d,col=get_col(types[k]),border=F)
		rect((x_n+pct[k]*wid_frq)/prop+2*sp+5,(-i+1-gp)/d,(x_n+pct[k+1]*wid_frq)/prop+2*sp+5,(-i+gp)/d,col=get_col(types[k]),border=F)
	}
}

text((x_n+0.07*wid_frq)/prop+2*sp+5+5,4/d,"Frequency (%)",cex=0.4,font=2)
segments((x_n+0*wid_frq)/prop+2*sp+5,0.5/d,(x_n+0.15*wid_frq)/prop+2*sp+5,0.5/d,lwd=1)
lat = c(0,0.05,0.1,0.15)
lat2 = 0.01*c(1:4,6:9,11:14)
for(i in 1:length(lat)){
	x = (x_n+lat[i]*wid_frq)/prop+2*sp+5
	segments(x,0.5/d,x,0.80/d,lwd=1)
}

x = (x_n+lat[1]*wid_frq)/prop+2*sp+5
text(x,2/d,paste0((lat[1]*100)),adj=0,cex=0.4,font=2)
x = (x_n+lat[2]*wid_frq)/prop+2*sp+5
text(x,2/d,paste0((lat[2]*100)),adj=0.5,cex=0.4,font=2)
x = (x_n+lat[3]*wid_frq)/prop+2*sp+5
text(x,2/d,paste0((lat[3]*100)),adj=0.5,cex=0.4,font=2)
x = (x_n+lat[4]*wid_frq)/prop+2*sp+5
text(x,2/d,paste0((lat[4]*100)),adj=0.5,cex=0.4,font=2)

for(i in 1:length(lat2)){
	x = (x_n+lat2[i]*wid_frq)/prop+2*sp+5
	segments(x,0.5/d,x,0.65/d,lwd=1)
}

dev.off()


# legend
pdf("landscape_cbc_sorted_legend.pdf")

plot(NULL, NULL,
xlim=c(0,15),ylim=c(0,22),
axes=FALSE,
xlab="", ylab="",
main=NULL
)   

# category
text(1,21.7,"CBC category",font=2,adj=0)
rect(1,20,1+1/4,21,col="cadetblue4",border=F); text(1+1/2,20.5,"Normal",adj=0,font=2,cex=0.8)
rect(1,19,1+1/4,20,col="coral1",border=F); text(1+1/2,19.5,"Abnormal",adj=0,font=2,cex=0.8)
rect(1,18,1+1/4,19,col="burlywood3",border=F); text(1+1/2,18.5,"Unknown",adj=0,font=2,cex=0.8)

# age
text(1,16.7,"Age",font=2,adj=0)
for(i in 1:40){
	cur_age_2 = 100 - round(i*(10/4))
	cur_grey = paste0("grey",cur_age_2)
	rect(1,i*0.05+14,1+1/4,(i+1)*0.05+14,col=cur_grey,border=F)
}
text(1+1/2,16-0.25,"100",adj=0,font=2,cex=0.8)
text(1+1/2,14+0.25,"60",adj=0,font=2,cex=0.8)

# mutation
text(1,12.7,"Mutation",font=2,adj=0)
rect(1,11,1+1/4,12,col=get_col(7),border=F); text(1+1/2,11.5,"Missense",adj=0,font=2,cex=0.8)
rect(1,10,1+1/4,11,col=get_col(6),border=F); text(1+1/2,10.5,"Inframe indel",adj=0,font=2,cex=0.8)
rect(1,9,1+1/4,10,col=get_col(5),border=F); text(1+1/2,9.5,"Splice-site",adj=0,font=2,cex=0.8)
rect(1,8,1+1/4,9,col=get_col(4),border=F); text(1+1/2,8.5,"Frameshift indel",adj=0,font=2,cex=0.8)
rect(1,7,1+1/4,8,col=get_col(3),border=F); text(1+1/2,7.5,"Stop-gain",adj=0,font=2,cex=0.8)
rect(1,6,1+1/4,7,col=get_col(2),border=F); text(1+1/2,6.5,"Multiple",adj=0,font=2,cex=0.8)

# cna
text(1,4.7,"CNA",font=2,adj=0)
rect(1,3,1+1/4,4,col=get_col(-3),border=F); text(1+1/2,3.5,"UPD",adj=0,font=2,cex=0.8)
rect(1,2,1+1/4,3,col=get_col(-2),border=F); text(1+1/2,2.5,"Deletion",adj=0,font=2,cex=0.8)
rect(1,1,1+1/4,2,col=get_col(-1),border=F); text(1+1/2,1.5,"Duplication",adj=0,font=2,cex=0.8)
rect(1,0,1+1/4,1,col=get_col(-5),border=F); text(1+1/2,0.5,"Allele imbalance (not classified)",adj=0,font=2,cex=0.8)

dev.off()


