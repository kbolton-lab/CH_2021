# library
library(tidyverse)

id_path = as.character(1:11234)
mut_path = "MatMed_2021_mut_table.txt"
cna_path = "MatMed_2021_cna_table.csv"  


# get info of chr length
arminfo = read.table(file="cytoBand/arm_chr_length.txt")
cent = arminfo[1:22,4]


mut = read.table(mut_path,header=T,stringsAsFactor=F,sep="\t",quote="",colClass="character")
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
	}else if(mut$Merge_Func[i] == "focal deletion"){
		mat_mut[cur_id,cur_gene] = 1
	}
}


cna = read.csv(cna_path,header=T,stringsAsFactor=F,quote="",comment.char="")
cna = cna[which(!is.na(cna$chr)),]
label_cna = paste0(sort(rep(1:22,2)),":",c("p","q"))
mat_cna =  matrix(0,ncol=length(label_cna),nrow=length(id))
colnames(mat_cna) = label_cna; rownames(mat_cna) = id

cna_count = matrix(0,ncol=4,nrow=nrow(mat_cna))
rownames(cna_count) = rownames(mat_cna)

for(i in 1:nrow(cna)){
	cur_id = cna$id[i]
	if(! is.element(cur_id,id)){
		next
	}

 	cur_type = "unknown"
	if(cna$COPY_CHANGE[i] == "loss"){
		cur_type = -2
	}else if(cna$COPY_CHANGE[i] == "neutral"){
		cur_type = -3
	}else if(cna$COPY_CHANGE[i] == "gain"){
		cur_type = -1
	}else{
		cur_type = -5
	}

	cur_cna = cna$chr[i]
	
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
	n_row = ifelse(cur_type == -5,1,ifelse(cur_type == -1,2,ifelse(cur_type == -2,3,ifelse(cur_type == -3,4,NA))))
	cna_count[cur_id,n_row] = cna_count[cur_id,n_row] + 1
	mat_cna[cur_id,cur_cna] = cur_type
}



mat_cna = mat_cna[,which(apply(mat_cna,2,sum)!=0)]

a1 = mat_mut[,"DNMT3A"]!=0 & (mat_cna[,"2:p"]==-5 | mat_cna[,"2:p"]==-2 | mat_cna[,"2:p"]==-3)
a2 = mat_mut[,"TET2"]!=0 & (mat_cna[,"4:q"]==-5 | mat_cna[,"4:q"]==-2 | mat_cna[,"4:q"]==-3)
a3 = mat_mut[,"JAK2"]!=0 & mat_cna[,"9:p"]==-3
a4 = mat_mut[,"TP53"]!=0 & (mat_cna[,"17:p"]==-5 | mat_cna[,"17:p"]==-2 | mat_cna[,"17:p"]==-3)
aa = sum(a1 | a2 | a3 | a4)

sum(a1)/sum(mat_mut[,"DNMT3A"]!=0)
sum(a1)
sum(mat_mut[,"DNMT3A"]!=0)
sum(a2)/sum(mat_mut[,"TET2"]!=0)
sum(a2)
sum(mat_mut[,"TET2"]!=0)
sum(a3)/sum(mat_mut[,"JAK2"]!=0)
sum(a3)
sum(mat_mut[,"JAK2"]!=0)
sum(a4)/sum(mat_mut[,"TP53"]!=0)
sum(a4)
sum(mat_mut[,"TP53"]!=0)

if(T){
N = length(id)
#N = N - aa
m = sum(apply(mat_mut,1,sum) != 0)
#m = m - aa
c = sum(apply(mat_cna,1,sum) != 0)
#c = c - aa
its = sum(apply(mat_mut,1,sum) != 0 & apply(mat_cna,1,sum) != 0)
#its = its - aa

print(m)
print(c)
print(m*100/10612)
print(c*100/10612)
print(its)
print(m+c-its)
print(its*100/10612)
print(its*100/(m+c-its))

print(fisher.test(matrix(c(its,m-its,c-its,N-m-c+its),ncol=2)))

print(table((apply(mut_count,1,sum)>=2)+(apply(cna_count,1,sum)>=2)))
print(table(apply(mut_count,1,sum)+apply(cna_count,1,sum)))

area = function(r1,r2,d) {
	theta1 = acos(min(max((r1^2+d^2-r2^2)/(2*r1*d),-1),1))
	theta2 = acos(min(max((r2^2+d^2-r1^2)/(2*r2*d),-1),1))
	res = r1^2*(pi-theta1+sin(theta1)*cos(theta1))+r2^2*(pi-theta2+sin(theta2)*cos(theta2))
	return(res)
}

r1 = sqrt(m/pi)
r2 = sqrt(c/pi)
a = m+c-its
d = uniroot(function(d) area(r1,r2,d)-a, c(abs(r1-r2),r1+r2))$root

pdf("Venn_mutation_cna.pdf")
plot(NULL, xlim=c(-r1,d+r2), ylim=c(-r1,r1), axes=FALSE, asp=1, xlab="", ylab="")
theta = seq(0, 2*pi, length.out=100)
polygon(r1*cos(theta),r1*sin(theta),col=rgb(1,0,0,alpha=0.65),border = "transparent")
polygon(r2*cos(theta)+d,r2*sin(theta),col=rgb(0,0,1,alpha=0.65),border = "transparent")

text((d-r2-r1)/2,(1/8)*r1,"Mutation",font=2,cex=1.8,col="white")
text((d-r2-r1)/2,-(1/8)*r1,paste0(round(100*(m-its)/N)," %"),font=2,cex=1.8,col="white")

text((r1+d+r2)/2,(1/8)*r1,"CNA",font=2,cex=1.8,col="white")
text((r1+d+r2)/2,-(1/8)*r1,paste0(round(100*(c-its)/N)," %"),font=2,cex=1.8,col="white")

text((r1+d-r2)/2,(1/8)*r1,"Both",font=2,cex=1.8,col="white")
text((r1+d-r2)/2,-(1/8)*r1,paste0(round(100*(its)/N)," %"),font=2,cex=1.8,col="white")
dev.off()
}

mut_label = names(sort(apply(mat_mut!=0,2,sum),decreasing=T))
#mut_label = setdiff(names(sort(apply(mat_mut!=0,2,sum),decreasing=T)),c("RUNX1","ETV6","GATA2","CEBPA","DDX41"))
#c(
#"DNMT3A","TET2","IDH2","IDH1",
#"ASXL1","EZH2",
#"SF3B1","SRSF2","U2AF1;U2AF1L5",
#"PPM1D","TP53",
#"CBL","NRAS","KRAS",
#"JAK2","MYD88",
#"GNB1","GNAS"
#,"RUNX1","ETV6","GATA2","CEBPA"
#,"DDX41"
#)

cna_label = names(sort(apply(mat_cna!=0,2,sum),decreasing=T)) #[1:20]
# cna_label = c(
# "21:Duplication:q","20:Deletion:q","11:CNN-LOH:p","11:CNN-LOH:q",
# "9:CNN-LOH:p","12:Duplication:q","14:CNN-LOH:q",
# "1:CNN-LOH:p","1:CNN-LOH:q","13:Deletion:q",
# "22:CNN-LOH:q","8:Duplication:q","13:CNN-LOH:q",
# "15:Duplication:q","17:CNN-LOH:p","17:CNN-LOH:q",
# "4:CNN-LOH:q","7:CNN-LOH:q","15:CNN-LOH:q","20:Duplication:q"
# )

mat_mut = mat_mut[,mut_label]
mat_cna = mat_cna[,cna_label]

cur_order = id

mat_cna_2 = mat_cna
mat_mut_2 = mat_mut

for(i in ncol(mat_cna_2):1){
	cur_order = cur_order[order(1 / mat_cna_2[,i],decreasing=F)]
	mat_cna_2 = mat_cna_2[cur_order,]
	mat_mut_2 = mat_mut_2[cur_order,]
}

for(i in ncol(mat_mut_2):1){
	cur_order = cur_order[order(mat_mut_2[,i],decreasing=T)]
	mat_cna_2 = mat_cna_2[cur_order,]
	mat_mut_2 = mat_mut_2[cur_order,]
}

num_cutoff_cna = 50
mat_cna_3 = mat_cna[cur_order,which(apply(mat_cna[cur_order,]!=0,2,sum) >= num_cutoff_cna)]
other_cna = as.numeric(apply(mat_cna[cur_order,which(apply(mat_cna[cur_order,]!=0,2,sum) < num_cutoff_cna)],1,sum)!=0) 

mat = cbind(mat_mut[cur_order,],mat_cna_3,other_cna)
mat = mat[which(apply(mat!=0,1,sum)!=0),]

mut_count = mut_count[cur_order,]
cna_count = cna_count[cur_order,]


label_y = c(mut_label,gsub(":","",colnames(mat_cna_2)),"others")
label_y[which(label_y == "U2AF1;U2AF1L5")] = "U2AF1"

col1 = "gray80"
col2 = "gray85"
label_col = rep(c(col1,col2),100)
#label_col = c(
#rep(col1,18),
#rep(col2,length(cna_label))
#)

label_font = c(
rep(4,length(mut_label)),
rep(2,length(label_y)-length(mut_label))
)

prop = 10
d = 1.1
x_n = nrow(mat); x_wid = nrow(mat)/prop
y_n = ncol(mat); y_wid = ncol(mat)


tag = "riken_landscape_unsorted.pdf"
pdf(tag,width=8,height=6)

plot(NULL, NULL,
xlim=c(-15,x_wid+70),ylim=c(-y_wid-1,5),
axes=FALSE,
xlab="", ylab="",
main=NULL
)

#### editing ###

get_col = function(n){
	if(n==-1){ # 1: Duplication
		col = "firebrick3"
	}else if(n==-2){ # 2: Deletion
		col = "dodgerblue3"
	}else if(n==-3){ # 3: CNN-LOH
		col = "chartreuse3"
	}else if(n==-5){
		col = "tan"
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
for (i in 1:x_n){
	# mutations
	for (j in 1:ncol(mat_mut)){
		cur_var = mat[i,j]
		col = get_col(cur_var)
		
		if(cur_var!=0){
			rect((i-1)/prop,(-j+1)/d,i/prop,-j/d,col=col,border=F)
		}else{
			rect((i-1)/prop,(-j+1)/d,i/prop,-j/d,col=label_col[j],border=F)
		}
#		rect((i-1)/prop,-j+1,i/prop,-j,density=0,border="white",lwd = 0.00001)
	}
	
	# CNAs
	for (j in (ncol(mat_mut)+1):y_n){
		cur_var = mat[i,j]
		if(j != y_n){
			col = get_col(cur_var)
		}else{
			cur_type = "other"
			n_col = -4
			col = get_col(n_col)
		}
		
		if(cur_var!=0){
			rect((i-1)/prop,(-j+1)/d,i/prop,-j/d,col=col,border=F)
			#mat[i,j] = n_col
		}else{
			rect((i-1)/prop,(-j+1)/d,i/prop,-j/d,col=label_col[j],border=F)
		}
#		rect((i-1)/prop,-j+1,i/prop,-j,density=0,border="white",lwd=0.00001)
	}
}

## name of alterations ##
for (i in 1:y_n){
	text(-2.5,(-i+0.5)/d,cex=0.4,label_y[i],adj=1,srt=0,font=label_font[i])
	#segments(0,(-i)/d,x_n/prop,(-i)/d,col="white",lwd=0.01)
}

## upper side, number of alterations ##
gp = 1.0
vat = 1.5
types = c(1:7,-5,(-1):(-3))
for(i in 1:x_n){
	counts = c(mut_count[i,],cna_count[i,])
	c_sum = c(0,cumsum(counts))
	for(k in 1:(length(c_sum)-1)){
		rect((i-1)/prop,c_sum[k]/(vat*d)+gp,i/prop,c_sum[k+1]/(vat*d)+gp,col=get_col(types[k]),border=F)
	}
}

text(-7,3.5/(vat*d)+gp,"Count",adj=1,cex=0.4,font=2)
segments(-2,0/(vat*d)+gp,-2,8/(vat*d)+gp,lwd=0.3)
for(i in 0:8){
	segments(-3,i/(vat*d)+gp,-2,i/(vat*d)+gp,lwd=0.3)
}

## right side, frequency graph ##
sp = 5
gp = 0.1
types = c(8:1,(-1):(-4),-5)
for(i in 1:(y_n-1)){
	counts = rep(NA,12)
	for(j in 1:length(types)){
		counts[j] = sum(mat[,i] == types[j])
	}
	pct = c(0,cumsum(counts/length(id)))
	for(k in 1:(length(pct)-1)){
		rect((x_n+pct[k]*5000)/prop+sp,(-i+1-gp)/d,(x_n+pct[k+1]*5000)/prop+sp,(-i+gp)/d,col=get_col(types[k]),border=F)
	}
}

text((x_n+0.07*5000)/prop+sp,4/d,"Count",cex=0.4,font=2)
segments((x_n+0*5000)/prop+sp,0.5/d,(x_n+0.15*5000)/prop+sp,0.5/d,lwd=0.3)
lat = c(0,500,1000,1500)
lat2 = 0.01*c(1:4,6:9,11:14,16)
for(i in 1:length(lat)){
	x = (x_n+(lat[i]/length(id))*5000)/prop+sp
	segments(x,0.5/d,x,0.80/d,lwd=0.3)
}

x = (x_n+(lat[1]/length(id))*5000)/prop+sp
#text(x,2/d,paste0((lat[1]*100),"%"),adj=0,cex=0.4,font=2)
x = (x_n+(lat[2]/length(id))*5000)/prop+sp
#text(x,2/d,paste0((lat[2]*100),"%"),adj=0.5,cex=0.4,font=2)
x = (x_n+(lat[3]/length(id))*5000)/prop+sp
#text(x,2/d,paste0((lat[3]*100),"%"),adj=0.5,cex=0.4,font=2)
x = (x_n+(lat[4]*5000))/prop+sp
#text(x,2/d,paste0((lat[4]*100),"%"),adj=0.5,cex=0.4,font=2)

for(i in 1:length(lat2)){
	x = (x_n+(lat2[i]/length(id))*5000)/prop+sp
	segments(x,0.5/d,x,0.65/d,lwd=0.3)
}

dev.off()

# legend
pdf("landscape/legend.pdf")

plot(NULL, NULL,
xlim=c(0,10),ylim=c(0,15),
axes=FALSE,
xlab="", ylab="",
main=NULL
)   

# mutation
text(1,12.7,"Mutation",font=2,adj=0)
rect(1,11,1+1/4,12,col=get_col(7),border=F); rect(1,11,1+1/4,12,density=0,border="white",lwd=0.01); text(1+1/2,11.5,"Missense",adj=0)
rect(1,10,1+1/4,11,col=get_col(6),border=F); rect(1,10,1+1/4,11,density=0,border="white",lwd=0.01); text(1+1/2,10.5,"Inframe indel",adj=0)
rect(1,9,1+1/4,10,col=get_col(5),border=F); rect(1,9,1+1/4,10,density=0,border="white",lwd=0.01); text(1+1/2,9.5,"Splice-site",adj=0)
rect(1,8,1+1/4,9,col=get_col(4),border=F); rect(1,8,1+1/4,9,density=0,border="white",lwd=0.01); text(1+1/2,8.5,"Frameshift indel",adj=0)
rect(1,7,1+1/4,8,col=get_col(3),border=F); rect(1,7,1+1/4,8,density=0,border="white",lwd=0.01); text(1+1/2,7.5,"Stop-gain",adj=0)
rect(1,6,1+1/4,7,col=get_col(2),border=F); rect(1,6,1+1/4,7,density=0,border="white",lwd=0.01); text(1+1/2,6.5,"Multiple",adj=0)

# cna
text(1,4.7,"CNA",font=2,adj=0)
rect(1,3,1+1/4,4,col=get_col(-1),border=F); rect(1,3,1+1/4,4,density=0,border="white",lwd=0.01); text(1+1/2,3.5,"Duplication",adj=0)
rect(1,2,1+1/4,3,col=get_col(-2),border=F); rect(1,2,1+1/4,3,density=0,border="white",lwd=0.01); text(1+1/2,2.5,"Deletion",adj=0)
rect(1,1,1+1/4,2,col=get_col(-3),border=F); rect(1,1,1+1/4,2,density=0,border="white",lwd=0.01); text(1+1/2,1.5,"UPD",adj=0)
rect(1,0,1+1/4,1,col=get_col(-5),border=F); rect(1,0,1+1/4,1,density=0,border="white",lwd=0.01); text(1+1/2,0.5,"Unclassifiable",adj=0)

dev.off()




