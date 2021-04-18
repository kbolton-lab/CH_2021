library(tidyverse)

chr = as.numeric(commandArgs(trailingOnly = T)[1])
wid = 0.2

mut_path = "MatMed_2021_mut_table.txt"
cna_path = "MatMed_2021_cna_table.csv"
id = as.character(1:11234)

cna = read.csv(cna_path,header=T,stringsAsFactor=F,quote="",comment.char="")
cna = cna[which(!is.na(cna$chr)),]
cna

cytoinfo <- read.table(file="cytoBand/cytoBand.txt")
arminfo <- read.table(file="cytoBand/arm_chr_length.txt")
cytoinfo$col <- ifelse(cytoinfo$V5=="gneg","white",ifelse(cytoinfo$V5=="gpos25","gray25",ifelse(cytoinfo$V5=="gpos50","gray",
ifelse(cytoinfo$V5=="gpos75","darkgray",ifelse(cytoinfo$V5=="gpos100","black",ifelse(cytoinfo$V5=="acen","red","blue4")))
)))
cytoinfo$col = as.character(cytoinfo$col)
cyto_chr = cytoinfo[which(cytoinfo$V1==paste("chr",chr,sep="")),]
arm_chr = arminfo[which(arminfo$V1==paste("chr",chr,sep="")),]
len_chr1 = arminfo[which(arminfo$V1==paste("chr",1,sep="")),]$V2
chr_scale = arm_chr$V2/len_chr1

cna = cna %>% 
mutate(TYPE = if_else(COPY_CHANGE == "loss","Deletion",if_else(COPY_CHANGE == "neutral","CNN-LOH",if_else(COPY_CHANGE == "gain","Duplication","unknown")))) %>%
mutate(BEG_GRCh37 = as.numeric(start) * (10^6)) %>% mutate(END_GRCh37 = as.numeric(end) * (10^6)) %>%
mutate(CHROM = as.numeric(chr)) %>% mutate(ID = id)


cna$COLOR = ifelse(cna$TYPE=="Duplication","firebrick",ifelse(cna$TYPE=="Deletion","dodgerblue",ifelse(cna$TYPE=="CNN-LOH","chartreuse","gray")))
cna$DENSE = as.numeric(cna$CELL_FRAC)

cna_chr <- cna[which(cna$CHROM==chr),]
cna_chr_a <- cna_chr[which(cna_chr$TYPE=="Duplication"),]
cna_chr_d <- cna_chr[which(cna_chr$TYPE=="Deletion"),]
cna_chr_u <- cna_chr[which(cna_chr$TYPE=="CNN-LOH"),]
cna_chr_o = cna_chr[which(cna_chr$TYPE=="unknown"),]

sort_a = order(cna_chr_a$END_GRCh37,decreasing=T)
cna_chr_a = cna_chr_a[sort_a,]
sort_a = order(cna_chr_a$BEG_GRCh37)
cna_chr_a = cna_chr_a[sort_a,]
list_a = unique(cna_chr_a$ID)
list_cna_chr_a = rep(NA,length(list_a))
names(list_cna_chr_a) = list_a
for(i in list_a){
	cur_cna_chr_a = cna_chr_a[which(cna_chr_a$ID==i),]
	list_cna_a = rep(NA,nrow(cur_cna_chr_a))
 	for(k in 1:length(list_cna_a)){
		list_cna_a[k] = paste(cur_cna_chr_a[k,c("BEG_GRCh37","END_GRCh37","COLOR","DENSE")],collapse="-")
	}
	list_cna_chr_a[which(list_a==i)] = paste(c(i,paste(list_cna_a,collapse=",")),collapse=":")
}

sort_d = order(cna_chr_d$END_GRCh37,decreasing=T)
cna_chr_d = cna_chr_d[sort_d,]
sort_d = order(cna_chr_d$BEG_GRCh37)
cna_chr_d = cna_chr_d[sort_d,]
list_d = unique(cna_chr_d$ID)
list_cna_chr_d = rep(NA,length(list_d))
names(list_cna_chr_d) = list_d
for(i in list_d){
	cur_cna_chr_d = cna_chr_d[which(cna_chr_d$ID==i),]
	list_cna_d = rep(NA,nrow(cur_cna_chr_d))
	for(k in 1:length(list_cna_d)){
		list_cna_d[k] = paste(cur_cna_chr_d[k,c("BEG_GRCh37","END_GRCh37","COLOR","DENSE")],collapse="-")
 	}
	list_cna_chr_d[which(list_d==i)] = paste(c(i,paste(list_cna_d,collapse=",")),collapse=":")
}

sort_u = order(cna_chr_u$END_GRCh37,decreasing=T)
cna_chr_u = cna_chr_u[sort_u,]
sort_u = order(cna_chr_u$BEG_GRCh37)
cna_chr_u = cna_chr_u[sort_u,]
list_u = unique(cna_chr_u$ID)
list_cna_chr_u = rep(NA,length(list_u))
names(list_cna_chr_u) = list_u
for(i in list_u){
	cur_cna_chr_u = cna_chr_u[which(cna_chr_u$ID==i),]
	list_cna_u = rep(NA,nrow(cur_cna_chr_u))
	for(k in 1:length(list_cna_u)){
		list_cna_u[k] = paste(cur_cna_chr_u[k,c("BEG_GRCh37","END_GRCh37","COLOR","DENSE")],collapse="-")
	}
	list_cna_chr_u[which(list_u==i)] = paste(c(i,paste(list_cna_u,collapse=",")),collapse=":")
}

sort_o = order(cna_chr_o$END_GRCh37,decreasing=T)
cna_chr_o = cna_chr_o[sort_o,]
sort_o = order(cna_chr_o$BEG_GRCh37)
cna_chr_o = cna_chr_o[sort_o,]
list_o = unique(cna_chr_o$ID)
list_cna_chr_o = rep(NA,length(list_o))
names(list_cna_chr_o) = list_o
for(i in list_o){
	cur_cna_chr_o = cna_chr_o[which(cna_chr_o$ID==i),]
	list_cna_o = rep(NA,nrow(cur_cna_chr_o))
	for(k in 1:length(list_cna_o)){
		list_cna_o[k] = paste(cur_cna_chr_o[k,c("BEG_GRCh37","END_GRCh37","COLOR","DENSE")],collapse="-")
	}
	list_cna_chr_o[which(list_o==i)] = paste(c(i,paste(list_cna_o,collapse=",")),collapse=":")
}

#sep = "AAAA:1-1-white"
#list_cna_chr_adu = c(list_cna_chr_a,rep(sep,1),list_cna_chr_d,rep(sep,1),list_cna_chr_u,rep(sep,1),list_cna_chr_o)
list_cna_chr_adu = c(list_cna_chr_a,list_cna_chr_d,list_cna_chr_u,list_cna_chr_o)

list_cna_chr_adu_unq = paste0(unique(names(list_cna_chr_adu)),":")
names(list_cna_chr_adu_unq) = unique(names(list_cna_chr_adu))
for(i in 1:length(names(list_cna_chr_adu))){
	#print(unique(names(list_cna_chr_adu))[i])
	#print(list_cna_chr_adu[unique(names(list_cna_chr_adu))[i]])
	cur = gsub("^.*:","",list_cna_chr_adu[[i]])
	list_cna_chr_adu_unq[names(list_cna_chr_adu)[i]] = paste0(list_cna_chr_adu_unq[names(list_cna_chr_adu)[i]],cur,",")
}
list_cna_chr_adu = gsub("^,","",list_cna_chr_adu_unq)
print(list_cna_chr_adu_unq)


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

get_col = function(n){
	if(n==-1){ # 1: Duplication
		col = "firebrick3"
	}else if(n==-2){ # 2: Deletion
		col = "dodgerblue3"
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


## Draw ideogram ##
tag = paste0("CH_ideogram_chr",chr,".pdf")
pdf(file=tag,width=10,height=48)

plot(
NULL,
NULL,
xlim=c(-1.2,10.2),
ylim=c(-42,6),
axes=FALSE,
xlab="",
ylab=""
)

t = 1/6
for(i in 1:nrow(cyto_chr)){
	x_chr <- c(cyto_chr[i,"V2"]/arm_chr$V2*9,cyto_chr[i,"V2"]/arm_chr$V2*9,cyto_chr[i,"V3"]/arm_chr$V2*9,cyto_chr[i,"V3"]/arm_chr$V2*9)
	y_chr <- c(t-wid/2, t+wid/2, t+wid/2, t-wid/2)
	polygon(x_chr*chr_scale,y_chr,col=cyto_chr[i,"col"])
}

if(length(list_cna_chr_adu)!=0){
	for(i in 1:length(list_cna_chr_adu)){
		cur_id = unlist(strsplit(list_cna_chr_adu[i],":"))[1]
		cur_segs = unlist(strsplit(unlist(strsplit(list_cna_chr_adu[i],":"))[2],","))
		
		text(cur_id,,y_chr)
		
		for(s in cur_segs){
			stt = as.numeric(unlist(strsplit(s,"-"))[1])
			end = as.numeric(unlist(strsplit(s,"-"))[2])
			col = unlist(strsplit(s,"-"))[3]
			dense = as.numeric(unlist(strsplit(s,"-"))[4])
			if(is.na(dense)){
				dense = 0.1
			}
			print(dense)

			if(col == "firebrick"){
				col1 = "lightpink" #"lightsalmon"
				col2 = "firebrick4"
				cp_func = colorRampPalette(c("white",col1,col2))
				#cp_func = colorRampPalette(c("white",paste0(col,"1"),paste0(col,"2"),paste0(col,"3"),paste0(col,"4")))
				cp = cp_func(110)
				#col_2 = cp[min(round(log10(dense)*66+200)+10,110)]
				col_3 = cp[round(log10(dense)*33+100)+10]
			}else if(col == "dodgerblue"){
				col1 = "lightskyblue"
				col2 = "dodgerblue4"
				cp_func = colorRampPalette(c("white",col1,col2))
				#cp_func = colorRampPalette(c("white",paste0(col,"1"),paste0(col,"2"),paste0(col,"3"),paste0(col,"4")))
				cp = cp_func(110)
				#col_2 = cp[min(round(log10(dense)*66+200)+10,110)]
				col_3 = cp[round(log10(dense)*33+100)+10]
			}else if(col == "chartreuse"){
				col1 = "darkseagreen1" #"chartreuse"
				col2 = "chartreuse4"
				cp_func = colorRampPalette(c("white",col1,col2))
				#cp_func = colorRampPalette(c("white",paste0(col,"1"),paste0(col,"2"),paste0(col,"3"),paste0(col,"4")))
				cp = cp_func(110)
				#col_2 = cp[min(round(log10(dense)*66+200)+10,110)]
				col_3 = cp[round(log10(dense)*33+100)+10]
			}else{
				col_3 = col
			}
			l_width = 0.2*6 + 2 # (dense)*6 + 2
			segments(stt*chr_scale/arm_chr$V2*9,-(i+3)*wid/2,end*chr_scale/arm_chr$V2*9,-(i+3)*wid/2,col=col_3,lwd=l_width)
			
			x_range = 2*(10^7)			
			cur_mut = mut %>% filter(id == cur_id & as.numeric(Chr) == chr & as.numeric(Start) >= stt - x_range & as.numeric(End) <= end + x_range)
			
			if(nrow(cur_mut) >= 2){
				in_col = get_col(2)
			}else if(nrow(cur_mut) == 1){
				cur_type = cur_mut$Merge_Func
				if(cur_type == "nonsynonymous SNV" | cur_type == "stoploss"){
					in_col = get_col(7)
				}else if(cur_type == "splicing"){
					in_col = get_col(5)
				}else if(cur_type == "frameshift deletion" | cur_type == "frameshift insertion"){
					in_col = get_col(4)
				}else if(cur_type == "stopgain"){
					in_col = get_col(3)
				}else if(cur_type == "nonframeshift deletion" | cur_type == "nonframeshift insertion"){
					in_col = get_col(6)
				}else{
					in_col = "white"
				}
			}

			print(cur_id)
			print(s)

			if(nrow(cur_mut) >= 1){
				print("mut+")
				print(cur_mut)
				for (j in 1:nrow(cur_mut)){
					x_c = as.numeric(cur_mut$Start[j])*chr_scale/arm_chr$V2*9
					y_c = -(i+3)*wid/2
					#theta <- seq(-pi, pi, length=100)
					#polygon(x_c+0.06*cos(theta),y_c+0.06*sin(theta),col="orange",border="transparent")
					#polygon(x_c+0.045*cos(theta),y_c+0.045*sin(theta),col=in_col,border="transparent")
					mut_col = "red"
					segments(x_c,y_c-0.075,x_c,y_c+0.075,col=mut_col,lwd=l_width/4)
					segments(x_c-0.05,y_c-0.05,x_c+0.05,y_c+0.05,col=mut_col,lwd=l_width/4)
					segments(x_c+0.05,y_c-0.05,x_c-0.05,y_c+0.05,col=mut_col,lwd=l_width/4)
					text(-0.2,y_c,cur_id,cex=0.5,adj=1)
				}
			}else{
				print("mut-")
				y_c = -(i+3)*wid/2
				text(-0.2,y_c,cur_id,cex=0.5,adj=1)
			}
		}
	}
}

dev.off()



## draw legend ##
cp_func_1 = colorRampPalette(c("white", "lightpink","firebrick4"))
cp_func_2 = colorRampPalette(c("white", "lightskyblue","dodgerblue4"))
cp_func_3 = colorRampPalette(c("white", "darkseagreen1","chartreuse4"))
cp_1 = cp_func_1(110)
cp_2 = cp_func_2(110)
cp_3 = cp_func_3(110)
log_scale = seq(-3,0,length=100)
cs_1 = cp_1[round(log_scale*33+100)+10]
cs_2 = cp_2[round(log_scale*33+100)+10]
cs_3 = cp_3[round(log_scale*33+100)+10]
ws = ((1:100)/100)*6 + 2

pdf("ideogram_legend.pdf",width=10,height=20)
plot(NULL,NULL,
xlim=c(-1.2,10.2),
ylim=c(-15,5),
axes=FALSE,
xlab="",ylab="") 

ws = 0.2*6 + 2 # (dense)*6 + 2                                                                                                                           
for(i in 2:100){
	ls1 = log_scale[i-1]
	ls2 = log_scale[i]
	segments(ls1+3,-10,ls2+3,-10,lwd=ws,col=cs_1[i],border=NULL)
	segments(ls1+3,-10-0.3,ls2+3,-10-0.3,lwd=ws,col=cs_2[i],border=NULL)
	segments(ls1+3,-10-0.6,ls1+3,-10-0.6,lwd=ws,col=cs_3[i],border=NULL)
}

text(0+3.2,-10,"dupulication",adj=0,cex=1.2,font=2)
text(0+3.2,-10-0.3,"deletion",adj=0,cex=1.2,font=2)
text(0+3.2,-10-0.6,"UPD",adj=0,cex=1.2,font=2)
 
text(-1.5+3,-10+0.6,"Cell fraction (%)",cex=1.2,font=2)
text(-3+3,-10+0.3,"0.1",adj=0.4,cex=1,font=2)
text(-2+3,-10+0.3,"1",adj=0.4,cex=1,font=2)
text(-1+3,-10+0.3,"10",adj=0.4,cex=1,font=2)
text(0+3,-10+0.3,"100",adj=0.4,cex=1,font=2)
  

text(0,-11.5+0.4,"Type of mutations",adj=0,cex=1.2,font=2)
theta <- seq(-pi, pi, length=100)
#missense
x_c = 0
y_c = -11.5
polygon(x_c+0.06*cos(theta),y_c+0.06*sin(theta),col="orange",border="transparent")
polygon(x_c+0.045*cos(theta),y_c+0.045*sin(theta),col="deepskyblue1",border="transparent")
text(x_c+0.2,y_c,"Missense",adj=0,cex=1.2,font=2)
#missense
x_c = 0
y_c = -11.5-0.3*1
polygon(x_c+0.06*cos(theta),y_c+0.06*sin(theta),col="orange",border="transparent")
polygon(x_c+0.045*cos(theta),y_c+0.045*sin(theta),col="yellow",border="transparent")
text(x_c+0.2,y_c,"Inframe indel",adj=0,cex=1.2,font=2)
#missense
x_c = 0
y_c = -11.5-0.3*2
polygon(x_c+0.06*cos(theta),y_c+0.06*sin(theta),col="orange",border="transparent")
polygon(x_c+0.045*cos(theta),y_c+0.045*sin(theta),col="purple",border="transparent")
text(x_c+0.2,y_c,"Splice-site",adj=0,cex=1.2,font=2)
#missense
x_c = 0
y_c = -11.5-0.3*3
polygon(x_c+0.06*cos(theta),y_c+0.06*sin(theta),col="orange",border="transparent")
polygon(x_c+0.045*cos(theta),y_c+0.045*sin(theta),col="salmon",border="transparent")
text(x_c+0.2,y_c,"Frameshift indel",adj=0,cex=1.2,font=2)
#missense
x_c = 0
y_c = -11.5-0.3*4
polygon(x_c+0.06*cos(theta),y_c+0.06*sin(theta),col="orange",border="transparent")
polygon(x_c+0.045*cos(theta),y_c+0.045*sin(theta),col="red",border="transparent")
text(x_c+0.2,y_c,"Stop-gain",adj=0,cex=1.2,font=2)
#missense
x_c = 0
y_c = -11.5-0.3*5
polygon(x_c+0.06*cos(theta),y_c+0.06*sin(theta),col="orange",border="transparent")
polygon(x_c+0.045*cos(theta),y_c+0.045*sin(theta),col="brown",border="transparent")
text(x_c+0.2,y_c,"Multiple",adj=0,cex=1.2,font=2)

dev.off()



