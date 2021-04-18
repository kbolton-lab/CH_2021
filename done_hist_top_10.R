
## Draw ED Fig.2a,b ##

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
label_cna = paste0(sort(rep(1:22,4*2)),":",sort(rep(c("CNN-LOH","Duplication","Deletion","Unknown"),2)),":",c("p","q"))
mat_cna =  matrix(0,ncol=length(label_cna),nrow=length(id))
colnames(mat_cna) = label_cna; rownames(mat_cna) = id

for(i in 1:nrow(cna)){
	cur_id = cna$id[i]
	if(! is.element(cur_id,id)){
		next
	}

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
cna_count[,1] = apply(mat_cna[,unk],1,sum) # Unknown

mat_cna = mat_cna[,which(apply(mat_cna,2,sum)!=0)]
mat_cna = mat_cna[,-grep("Unknown",colnames(mat_cna))]

mut_10 = sort(apply(mat_mut!=0,2,sum),decreasing=T)[1:20]
cna_10 = sort(apply(mat_cna!=0,2,sum),decreasing=T)[1:20]

getCnaLabel = function(s){
	splited = unlist(strsplit(s,":"))
	type=NA
	prf = ifelse(splited[2]=="CNN-LOH","",ifelse(splited[2]=="Duplication","+",ifelse(splited[2]=="Deletion","del(","Unknown")))
	suf = ifelse(splited[2]=="CNN-LOH","UPD",ifelse(splited[2]=="Duplication","",ifelse(splited[2]=="Deletion",")","Unknown")))
	return(paste0(prf,splited[1],splited[3],suf))
}
names(cna_10) = sapply(names(cna_10),getCnaLabel)
names(mut_10) = gsub("U2AF1;U2AF1L5","U2AF1",names(mut_10))
#mut_10*100/length(id)
#cna_10*100/length(id)

pdf("hist_top10_mutation_cna.pdf",width=10,height=12)
par(mfrow=c(2,1))

plot(NULL,NULL,xlim=c(-10,34),ylim=c(-15,33),axes=FALSE,xlab="",ylab="")
lat1 = 500*(0:3)/50
lat2 = setdiff(100*(0:16)/50,lat1)
for(i in 1:length(lat1)){
	segments(0,lat1,-0.3,lat1)
	text(-0.5,lat1,lat1*50,adj=1,cex=0.5)
}
for(i in 1:length(lat2)){
	segments(0,lat2,-0.15,lat2)
}
segments(0,0,0,32)
text(-2.5,15,"Number of subjects",srt=90,cex=0.7)

for(i in 1:length(mut_10)){
	x= i*3/2
	text(x-0.5,-1.5,names(mut_10)[i],srt=45,adj=1,cex=0.5,font=4)
	rect(x-1,0,x,mut_10[i]/50,col="gray",border="transparent")
}

plot(NULL,NULL,xlim=c(-10,34),ylim=c(-15,33),axes=FALSE,xlab="",ylab="")
lat1 = 500*(0:3)/50
lat2 = setdiff(100*(0:16)/50,lat1)
for(i in 1:length(lat1)){
	segments(0,lat1,-0.3,lat1)
	text(-0.5,lat1,lat1*50,adj=1,cex=0.5)
}
for(i in 1:length(lat2)){
	segments(0,lat2,-0.15,lat2)
	}
segments(0,0,0,32)
text(-2.5,15,"Number of subjects",srt=90,cex=0.7)

for(i in 1:length(cna_10)){
	x= i*3/2
	text(x-0.5,-1.5,names(cna_10)[i],srt=45,adj=1,cex=0.5,font=2)
	rect(x-1,0,x,cna_10[i]/50,col="gray",border="transparent")
}

dev.off()

## End of analysis
q()



