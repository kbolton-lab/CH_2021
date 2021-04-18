
## import CH data ##
mut_path = ~~ ## please import mutation data
cna_path = ~~ ## please import cna data
mut = read.table(mut_path,header=T,stringsAsFactor=F,sep="\t",quote="",colClass="character")
cna = read.table(cna_path,header=T,stringsAsFactor=F,sep="\t",quote="",colClass="character")

id = as.character(1:11234)

## get centromere position ##
arminfo = read.table(file=paste0(path,"/bbj/arm_chr_length.txt"))
cent = arminfo[1:22,4]


## make mutation table ##
mut = read.table(mut_path,header=T,stringsAsFactor=F,quote="",sep="\t",comment.char="")
mat_mut = matrix(0,ncol=length(unique(mut[,8])),nrow=length(id))
colnames(mat_mut) = unique(mut[,8]); rownames(mat_mut) = id
n_mut = rep(0,length(id)); names(n_mut) = id

for(i in 1:nrow(mut)){
	cur_id = mut[i,1]
	cur_gene = mut[i,8]
	cur_vaf = mut[i,57]
	if(! is.element(cur_id,id)){
		next
	}
	if(cur_vaf > mat_mut[cur_id,cur_gene]){
		mat_mut[cur_id,cur_gene] = cur_vaf
	}
	n_mut[cur_id] = n_mut[cur_id] + 1
}


## make cna table ##
cna = read.table(cna_path,header=T,stringsAsFactor=F,quote="",sep="\t",comment.char="")
label_cna = paste0(sort(rep(1:22,4*2)),"_",sort(rep(c("CNN-LOH","Duplication","Deletion","unknown"),2)),"_",c("p","q"))
mat_cna =  matrix(0,ncol=length(label_cna),nrow=length(id))
colnames(mat_cna) = label_cna; rownames(mat_cna) = id
n_cna = rep(0,length(id)); names(n_cna) = id

for(i in 1:nrow(cna)){
	cur_id = cna$IDY[i]
	if(! is.element(cur_id,id)) next
	
	## get type of cna ##
	cur_type = "unknown"
	if(cna$COPY_CHANGE[i] == "loss"){
		cur_type = "Deletion"
	}else if(cna$COPY_CHANGE[i] == "neutral"){
		cur_type = "CNN-LOH"
	}else if(cna$COPY_CHANGE[i] == "gain"){
		cur_type = "Duplication"
	}
	
	## ignore unclassifiable cna
	if(cur_type == "unknown") next
	
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

	cur_cf = cna$CELL_FRAC[i]
	if(mat_cna[cur_id,cur_cna] < cur_cf & !is.na(cur_cf)){
		mat_cna[cur_id,cur_cna] = cur_cf
	}
	n_cna[cur_id] = n_cna[cur_id] + 1
}


## who has mutation
## who has cna
## who has both
mut_mat = (mat_mut != 0)
cna_mat = (mat_cna != 0)
is_mut = as.numeric(apply(mat_mut!=0,1,sum)!=0)
is_cna = as.numeric(apply(mat_cna!=0,1,sum)!=0)
is_both = is_mut*is_cna


## get largerst clone size
cs_mut = as.numeric(apply(mat_mut,1,max))
cs_cna = as.numeric(apply(mat_cna,1,max))
cs_size = cs_mut*2*(cs_mut*2>cs_cna) + cs_cna*(cs_mut*2<=cs_cna)
cs_size[which(cs_size>1)] = 1


## get total number of alterations
n_alt = n_mut + n_cna


## draw Fig.1e ##
pdf("mnax_clone_size_by n_alteration.pdf",width=10,height=5)
plot(NULL,NULL,xlim=c(-0.5,10.5),ylim=c(-3,0),axes=F,xlab="",ylab="")

cs1_list = list()
cs2_list = list()
for(i in 1:4){
cs1 = cs_size[which(n_mut + n_cna==i & ((n_mut>=1 & n_cna==0) | (n_mut==0 & n_cna>=1)))]
cs2 = cs_size[which(n_mut + n_cna==i & n_mut>=1 & n_cna>=1)]
cs1[which(cs1>1)] = cs2[which(cs2>1)] = 1

cs1_list[[i]] = cs1
cs2_list[[i]] = cs2

cs1_qt = quantile(cs1)
cs1_lw = cs1_qt[2] - 1.5*(cs1_qt[4]-cs1_qt[2])
cs1_lw = max(min(cs1[which(cs1>=cs1_lw)]),cs1_lw)
cs1_up = cs1_qt[4] + 1.5*(cs1_qt[4]-cs1_qt[2])
cs1_up = min(max(cs1[which(cs1<=cs1_up)]),cs1_up)
cs1_x = i + runif(length(cs1),-.15,.15)-0.2

cs2_qt = quantile(cs2)
cs2_lw = cs2_qt[2] - 1.5*(cs2_qt[4]-cs2_qt[2])
cs2_lw = max(min(cs2[which(cs2>=cs1_lw)]),cs2_lw)
cs2_up = cs2_qt[4] + 1.5*(cs2_qt[4]-cs2_qt[2])
cs2_up = min(max(cs2[which(cs2<=cs2_up)]),cs2_up)
cs2_x = i + runif(length(cs2),-.15,.15)+0.2

i1 = i-0.2
points(cs1_x,log10(cs1),pch=20,cex=0.5,col=rgb(0,0,1,alpha=0.15))
rect(i1-0.1,log10(cs1_qt[2]),i1+0.1,log10(cs1_qt[4]),lwd=2,border="black")
segments(i1-0.1,log10(cs1_qt[3]),i1+0.1,log10(cs1_qt[3]),lwd=2,col="black")
segments(i1,log10(cs1_up),i1,log10(cs1_qt[4]),lwd=2,col="black")
segments(i1,log10(cs1_lw),i1,log10(cs1_qt[2]),lwd=2,col="black")
segments(i1-0.05,log10(cs1_up),i1+0.05,log10(cs1_up),lwd=2,col="black")
segments(i1-0.05,log10(cs1_lw),i1+0.05,log10(cs1_lw),lwd=2,col="black")

i2 = i+0.2
points(cs2_x,log10(cs2),pch=20,cex=0.5,col=rgb(1,0,1,alpha=0.3))
rect(i2-0.1,log10(cs2_qt[2]),i2+0.1,log10(cs2_qt[4]),lwd=2,border="black")
segments(i2-0.1,log10(cs2_qt[3]),i2+0.1,log10(cs2_qt[3]),lwd=2,col="black")
segments(i2,log10(cs2_up),i2,log10(cs2_qt[4]),lwd=2,col="black")
segments(i2,log10(cs2_lw),i2,log10(cs2_qt[2]),lwd=2,col="black")
segments(i2-0.05,log10(cs2_up),i2+0.05,log10(cs2_up),lwd=2,col="black")
segments(i2-0.05,log10(cs2_lw),i2+0.05,log10(cs2_lw),lwd=2,col="black")

}

lat1 = c(0.001,0.01,0.1,1)
lat2 = c(0.001*(2:9),0.01*(2:9),0.1*(2:9),1)
segments(-0.05,log10(lat1),0,log10(lat1))
segments(-0.025,log10(lat2),0,log10(lat2))
segments(5.05,log10(lat1),5,log10(lat1))
segments(5.025,log10(lat2),5,log10(lat2))
text(-0.1,log10(lat1),lat1,cex=0.8,adj=1)
text(5.1,log10(lat1),lat1,cex=0.8,adj=0)
segments(0,log10(0.001),0,log10(1))
segments(5,log10(0.001),5,log10(1))
segments(0,log10(0.001),5,log10(0.001))

dev.off()


## comparison: different number of alterations ##
wilcox.test(c(cs1_list[[1]],cs2_list[[1]]),c(cs1_list[[2]],cs2_list[[2]]))
wilcox.test(c(cs1_list[[2]],cs2_list[[2]]),c(cs1_list[[3]],cs2_list[[3]]))
wilcox.test(c(cs1_list[[3]],cs2_list[[3]]),c(cs1_list[[4]],cs2_list[[4]]))


## comparison: either vs both ##
wilcox.test(cs1_list[[1]],cs2_list[[1]])
wilcox.test(cs1_list[[2]],cs2_list[[2]])
wilcox.test(cs1_list[[3]],cs2_list[[3]])
wilcox.test(cs1_list[[4]],cs2_list[[4]])


## End of analysis ##
q()



