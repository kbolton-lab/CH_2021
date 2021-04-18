
## import data
mut_path = "MatMed_2021_mut_table.txt"
cna_path = "MatMed_2021_cna_table.csv"  
id = as.character(1:11234)


# get info of chr length
arminfo = read.table(file="cytoBand/arm_chr_length.txt")
cent = arminfo[1:22,4]

mut = read.table(mut_path,header=T,stringsAsFactor=F,sep="\t",quote="",colClass="character")

age = rep(NA,length(id))
for(i in 1:length(id)){
age[i] = mut$age[which(mut$id==id[i])[1]]
}
age = as.numeric(age); names(age) = id
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
label_cna = paste0(sort(rep(1:22,4*2)),":",sort(rep(c("CNN-LOH","Duplication","Deletion","Unknown"),2)),":",c("p","q"))
mat_cna =  matrix(0,ncol=length(label_cna),nrow=length(id))
colnames(mat_cna) = label_cna; rownames(mat_cna) = id
n_cna = rep(0,length(id)); names(n_cna) = rownames(mat_cna)

for(i in 1:nrow(cna)){
	cur_id = cna$id[i]
	
	cur_type = "Unknown"
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
	if(cur_type == "Unknown") mat_cna[cur_id,cur_cna] = 0.000001
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



n_obs = sum(all_n_mut + all_n_cna >= 2)

bin = as.numeric(names(table(age)))
N = 50000
n_sim_bin = rep(0,N)
for(i in 1:length(bin)){
bin_n = sum(age==bin[i])
alt_n = sum(all_n_mut[which(age==bin[i])] + all_n_cna[which(age==bin[i])])
for(j in 1:N){
n_sim_bin[j] = n_sim_bin[j] + sum(table(sample(bin_n,alt_n,replace=T))>=2)
}
print(bin[i])
}


## P value for cooccurrences of multiple alterations : 0.01 ~ 0.005
sum(n_sim_bin>=n_obs)/length(n_sim_bin)



mat_mut_2 = apply(mat_mut,c(1,2),function(x){as.numeric(x!=0)})
mat_cna_2 = apply(mat_cna,c(1,2),function(x){as.numeric(x!=0)})


A = matrix(1,nrow=ncol(mat_mut_2),ncol=ncol(mat_cna_2))
rownames(A) = colnames(mat_mut_2); colnames(A) = colnames(mat_cna_2)
B = A

A["TET2","4:CNN-LOH:q"] = 0; A["TET2","4:Deletion:q"] = 0; A["TET2","4:Unknown:q"] = 0
A["DNMT3A","2:CNN-LOH:p"] = 0; A["DNMT3A","2:Deletion:p"] = 0; A["DNMT3A","2:Unknown:p"] = 0
A["TP53","17:CNN-LOH:p"] = 0; A["TP53","17:Deletion:p"] = 0; A["TP53","17:Unknown:p"] = 0
A["JAK2","9:CNN-LOH:p"] = 0; A["JAK2","9:Duplication:p"] = 0; A["JAK2","9:Unknown:p"] = 0

mat_A = mat_mut_2 %*% A %*% t(mat_cna_2)
res_A = sum(diag(mat_A)>=1)

bin = as.numeric(names(table(age)))
N = 1000
n_sim_bin = rep(0,N)
for(i in 1:length(bin)){
	bin_n = sum(age==bin[i])
	cur_mat_mut = matrix(mat_mut_2[which(age==bin[i]),],ncol=ncol(mat_mut_2))
	cur_mat_cna = matrix(mat_cna_2[which(age==bin[i]),],ncol=ncol(mat_cna_2))
	
	for(j in 1:N){
		cur_mat_mut = matrix(cur_mat_mut[sample(1:nrow(cur_mat_mut),nrow(cur_mat_mut),replace=F),],ncol=ncol(cur_mat_mut))
		cur_mat_cna = matrix(cur_mat_cna[sample(1:nrow(cur_mat_cna),nrow(cur_mat_cna),replace=F),],ncol=ncol(cur_mat_cna))
		n_sim_bin[j] = n_sim_bin[j] + sum(diag(cur_mat_mut %*% A %*% t(cur_mat_cna))>=1)
	}

print(bin[i])
}

## P value for cooccurrences of SNV and CNA (excluding significant combinations) : ~ 0.006
sum(res_A<=n_sim_bin)/length(n_sim_bin)



mat_B = mat_mut_2 %*% B %*% t(mat_cna_2)
res_B = sum(diag(mat_B)>=1)

bin = as.numeric(names(table(age)))
N = 1000
n_sim_bin = rep(0,N)
for(i in 1:length(bin)){
	bin_n = sum(age==bin[i])
	cur_mat_mut = matrix(mat_mut_2[which(age==bin[i]),],ncol=ncol(mat_mut_2))
	cur_mat_cna = matrix(mat_cna_2[which(age==bin[i]),],ncol=ncol(mat_cna_2))
	
	for(j in 1:N){
		cur_mat_mut = matrix(cur_mat_mut[sample(1:nrow(cur_mat_mut),nrow(cur_mat_mut),replace=F),],ncol=ncol(cur_mat_mut))
		cur_mat_cna = matrix(cur_mat_cna[sample(1:nrow(cur_mat_cna),nrow(cur_mat_cna),replace=F),],ncol=ncol(cur_mat_cna))
		n_sim_bin[j] = n_sim_bin[j] + sum(diag(cur_mat_mut %*% B %*% t(cur_mat_cna))>=1)
	}

print(bin[i])
}

## P value for cooccurrences of multiple alterations (including significant combinations) : <0.001
sum(res_B<=n_sim_bin)/length(n_sim_bin)



