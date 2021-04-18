
## import data ##
cna_path = "Loh_2018_CNA.txt" ## from Loh et al. (2020) nature.
cna = read.table(cna_path,header=T,sep="\t",stringsAsFactor=F,quote="",colClass="character")
cna = cna[which(is.element(cna$AGE,c("60-65","65-70","70-75"))),]
id = unique(cna$ID)


## import marker ##
mrk_path = "array_marker.txt"
mrk = read.table(mrk_path,header=T,sep="\t",stringsAsFactor=F,quote="")


lrr_tmp = unlist(strsplit(cna$LRR..SE.," "))  
lrr_tmp = lrr_tmp[(1:(length(lrr_tmp)/2))*2-1]
cna$lrr = as.numeric(lrr_tmp)

baf_tmp = unlist(strsplit(cna$BAF..SE.," "))  
baf_tmp = baf_tmp[(1:(length(baf_tmp)/2))*2-1]
cna$baf = as.numeric(baf_tmp)

cna = cna[which(cna$CHR != "X"),]
cna = cna[order(as.numeric(cna$End..GRCh37.),decreasing=F),]
cna = cna[order(as.numeric(cna$Start..GRCh37.),decreasing=F),]
cna = cna[order(as.numeric(cna$CHR),decreasing=F),]

gain_loss = cna[(cna$COPY_CHANGE == "gain" | cna$COPY_CHANGE == "loss"),]
upd = cna[(cna$COPY_CHANGE == "neutral"),]


chr_def = matrix(nrow=0,ncol=5)
for(i in 1:22){
	chr_def = rbind(chr_def,
	c(
	i,
	min(as.numeric(mrk$Position[which(mrk$Chromosome==i)])),
	max(as.numeric(mrk$Position[which(mrk$Chromosome==i)])),
	sum(mrk$Chromosome==i),
	0
	)
	)
}

gain_loss_res = matrix(nrow=0,ncol=6)

for(s in 1:length(id)){
	print(paste0(s,"/",length(id)))

	c_id =id[s]
	c_gain_loss = gain_loss[which(gain_loss$ID == c_id),]
	c_res = matrix(nrow=0,ncol=5)

	for(i in 1:22){
		chr_c_gain_loss = c_gain_loss[which(c_gain_loss$CHR == i),]

		if(nrow(chr_c_gain_loss) == 0){

			c_res = rbind(c_res,chr_def[i,])

		}else{

			c_stt = as.numeric(chr_c_gain_loss$Start..GRCh37.)
			c_end = as.numeric(chr_c_gain_loss$End..GRCh37.)+1
			c_bks = sort(c(0,c_stt,c_end,chr_def[i,3]+1),decreasing=F)

			bks_mat = cbind(
			rep(i,length(c_bks)-1),
			c_bks[-length(c_bks)],
			(c_bks-1)[-1],
			rep(NA,length(c_bks)-1),
			rep(NA,length(c_bks)-1)
			)

			for(j in 1:nrow(bks_mat)){
				bks_mat[j,4] = sum(mrk$Chromosome ==i & mrk$Position >= bks_mat[j,2] & mrk$Position <= bks_mat[j,3])

				c_lrr = chr_c_gain_loss$lrr[which(chr_c_gain_loss$Start..GRCh37. == bks_mat[j,2])]
				if(length(c_lrr)==0){c_lrr = 0}
				if(c_lrr > 0){
					c_lrr = 1
				}else if(c_lrr < 0){
					c_lrr = -1
				}
				bks_mat[j,5] = c_lrr
				
				
				if(bks_mat[j,4] != 0){
					bks_mat[j,2] = min(mrk$Position[which(mrk$Chromosome == i & mrk$Position >= bks_mat[j,2])])
					bks_mat[j,3] = max(mrk$Position[which(mrk$Chromosome == i & mrk$Position <= bks_mat[j,3])])
				}
			}

			c_res = rbind(c_res,bks_mat)

		}

	}

	gain_loss_res = rbind(gain_loss_res,
	cbind(rep(c_id,nrow(c_res)),c_res)
	)

}

gain_loss_res = gain_loss_res[which(gain_loss_res[,5] != "0"),]
gain_res = gain_loss_res[which(gain_loss_res[,6] == 1),]
loss_res = gain_loss_res[which(gain_loss_res[,6] == -1),]
loss_res[,6] = -1 * as.numeric(loss_res[,6])
write.table(gain_res,file="gain_loh_for_part.txt",col.name=F,row.name=F,sep="\t",quote=F)
write.table(loss_res,file="loss_loh_for_part.txt",col.name=F,row.name=F,sep="\t",quote=F)



upd_res = matrix(nrow=0,ncol=6)

for(s in 1:length(id)){
	print(paste0(s,"/",length(id)))

	c_id =id[s]
	c_upd = upd[which(upd$ID == c_id),]
	c_res = matrix(nrow=0,ncol=5)

	for(i in 1:22){
		chr_c_upd = c_upd[which(c_upd$CHR == i),]

		if(nrow(chr_c_upd) == 0){

			c_res = rbind(c_res,chr_def[i,])

		}else{

			c_stt = as.numeric(chr_c_upd$Start..GRCh37.)
			c_end = as.numeric(chr_c_upd$End..GRCh37.)+1
			c_bks = sort(c(0,c_stt,c_end,chr_def[i,3]+1),decreasing=F)

			bks_mat = cbind(
			rep(i,length(c_bks)-1),
			c_bks[-length(c_bks)],
			(c_bks-1)[-1],
			rep(NA,length(c_bks)-1),
			rep(NA,length(c_bks)-1)
			)

			for(j in 1:nrow(bks_mat)){
				bks_mat[j,4] = sum(mrk$Chromosome ==i & mrk$Position >= bks_mat[j,2] & mrk$Position <= bks_mat[j,3])

				c_lrr = chr_c_upd$baf[which(chr_c_upd$Start..GRCh37. == bks_mat[j,2])]
				if(length(c_lrr)==0){c_lrr = 0}
				if(c_lrr != 0){
					c_lrr = 1
				}

				bks_mat[j,5] = c_lrr

				
				if(bks_mat[j,4] != 0){
					bks_mat[j,2] = min(mrk$Position[which(mrk$Chromosome == i & mrk$Position >= bks_mat[j,2])])
					bks_mat[j,3] = max(mrk$Position[which(mrk$Chromosome == i & mrk$Position <= bks_mat[j,3])])
				}
			}

			c_res = rbind(c_res,bks_mat)

		}

	}

	upd_res = rbind(upd_res,
	cbind(rep(c_id,nrow(c_res)),c_res)
	)

}

upd_res = upd_res[which(upd_res[,5] != "0"),]
upd_res = upd_res[which(upd_res[,6] == 1),]
write.table(upd_res,file="upd_loh_for_part.txt",col.name=F,row.name=F,sep="\t",quote=F)


