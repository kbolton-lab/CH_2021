
## Draw ED Fig.2e-i ##

mut_file = ## import mutation data on JGA

id = as.character(1:11234)

dat = read.table(mut_file,header=T,stringsAsFactor=F,sep="\t",quote="")
dat = dat[which(is.element(dat$id,id)),]

dat_tet2 = dat[which(dat$Gene.refGene=="TET2"),]
dat_dnmt3a = dat[which(dat$Gene.refGene=="DNMT3A"),]
dat_asxl1 = dat[which(dat$Gene.refGene=="ASXL1"),]
dat_srsf2 = dat[which(dat$Gene.refGene=="SRSF2"),]
dat_cbl = dat[which(dat$Gene.refGene=="CBL"),]

id_tet2_asxl1 = unique(intersect(dat_tet2$id,dat_asxl1$id))
id_tet2_dnmt3a = unique(intersect(dat_tet2$id,dat_dnmt3a$id))
id_tet2_srsf2 = unique(intersect(dat_tet2$id,dat_srsf2$id))
id_asxl1_srsf2 = unique(intersect(dat_asxl1$id,dat_srsf2$id))
id_asxl1_cbl = unique(intersect(dat_asxl1$id,dat_cbl$id))


vaf_tet2_dnmt3a = matrix(NA,nrow=length(id_tet2_dnmt3a),ncol=2)
for(i in 1:length(id_tet2_dnmt3a)){
	cur_id = id_tet2_dnmt3a[i]
	vafs_dnmt3a = as.numeric(dat_dnmt3a[which(dat_dnmt3a$id == cur_id),"misRate"])
	vafs_tet2 = as.numeric(dat_tet2[which(dat_tet2$id == cur_id),"misRate"])
	if(length(vafs_dnmt3a != 0)){
		vaf_dnmt3a = max(vafs_dnmt3a)
	}else{
		vaf_dnmt3a = 0
	}
	if(length(vafs_tet2 != 0)){
		vaf_tet2 = max(vafs_tet2)
	}else{
		vaf_tet2 = 0
	}
	vaf_tet2_dnmt3a[i,] = c(vaf_tet2,vaf_dnmt3a)
}

vaf_tet2_asxl1 = matrix(NA,nrow=length(id_tet2_asxl1),ncol=2)
for(i in 1:length(id_tet2_asxl1)){
	cur_id = id_tet2_asxl1[i]
	vafs_asxl1 = as.numeric(dat_asxl1[which(dat_asxl1$id == cur_id),"misRate"])
	vafs_tet2 = as.numeric(dat_tet2[which(dat_tet2$id == cur_id),"misRate"])
	if(length(vafs_asxl1 != 0)){
		vaf_asxl1 = max(vafs_asxl1)
	}else{
		vaf_asxl1 = 0
	}
	if(length(vafs_tet2 != 0)){
		vaf_tet2 = max(vafs_tet2)
	}else{
		vaf_tet2 = 0
	}
	vaf_tet2_asxl1[i,] = c(vaf_tet2,vaf_asxl1)
}

vaf_tet2_srsf2 = matrix(NA,nrow=length(id_tet2_srsf2),ncol=2)
for(i in 1:length(id_tet2_srsf2)){
	cur_id = id_tet2_srsf2[i]
	vafs_srsf2 = as.numeric(dat_srsf2[which(dat_srsf2$id == cur_id),"misRate"])
	vafs_tet2 = as.numeric(dat_tet2[which(dat_tet2$id == cur_id),"misRate"])
	if(length(vafs_srsf2 != 0)){
		vaf_srsf2 = max(vafs_srsf2)
	}else{
		vaf_srsf2 = 0
	}
	if(length(vafs_tet2 != 0)){
		vaf_tet2 = max(vafs_tet2)
	}else{
		vaf_tet2 = 0
	}
	vaf_tet2_srsf2[i,] = c(vaf_tet2,vaf_srsf2)
}


vaf_asxl1_srsf2 = matrix(NA,nrow=length(id_asxl1_srsf2),ncol=2)
for(i in 1:length(id_asxl1_srsf2)){
	cur_id = id_asxl1_srsf2[i]
	vafs_asxl1 = as.numeric(dat_asxl1[which(dat_asxl1$id == cur_id),"misRate"])
	vafs_srsf2 = as.numeric(dat_srsf2[which(dat_srsf2$id == cur_id),"misRate"])
	if(length(vafs_asxl1 != 0)){
		vaf_asxl1 = max(vafs_asxl1)
	}else{
		vaf_asxl1 = 0
	}
	if(length(vafs_srsf2 != 0)){
		vaf_srsf2 = max(vafs_srsf2)
	}else{
		vaf_srsf2 = 0
	}
	vaf_asxl1_srsf2[i,] = c(vaf_asxl1,vaf_srsf2)
}


vaf_asxl1_cbl = matrix(NA,nrow=length(id_asxl1_cbl),ncol=2)
for(i in 1:length(id_asxl1_cbl)){
	cur_id = id_asxl1_cbl[i]
	vafs_asxl1 = as.numeric(dat_asxl1[which(dat_asxl1$id == cur_id),"misRate"])
	vafs_cbl = as.numeric(dat_cbl[which(dat_cbl$id == cur_id),"misRate"])
	if(length(vafs_asxl1 != 0)){
		vaf_asxl1 = max(vafs_asxl1)
	}else{
		vaf_asxl1 = 0
	}
	if(length(vafs_cbl != 0)){
		vaf_cbl = max(vafs_cbl)
	}else{
		vaf_cbl = 0
	}
	vaf_asxl1_cbl[i,] = c(vaf_asxl1,vaf_cbl)
}


pdf("tet2_dnmt3a.pdf")
m = 0.6
plot(NULL,NULL,xlim=c(-0.15,m+0.15),ylim=c(-0.15,m+0.15),axes=FALSE,xlab="",ylab="")
points(vaf_tet2_dnmt3a[,1],vaf_tet2_dnmt3a[,2],pch=20,,col=rgb(0.5,0,0.5, alpha=0.5))
segments(0,0,0,m); segments(m,m,0,m); segments(0,0,m,0); segments(m,m,m,0)
segments(0.5,0,0,0.5,col="gray",lty="dashed")
lat = (0:6)*0.1; for(i in 1:length(lat)){segments(lat[i],0,lat[i],-0.01)}
text(-0.04,lat,lat,cex=0.8); text(lat,-0.04,lat,cex=0.8)
lat = (0:6)*0.1; for(i in 1:length(lat)){segments(0,lat[i],-0.01,lat[i])}
lat = (1:6)*0.1; for(i in 1:length(lat)){segments(lat[i]-0.05,0,lat[i]-0.05,-0.005)}
lat = (1:6)*0.1; for(i in 1:length(lat)){segments(0,lat[i]-0.05,-0.005,lat[i]-0.05)}
text(-0.1,0.3,"VAF of DNMT3A mutation",srt=90,cex=0.8)
text(0.3,-0.1,"VAF of TET2 mutation",cex=0.8)
dev.off()

pdf("tet2_asxl1.pdf")
m = 0.6
plot(NULL,NULL,xlim=c(-0.15,m+0.15),ylim=c(-0.15,m+0.15),axes=FALSE,xlab="",ylab="")
points(vaf_tet2_asxl1[,1],vaf_tet2_asxl1[,2],pch=20,,col=rgb(0.5,0,0.5, alpha=0.5))
segments(0,0,0,m); segments(m,m,0,m); segments(0,0,m,0); segments(m,m,m,0)
segments(0.5,0,0,0.5,col="gray",lty="dashed")
lat = (0:6)*0.1; for(i in 1:length(lat)){segments(lat[i],0,lat[i],-0.01)}
text(-0.04,lat,lat,cex=0.8); text(lat,-0.04,lat,cex=0.8)
lat = (0:6)*0.1; for(i in 1:length(lat)){segments(0,lat[i],-0.01,lat[i])}
lat = (1:6)*0.1; for(i in 1:length(lat)){segments(lat[i]-0.05,0,lat[i]-0.05,-0.005)}
lat = (1:6)*0.1; for(i in 1:length(lat)){segments(0,lat[i]-0.05,-0.005,lat[i]-0.05)}
text(-0.1,0.3,"VAF of ASXL1 mutation",srt=90,cex=0.8)
text(0.3,-0.1,"VAF of TET2 mutation",cex=0.8)
dev.off()

pdf("tet2_srsf2.pdf")
m = 0.6
plot(NULL,NULL,xlim=c(-0.15,m+0.15),ylim=c(-0.15,m+0.15),axes=FALSE,xlab="",ylab="")
points(vaf_tet2_srsf2[,1],vaf_tet2_srsf2[,2],pch=20,,col=rgb(0.5,0,0.5, alpha=0.5))
segments(0,0,0,m); segments(m,m,0,m); segments(0,0,m,0); segments(m,m,m,0)
segments(0.5,0,0,0.5,col="gray",lty="dashed")
lat = (0:6)*0.1; for(i in 1:length(lat)){segments(lat[i],0,lat[i],-0.01)}
text(-0.04,lat,lat,cex=0.8); text(lat,-0.04,lat,cex=0.8)
lat = (0:6)*0.1; for(i in 1:length(lat)){segments(0,lat[i],-0.01,lat[i])}
lat = (1:6)*0.1; for(i in 1:length(lat)){segments(lat[i]-0.05,0,lat[i]-0.05,-0.005)}
lat = (1:6)*0.1; for(i in 1:length(lat)){segments(0,lat[i]-0.05,-0.005,lat[i]-0.05)}
text(-0.1,0.3,"VAF of SRSF2 mutation",srt=90,cex=0.8)
text(0.3,-0.1,"VAF of TET2 mutation",cex=0.8)
dev.off()

pdf("asxl1_srsf2.pdf")
m = 0.6
plot(NULL,NULL,xlim=c(-0.15,m+0.15),ylim=c(-0.15,m+0.15),axes=FALSE,xlab="",ylab="")
points(vaf_asxl1_srsf2[,1],vaf_asxl1_srsf2[,2],pch=20,,col=rgb(0.5,0,0.5, alpha=0.5))
segments(0,0,0,m); segments(m,m,0,m); segments(0,0,m,0); segments(m,m,m,0)
segments(0.5,0,0,0.5,col="gray",lty="dashed")
lat = (0:6)*0.1; for(i in 1:length(lat)){segments(lat[i],0,lat[i],-0.01)}
text(-0.04,lat,lat,cex=0.8); text(lat,-0.04,lat,cex=0.8)
lat = (0:6)*0.1; for(i in 1:length(lat)){segments(0,lat[i],-0.01,lat[i])}
lat = (1:6)*0.1; for(i in 1:length(lat)){segments(lat[i]-0.05,0,lat[i]-0.05,-0.005)}
lat = (1:6)*0.1; for(i in 1:length(lat)){segments(0,lat[i]-0.05,-0.005,lat[i]-0.05)}
text(-0.1,0.3,"VAF of ASXL1 mutation",srt=90,cex=0.8)
text(0.3,-0.1,"VAF of SRSF2 mutation",cex=0.8)
dev.off()

pdf("asxl1_cbl.pdf")
m = 0.6
plot(NULL,NULL,xlim=c(-0.15,m+0.15),ylim=c(-0.15,m+0.15),axes=FALSE,xlab="",ylab="")
points(vaf_asxl1_cbl[,1],vaf_asxl1_cbl[,2],pch=20,,col=rgb(0.5,0,0.5, alpha=0.5))
segments(0,0,0,m); segments(m,m,0,m); segments(0,0,m,0); segments(m,m,m,0)
segments(0.5,0,0,0.5,col="gray",lty="dashed")
lat = (0:6)*0.1; for(i in 1:length(lat)){segments(lat[i],0,lat[i],-0.01)}
text(-0.04,lat,lat,cex=0.8); text(lat,-0.04,lat,cex=0.8)
lat = (0:6)*0.1; for(i in 1:length(lat)){segments(0,lat[i],-0.01,lat[i])}
lat = (1:6)*0.1; for(i in 1:length(lat)){segments(lat[i]-0.05,0,lat[i]-0.05,-0.005)}
lat = (1:6)*0.1; for(i in 1:length(lat)){segments(0,lat[i]-0.05,-0.005,lat[i]-0.05)}
text(-0.1,0.3,"VAF of ASXL1 mutation",srt=90,cex=0.8)
text(0.3,-0.1,"VAF of CBL mutation",cex=0.8)
dev.off()


## End of analysis
q()
