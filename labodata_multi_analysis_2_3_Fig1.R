
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


## get largerst clone size
cs_mut = as.numeric(apply(mat_mut,1,max))
cs_cna = as.numeric(apply(mat_cna,1,max))
all_size = cs_mut*2*(cs_mut*2>cs_cna) + cs_cna*(cs_mut*2<=cs_cna)
all_size[which(all_size>1)] = 1


## get total number of alterations
all_n_alt = n_mut + n_cna



### cutoff on CBC ##

## please access clinical information provided by BBJ ##
Hb = ~~~
WBC = ~~
Plt = ~~~

Hb_l_thld = 11
Hb_h_m_thld = 16.5
Hb_h_f_thld = 16
Ht_h_thld = 49
WBC_l_thld = 3000
WBC_h_thld = 10000
Plt_l_thld = 10
Plt_h_thld = 45

Hb_l = as.numeric(all_dat$Hb <= Hb_l_thld)
WBC_l = as.numeric(all_dat$WBC <= WBC_l_thld)
Plt_l = as.numeric(all_dat$PLT<= Plt_l_thld)

Ht_h = as.numeric(all_dat$Ht >= Ht_h_thld)
Hb_h = as.numeric((all_dat$Hb >= Hb_h_m_thld & all_sex=="1") | (all_dat$Hb >= Hb_h_f_thld & all_sex=="2")) 
WBC_h = as.numeric(all_dat$WBC >= WBC_h_thld)
Plt_h = as.numeric(all_dat$PLT >= Plt_h_thld)

cbc_abnormal = as.numeric(Hb_l == 1 | WBC_l == 1 | Plt_l == 1 | Ht_h == 1 | Hb_h == 1 | WBC_h == 1 | Plt_h == 1)
cytopenia_all = as.numeric(Hb_l == 1 | WBC_l == 1 | Plt_l == 1)
cytopenia_m = as.numeric(Hb_l + WBC_l + Plt_l >= 2)



## CBC and number of alterations ##
## Draw Fig.1c ##

cbc_abnormal.n = cbc_normal.n = cytopenia_all.n = 
cytopenia_1.n = cytopenia_2.n = cytopenia_3.n =  cytopenia_m.n = 
Hb_l.n = WBC_l.n = Plt_l.n = Ht_h.n = Hb_h.n = WBC_h.n = Plt_h.n = rep(0,5)

for(i in 1:length(cbc_abnormal.n)){
	if(i!=length(cbc_abnormal.n)){
		cbc_abnormal.n[i] = sum(all_n_alt[which(cbc_abnormal==1)]==i-1)
		cbc_normal.n[i] = sum(all_n_alt[which(cbc_abnormal==0)]==i-1)
		cytopenia_all.n[i] = sum(all_n_alt[which(cytopenia_all==1)]==i-1)
		cytopenia_1.n[i] = sum(all_n_alt[which(cytopenia_1==1)]==i-1)
		cytopenia_2.n[i] = sum(all_n_alt[which(cytopenia_2==1)]==i-1)
		cytopenia_3.n[i] = sum(all_n_alt[which(cytopenia_3==1)]==i-1)
		cytopenia_m.n[i] = sum(all_n_alt[which(cytopenia_m==1)]==i-1)
		Hb_l.n[i] = sum(all_n_alt[which(Hb_l==1)]==i-1)
		WBC_l.n[i] = sum(all_n_alt[which(WBC_l==1)]==i-1)
		Plt_l.n[i] = sum(all_n_alt[which(Plt_l==1)]==i-1)
		Ht_h.n[i] = sum(all_n_alt[which(Ht_h==1)]==i-1)
		Hb_h.n[i] = sum(all_n_alt[which(Hb_h==1)]==i-1)
		WBC_h.n[i] = sum(all_n_alt[which(WBC_h==1)]==i-1)
		Plt_h.n[i] = sum(all_n_alt[which(Plt_h==1)]==i-1)
	}else{
		cbc_abnormal.n[i] = sum(all_n_alt[which(cbc_abnormal==1)]>=i-1)
		cbc_normal.n[i] = sum(all_n_alt[which(cbc_abnormal==0)]>=i-1)
		cytopenia_all.n[i] = sum(all_n_alt[which(cytopenia_all==1)]>=i-1)
		cytopenia_1.n[i] = sum(all_n_alt[which(cytopenia_1==1)]>=i-1)
		cytopenia_2.n[i] = sum(all_n_alt[which(cytopenia_2==1)]>=i-1)
		cytopenia_3.n[i] =  sum(all_n_alt[which(cytopenia_3==1)]>=i-1)
		cytopenia_m.n[i] =  sum(all_n_alt[which(cytopenia_m==1)]>=i-1)
		Hb_l.n[i] = sum(all_n_alt[which(Hb_l==1)]>=i-1)
		WBC_l.n[i] = sum(all_n_alt[which(WBC_l==1)]>=i-1)
		Plt_l.n[i] = sum(all_n_alt[which(Plt_l==1)]>=i-1)
		Ht_h.n[i] = sum(all_n_alt[which(Ht_h==1)]>=i-1)
		Hb_h.n[i] = sum(all_n_alt[which(Hb_h==1)]>=i-1)
		WBC_h.n[i] = sum(all_n_alt[which(WBC_h==1)]>=i-1)
		Plt_h.n[i] = sum(all_n_alt[which(Plt_h==1)]>=i-1)
	}
}

freq = function(x.n){sum(x.n[-1])/sum(x.n)}
cbc_abnormal.n.prop = freq(cbc_abnormal.n)
cbc_normal.n.prop = freq(cbc_normal.n)
cytopenia_all.n.prop = freq(cytopenia_all.n)
cytopenia_1.n.prop = freq(cytopenia_1.n)
cytopenia_2.n.prop = freq(cytopenia_2.n)
cytopenia_3.n.prop = freq(cytopenia_3.n)
cytopenia_m.n.prop = freq(cytopenia_m.n)
Hb_l.n.prop = freq(Hb_l.n)
WBC_l.n.prop = freq(WBC_l.n)
Plt_l.n.prop = freq(Plt_l.n)
Hb_h.n.prop = freq(Hb_h.n)
WBC_h.n.prop = freq(WBC_h.n)
Plt_h.n.prop = freq(Plt_h.n)

## CBC normal vs. CBC abnormal
wilcox.exact(all_n_alt[which(cbc_abnormal==1)],all_n_alt[which(cbc_abnormal==0)])$p.value
## CBC. normal vs. cytopenia
wilcox.exact(all_n_alt[which(cytopenia_all==1)],all_n_alt[which(cbc_abnormal==0)])$p.value
## CBC normal vs. multi-lineage cytopenia
wilcox.exact(all_n_alt[which(cytopenia_m==1)],all_n_alt[which(cbc_abnormal==0)])$p.value


bar = function(d_l,c){
plot(NULL,NULL,axes=F,xlab="",ylab="",
xlim=c(-1,length(d_l)+2),ylim=c(-0.5,1+0.5)
)
N = length(d_l)
for(i in 1:length(d_l)){
d = d_l[[i]]; d = d/sum(d)
rect(
(1:length(d))-0.3+0.6*(i-1)/N,
rep(0,length(d)),
(1:length(d))-0.3++0.6*i/N,
rep(d,length(d)),
border="transparent",col=c[i])
}

lat1 = (0:10)/10
lat2 = (1:100)/100
segments(0,0,0,1)
segments(0,lat1,-0.05,lat1)
segments(0,lat2,-0.025,lat2)
text(-0.1,lat1,lat1,adj=1,cex=0.5,font=2)
text(1:5,-.05,0:4,adj=.5,cex=0.7)
}

pdf("cbc_abnormal.n_all.pdf")
bar(list(cbc_normal.n,cbc_abnormal.n,cytopenia_all.n,cytopenia_m.n),
c("gray",adjustcolor("firebrick3",alpha=0.05),adjustcolor("firebrick3",alpha=0.2),
adjustcolor("firebrick3",alpha=0.9)))
dev.off()



## CBC abnormality and clone size ##
## Draw Fig.1d ##

pdf("size_cbc_plot.pdf",width=8,height=8)
plot(NULL,NULL,xlim=c(-1.5,5.5),ylim=c(log10(0.0005),1),axes=F,xlab="",ylab="")

cs1 = all_size[which(cbc_abnormal==0)][which(all_size[which(cbc_abnormal==0)]>1e-05)]
cs2 = all_size[which(cbc_abnormal==1)][which(all_size[which(cbc_abnormal==1)]>1e-05)]
cs3 = all_size[which(cytopenia_all==1)][which(all_size[which(cytopenia_all==1)]>1e-05)]
cs4 = all_size[which(cytopenia_2==1 | cytopenia_3==1)][which(all_size[which(cytopenia_2==1 | cytopenia_3==1)]>1e-05)]
cs1[which(cs1>1)] = cs2[which(cs2>1)] = cs3[which(cs3>1)] = cs4[which(cs4>1)] = 1
cs_l = list(cs1,cs2,cs3,cs4)

for(i in 1:length(cs_l)){
cs = cs_l[[i]]
cs_qt = quantile(cs)
cs_lw = cs_qt[2] - 1.5*(cs_qt[4]-cs_qt[2])
cs_lw = max(min(cs[which(cs>=cs_lw)]),cs_lw)
cs_up = cs_qt[4] + 1.5*(cs_qt[4]-cs_qt[2])
cs_up = min(max(cs[which(cs<=cs_up)]),cs_up)
cs_x = runif(length(cs),i-.3,i+.3)

points(cs_x,log10(cs),pch=20,cex=0.5,col=rgb(1,0,1,alpha=0.3))
rect(i-0.2,log10(cs_qt[2]),i+0.2,log10(cs_qt[4]),lwd=2,border="black")
segments(i-0.2,log10(cs_qt[3]),i+0.2,log10(cs_qt[3]),lwd=2,col="black")
segments(i,log10(cs_up),i,log10(cs_qt[4]),lwd=2,col="black")
segments(i,log10(cs_lw),i,log10(cs_qt[2]),lwd=2,col="black")
segments(i-0.1,log10(cs_up),i+0.1,log10(cs_up),lwd=2,col="black")
segments(i-0.1,log10(cs_lw),i+0.1,log10(cs_lw),lwd=2,col="black")

}

lat1 = c(0.01,0.1,1)
lat2 = c(0.001*(2:9),0.01*(2:9),0.1*(2:9),1)
segments(-0.05,log10(lat1),0,log10(lat1))
segments(-0.025,log10(lat2),0,log10(lat2))
text(-0.1,log10(lat1),lat1,cex=0.5,font=2,adj=1)
segments(0,log10(0.002),0,log10(1))
segments(0,log10(0.002),4.5,log10(0.002))
text(-0.5,(log10(1)+log10(0.002))/2,"Maximum clone size",srt=90,adj=.5,font=2,cex=0.8)

segments(1,.15,2,.15); text(1.5,.2,"**",font=2,adj=.5,cex=1.2)
segments(1,.3,3,.3); text(2,.35,"**",font=2,adj=.5,cex=1.2)

text(1,log10(0.002)-0.1,"Normal",adj=1,srt=45,cex=0.7,font=2)
text(2,log10(0.002)-0.1,"Abnormal",adj=1,srt=45,cex=0.7,font=2)
text(3,log10(0.002)-0.1,"Cytopenia\n(All)",adj=1,srt=45,cex=0.7,font=2)
text(4,log10(0.002)-0.1,"Cytopenia\n(multi)",adj=1,srt=45,cex=0.7,font=2)

dev.off()

## normal CBC vs. abnormal CBC
wilcox.test(cs1,cs2)
## normal CBC vs. cytopenia
wilcox.test(cs1,cs3)
## normal CBC vs. multi-lineage cytopenia
wilcox.test(cs1,cs4)



