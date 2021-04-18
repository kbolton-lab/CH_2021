
mut_path = "mut_all_210105.txt"   # edit 
cna_path = "cna_all_210105.txt"    # edit  
mut_all = read.table(mut_path,header=T,stringsAsFactor=F,sep="\t",quote="",colClass="character")
cna_all = read.table(cna_path,header=T,stringsAsFactor=F,sep="\t",quote="",colClass="character")

id_path = "id_all_210105.txt"
id = unlist(read.table(id_path,header=F,stringsAsFactor=F,sep="\t",quote="",colClass="character"))

## define HM positive cases ##
## please import clinical data
HM_ids = ~~~
is_HM = as.numeric(is.element(id,HM_ids))

## age and sex ##
## please import clinical data
basic_path = 
age = basic[id,"age"]
sex = basic[id,"sex"]

is_mut = as.numeric((is.element(id,mut_all$id)))
is_cna = as.numeric((is.element(id,cna_all$IDY)))
is_chip = as.numeric(is_mut + is_cna != 0)

age_category = age
age_category[which(age_category >= 60 & age_category < 65)] = 1
age_category[which(age_category >= 65 & age_category < 70)] = 2
age_category[which(age_category >= 70 & age_category < 75)] = 3
age_category[which(age_category >= 75 & age_category < 80)] = 4
age_category[which(age_category >= 80 & age_category < 85)] = 5
age_category[which(age_category >= 85 & age_category < 90)] = 6
# age_category[which(age_category >= 90 & age_category < 95)] = 7
age_category[which(age_category >= 90)] = 7

age_num_HM = age_mut_HM = age_cna_HM = numeric(7)
age_num_Other = age_mut_Other = age_cna_Other = numeric(7)
age_tet2_HM = age_dnmt3a_HM = numeric(7)
age_tet2_Other = age_dnmt3a_Other = numeric(7)

for (i in 1:7){
	age_num_HM[i] = sum(age_category == i & is_HM == 1)
	age_num_Other[i] = sum(age_category == i & is_HM == 0)
	age_mut_HM[i] = sum(is_mut[which(age_category == i & is_HM == 1)])
	age_mut_Other[i] = sum(is_mut[which(age_category == i & is_HM == 0)])
	age_cna_HM[i] = sum(is_cna[which(age_category == i & is_HM == 1)])
	age_cna_Other[i] = sum(is_cna[which(age_category == i & is_HM == 0)])
}

age_prop_mut_HM = age_mut_HM / age_num_HM
age_prop_mut_Other = age_mut_Other / age_num_Other
age_prop_cna_HM = age_cna_HM / age_num_HM
age_prop_cna_Other = age_cna_Other / age_num_Other

age_prop_mut_HM_lw = age_prop_mut_HM_up = age_prop_cna_HM_lw = age_prop_cna_HM_up = numeric(7)
age_prop_mut_Other_lw = age_prop_mut_Other_up = age_prop_cna_Other_lw = age_prop_cna_Other_up = numeric(7)

for (i in 1:7){
	age_prop_mut_HM_lw[i] = binom.test(age_mut_HM[i],age_num_HM[i])$conf.int[1]
	age_prop_mut_HM_up[i] = binom.test(age_mut_HM[i],age_num_HM[i])$conf.int[2]
	age_prop_mut_Other_lw[i] = binom.test(age_mut_Other[i],age_num_Other[i])$conf.int[1]
	age_prop_mut_Other_up[i] = binom.test(age_mut_Other[i],age_num_Other[i])$conf.int[2]
	
	age_prop_cna_HM_lw[i] = binom.test(age_cna_HM[i],age_num_HM[i])$conf.int[1]
	age_prop_cna_HM_up[i] = binom.test(age_cna_HM[i],age_num_HM[i])$conf.int[2]
	age_prop_cna_Other_lw[i] = binom.test(age_cna_Other[i],age_num_Other[i])$conf.int[1]
	age_prop_cna_Other_up[i] = binom.test(age_cna_Other[i],age_num_Other[i])$conf.int[2]
}

## mut and cna (iuncluding male and femal) ##
pdf("Fig1b_age_distribution_hm_or_not.pdf",width=10,height=10)

col_mut = rgb(1,0,0,alpha=0.65) #"mediumpurple4"
col_cna = rgb(0,0,1,alpha=0.65) # "springgreen2"
col_mut_trans = rgb(1,0,0,alpha=0.25) #"#9370DB35" # #9370DB + 20
col_cna_trans = rgb(0,0,1,alpha=0.25) #"#00FF7F35" # #00FF7F + 20

x = 1:7
plot(
x,age_prop_mut_HM,
xlim=c(-0.5,8),
ylim=c(-0.15,0.8),
type="l",
xlab="",ylab="",
xaxt="n",yaxt="n",
bty = "n",
col=col_mut,
lwd=2
)


par(new=T)
plot(
x,age_prop_mut_Other,
xlim=c(-0.5,8),
ylim=c(-0.15,0.8),
type="l",
xlab="",ylab="",
xaxt="n",yaxt="n",
bty = "n",
col=col_mut,
lwd=2,lty="dashed"
)
polygon(
c(x,rev(x)),
c(age_prop_mut_Other_lw,rev(age_prop_mut_Other_up)),
col=col_mut_trans,
border="transparent")

par(new=T)
plot(
x,age_prop_cna_HM,
xlim=c(-0.5,8),
ylim=c(-0.15,0.8),
type="l",
col=col_cna,
xlab="",ylab="",
xaxt="n",yaxt="n",
bty = "n",
lwd=2
)

par(new=T)
plot(
x,age_prop_cna_Other,
xlim=c(-0.5,8),
ylim=c(-0.15,0.8),
type="l",
col=col_cna,
xlab="",ylab="",
xaxt="n",yaxt="n",
bty = "n",
lwd=2,lty="dashed"
)
polygon(c(x,rev(x)),
c(age_prop_cna_Other_lw,rev(age_prop_cna_Other_up)),
col=col_cna_trans,
border="transparent")

segments(0.5,0,7.5,0)
segments(0.5,0,0.5,0.7)

x_labels = c("60-64","65-69","70-74","75-89","80-84","85-90","90-")
for(i in 1:length(x_labels)){
	segments(i,0,i,-0.01)
	text(i-0.2,-0.06,x_labels[i],srt=45,cex=1.5,font=2)
}
text(4,-0.14,"Age (Years)",cex=1.5,font=2)

lat_y_high = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7)
lat_y_small = c(0.05,0.15,0.25,0.35,0.45,0.55,0.65)

for(i in 1:length(lat_y_high)){
	segments(0.5,lat_y_high[i],0.43,lat_y_high[i])
	text(0.35,lat_y_high[i],lat_y_high[i],adj=1,cex=1.5,font=2)
}
for(i in 1:length(lat_y_small)){
	segments(0.5,lat_y_small[i],0.46,lat_y_small[i])	
}
text(-0.45,0.35,"Frequency",srt=90,cex=1.5,font=2)

segments(1.5,0.65,2,0.65,col=col_mut,lwd=2,lty="dashed")
polygon(c(1.5,2,2,1.5),c(0.64,0.64,0.66,0.66),col=col_mut_trans,border="transparent")
text(2.2,0.6,"Mutations",adj=0,cex=1.5,font=2)
segments(1.5,0.60,2,0.60,col=col_cna,lwd=2,lty="dashed")
polygon(c(1.5,2,2,1.5),c(0.59,0.59,0.61,0.61),col=col_cna_trans,border="transparent")
text(2.2,0.60,"CNAs",adj=0,cex=1.5,font=2)

segments(1.5,0.55,2,0.55,col=col_mut,lwd=2)
text(2.2,0.55,"Mutations",adj=0,cex=1.5,font=2)
segments(1.5,0.5,2,0.5,col=col_cna,lwd=2)
text(2.2,0.5,"CNAs",adj=0,cex=1.5,font=2)


dev.off()


## End analysis
q()

