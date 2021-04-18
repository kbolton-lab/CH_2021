
mut_path = "TCGA_SNV.txt"
cna_path = "TCGA_CNA.txt"
mut_all = read.table(mut_path,header=T,stringsAsFactor=F,sep="\t",quote="",colClass="character")
cna_all = read.table(cna_path,header=T,stringsAsFactor=F,sep="\t",quote="",colClass="character")
cna_all$SAMPLE = gsub(".-...-....-...CEL$","",cna_all$SAMPLE)

id_path = "TCGA_id.txt"
id = unlist(read.table(id_path,header=F,stringsAsFactor=F,sep="\t",quote="",colClass="character"))

basic_path = "TCGA_age.txt"
basic = read.table(basic_path,header=T,stringsAsFactor=F,sep="\t",quote="",colClass="character")
rownames(basic) = basic[,1]
age = basic[id,2]

is_mut = as.numeric((is.element(id,mut_all$id)))
is_cna = as.numeric((is.element(id,cna_all$SAMPLE)))
is_chip = as.numeric(is_mut + is_cna != 0)

age_category = age
age_category[which(age_category < 30)] = 1
age_category[which(age_category >= 30 & age_category < 40)] = 2
age_category[which(age_category >= 40 & age_category < 50)] = 3
age_category[which(age_category >= 50 & age_category < 60)] = 4
age_category[which(age_category >= 60 & age_category < 70)] = 5
age_category[which(age_category >= 70 & age_category < 80)] = 6
age_category[which(age_category >= 80 & age_category <= 90)] = 7

age_num = age_mut = age_cna = numeric(7)
for (i in 1:7){
	age_num[i] = sum(na.omit(age_category == i))
	age_mut[i] = sum(na.omit(is_mut[which(age_category == i)]))
	age_cna[i] = sum(na.omit(is_cna[which(age_category == i)]))
}

age_prop_mut = age_mut / age_num
age_prop_cna = age_cna / age_num

age_prop_mut_lw = age_prop_mut_up = age_prop_cna_lw = age_prop_cna_up = numeric(7)

for (i in 1:7){
	age_prop_mut_lw[i] = binom.test(age_mut[i],age_num[i])$conf.int[1]
	age_prop_mut_up[i] = binom.test(age_mut[i],age_num[i])$conf.int[2]
	age_prop_cna_lw[i] = binom.test(age_cna[i],age_num[i])$conf.int[1]
	age_prop_cna_up[i] = binom.test(age_cna[i],age_num[i])$conf.int[2]
}


## mut and cna (iuncluding male and femal) ##
pdf("TCGA_age_distribution.pdf",width=10,height=10)

col_mut = rgb(1,0,0,alpha=0.65) #"mediumpurple4"
col_cna = rgb(0,0,1,alpha=0.65) # "springgreen2"
col_mut_trans = rgb(1,0,0,alpha=0.25) #"#9370DB35" # #9370DB + 20
col_cna_trans = rgb(0,0,1,alpha=0.25) #"#00FF7F35" # #00FF7F + 20

x = 1:7
plot(
x,age_prop_mut,
xlim=c(-0.5,8),
ylim=c(-0.15,0.4),
type="l",
xlab="",ylab="",
xaxt="n",yaxt="n",
bty = "n",
col=col_mut,
lwd=2
)

polygon(
c(x,rev(x)),
c(age_prop_mut_lw,rev(age_prop_mut_up)),
col=col_mut_trans,
border="transparent")

par(new=T)
plot(
x,age_prop_cna,
xlim=c(-0.5,8),
ylim=c(-0.15,0.4),
type="l",
col=col_cna,
xlab="",ylab="",
xaxt="n",yaxt="n",
bty = "n",
lwd=2
)

polygon(c(x,rev(x)),
c(age_prop_cna_lw,rev(age_prop_cna_up)),
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

segments(1.5,0.60,2,0.60,col=col_mut,lwd=2)
polygon(c(1.5,2,2,1.5),c(0.59,0.59,0.61,0.61),col=col_mut_trans,border="transparent")
text(2.2,0.60,"Mutations",adj=0,cex=1.5,font=2)
segments(1.5,0.55,2,0.55,col=col_cna,lwd=2)
polygon(c(1.5,2,2,1.5),c(0.54,0.54,0.56,0.56),col=col_cna_trans,border="transparent")
text(2.2,0.55,"CNAs",adj=0,cex=1.5,font=2)

dev.off()

