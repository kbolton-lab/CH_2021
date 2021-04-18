
chr_TCRA = 14
stt_TCRA = 22.293663 - 1
end_TCRA = 23.019608 + 1

cna_file = ## import cna data on JGA
cna = read.table(cna_file,header=T,stringsAsFactor=F,sep="\t",quote="",colClass="character")

mut_file = ## import mutation data on JGA
mut = read.table(mut_file,header=T,stringsAsFactor=F,sep="\t",quote="",colClass="character",comment.char="")
mut = mut[which(mut$Gene.refGene=="TET2"),]

a = as.numeric(cna$start) - end_TCRA
b = as.numeric(cna$end) - stt_TCRA
c = (as.numeric(cna$chr) == chr_TCRA)
t = (cna$COPY_CHANGE == "loss")
sel = which(a*b <= 0 & c & t)

cna_sel = cna[sel,]
cna_sel = cna_sel[order(as.numeric(cna_sel$end),decreasing=F),]
cna_sel = cna_sel[order(as.numeric(cna_sel$start),decreasing=F),]


stt_cna = as.numeric(cna_sel$start)
end_cna = as.numeric(cna_sel$end)
id_cna = cna_sel$ID
d_cna = as.numeric(cna_sel$CELL_FRAC)

win1 = min(stt_cna,stt_TCRA) - (max(end_cna)-min(stt_cna))/2
win2 = max(end_cna,end_TCRA) + (max(end_cna)-min(stt_cna))/2

y_max = 5
y_min = -30
r = (win2-win1)/(y_max-y_min) #* 2


pdf("TRAC.pdf",width=10,height=10)

plot(
NULL,
NULL,
xlim=c(win1,win2),
ylim=c(y_min,y_max),
axes=FALSE,
xlab="",
ylab=""
)

text(stt_TCRA-0.2,1,"TCRA",font=4,adj=1)
polygon(c(stt_TCRA,stt_TCRA,end_TCRA,end_TCRA),c(1-1/5,1+1/5,1+1/5,1-1/5),col="gray",border="transparent")

segments(win1,3,win2,3,lwd=2.5)
lat1 = c(21,22,23,24)
for(i in 1:length(lat1)){
	segments(lat1[i],3,lat1[i],3+0.2,lwd=2.5)
	text(lat1[i],4,paste0(lat1[i]," Mb"),adj=0.5)
}

for(i in 1:length(stt_cna)){

	dense = d_cna[i]
	cp_func = colorRampPalette(c("white", "dodgerblue3"))
	cp = cp_func(100)
	col_2 = cp[(log10(dense*100)+2)*100/4]
	l_width = 5
	segments(stt_cna[i],-i,end_cna[i],-i,col=col_2,lwd=l_width)
		
		m = which(mut$id == id_cna[i])
		if(length(m) == 0){ next }
			
			if(length(m) >= 2) { in_col = "brown" }

			cur_type = mut$Merge_Func[m]
			if(length(m) >= 2){
				in_col = "brown"
			}else{
				if(cur_type == "nonsynonymous SNV" | cur_type == "stoploss"){ # 4: missense
					in_col = "deepskyblue1"
				}else if(cur_type == "splicing"){ # 5: splice site
				in_col = "purple"
				}else if(cur_type == "frameshift deletion" | cur_type == "frameshift insertion"){ # 6: frameshift indel
				in_col = "salmon"
				}else if(cur_type == "stopgain"){ # 7: stopgain
				in_col = "red"
				}else if(cur_type == "nonframeshift deletion" | cur_type == "nonframeshift insertion"){ # 8: inframe indel
				in_col = "yellow"
				}
			}

			x_c = stt_cna[i] - 0.2
			y_c = -i
			theta <- seq(-pi, pi, length=100)
			polygon(x_c+0.2*cos(theta)*r,y_c+0.2*sin(theta),col="red",border="transparent")

}

dev.off()

