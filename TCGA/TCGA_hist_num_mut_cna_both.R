# library
library(tidyverse)

id_path = "TCGA_id.txt"
mut_path = "TCGA_SNV.txt"
cna_path = "TCGA_CNA.txt"

id = unlist(read.table(id_path,header=F,stringsAsFactor=F,quote="",sep="\t",comment.char=""))

mut = read.table(mut_path,header=T,stringsAsFactor=F,quote="",sep="\t",comment.char="")
mut_count = rep(0,length(id)); names(mut_count) = id

mut = mut[which(is.element(mut$id,id)),]


for(i in 1:nrow(mut)){
	cur_id = mut[i,1]
	cur_gene = mut[i,8]
	if(! is.element(cur_id,id)){
		next
	}	
	mut_count[cur_id] = mut_count[cur_id] + 1
}

cna = read.table(cna_path,header=T,stringsAsFactor=F,quote="",sep="\t",comment.char="")
cna_count = rep(0,length(id)); names(cna_count) = id
cna$SAMPLE = gsub(".-...-....-...CEL$","",cna$SAMPLE)

for(i in 1:nrow(cna)){
	cur_id = cna$SAMPLE[i]
	if(! is.element(cur_id,id)){
		next
	}
	cna_count[cur_id] = cna_count[cur_id] + 1
}

total_count = mut_count + cna_count
up = 8
total_count[which(total_count>=up)] = up


#############
pdf("hist_num_mut_cna_both.pdf",width=15,height=10)

plot(NULL, NULL,
xlim=c(-5,15),ylim=c(-300,700),
axes=FALSE,
xlab="", ylab="",
main=NULL
)

for(i in 1:max(total_count)){
	total = sum(total_count  == i)
	none = sum(total_count  == i & mut_count == 0 & cna_count == 0)
	mut_only = sum(total_count  == i & mut_count >= 1 & cna_count == 0)
	cna_only = sum(total_count  == i & mut_count == 0 & cna_count >= 1)
	both = sum(total_count  == i & mut_count >= 1 & cna_count >= 1)

#	print(total == none + mut_only + cna_only + both)

	wid=0.5
	rect(i-wid/2,0,i+wid/2,none,col="gray",border="transparent")
	rect(i-wid/2,none,i+wid/2,none+mut_only,col=rgb(1,0,0,alpha=0.65),border="transparent")
	rect(i-wid/2,none+mut_only,i+wid/2,none+mut_only+cna_only,col=rgb(0,0,1,alpha=0.65),border="transparent")
	rect(i-wid/2,none+mut_only+cna_only,i+wid/2,none+mut_only+cna_only+both,col=rgb(1,0,1,alpha=0.65),border="transparent")
	text(i,none+mut_only+cna_only+both+30,none+mut_only+cna_only+both,font=2)

	text(i,-50,i,cex=1,font=2,adj=0.5)
}
text(4.5,-500,"Number of alterations",cex=1.3,font=2,adj=0.5)


segments(0,0,max(total_count)+1,0)
segments(0,0,0,600)
lat_y1 = 100*(0:6)
for(i in 1:length(lat_y1)){
	segments(0,lat_y1[i],-0.2,lat_y1[i])
	text(-0.4,lat_y1[i],lat_y1[i],adj=1,cex=1,font=2)
}
lat_y2 = setdiff(20*(1:30),lat_y1)
for(i in 1:length(lat_y2)){   
	segments(0,lat_y2[i],-0.1,lat_y2[i])
}
text(-2,300,"Number of subjects",adj=0.5,cex=1.3,font=2,srt=90)


dx = 0  
dy = 350
mag = 50

for(i in 4:max(total_count)){
	total = mag*sum(total_count  == i)
	none = mag*sum(total_count  == i & mut_count == 0 & cna_count == 0)
	mut_only = mag*sum(total_count  == i & mut_count >= 1 & cna_count == 0)
	cna_only = mag*sum(total_count  == i & mut_count == 0 & cna_count >= 1)
	both = mag*sum(total_count  == i & mut_count >= 1 & cna_count >= 1)

#	print(total == none + mut_only + cna_only + both)
	wid = 0.5
	rect(i-wid/2+dx,0+dy,i+wid/2+dx,none+dy,col="gray",border="transparent")
	rect(i-wid/2+dx,none+dy,i+wid/2+dx,none+mut_only+dy,col=rgb(1,0,0,alpha=0.65),border="transparent")
	rect(i-wid/2+dx,none+mut_only+dy,i+wid/2+dx,none+mut_only+cna_only+dy,col=rgb(0,0,1,alpha=0.65),border="transparent")
	rect(i-wid/2+dx,none+mut_only+cna_only+dy,i+wid/2+dx,none+mut_only+cna_only+both+dy,col=rgb(1,0,1,alpha=0.65),border="transparent")
#	text(i+dx,none+mut_only+cna_only+both+dy,(none+mut_only+cna_only+both)/mag)

	text(i+dx,-50+dy,i,cex=1,font=2,adj=0.5)
}
# text(6.5,-400+dy,"Number of alterations",cex=1,font=2,adj=0.5)

dx = 4.5
dy = 350

segments(-1+dx,0+dy,8,0+dy)
segments(-1+dx,0+dy,-1+dx,5*mag+dy)

lat_y2 = setdiff((0:5),lat_y1)
for(i in 1:length(lat_y2)){
	segments(-1+dx,lat_y2[i]*mag+dy,-1.1+dx,lat_y2[i]*mag+dy)
}
# text(-2.3+dx,mean(lat_y1)*mag+dy,"Number of subjects",adj=0.5,cex=1,font=2,srt=90)


rect(10,200,10.5,300,col=rgb(1,0,0,alpha=0.65),border="transparent"); text(11,250,"Mutation alone",adj=0,font=2,cex=1)
rect(10,400,10.5,500,col=rgb(0,0,1,alpha=0.65),border="transparent"); text(11,450,"CNA alone",adj=0,font=2,cex=1)
rect(10,600,10.5,700,col=rgb(1,0,1,alpha=0.65),border="transparent"); text(11,650,"Mutation + CNA",adj=0,font=2,cex=1)


segments(4,100,7,100)#,col="gray")
segments(4,75,4,100)#,col="gray")
segments(7,75,7,100)#,col="gray")
segments(5.5,100,5.5,250)#,col="gray")
polygon(c(5.2,5.8,5.5),c(225,225,250),col="black")#col="gray")
#arrows(6.5,350,6.5,1050,length=0.2,angle=55)
dev.off()





