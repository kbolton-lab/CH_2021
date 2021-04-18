
## draw Fig.1a ##
mut_path = ## please access mutation file on JGA
cna_path = ## please access cna file on JGA
id = as.character(1:11234)

## count mutation ##
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

## count cna ##
cna = read.table(cna_path,header=T,stringsAsFactor=F,quote="",sep="\t",comment.char="")
cna_count = rep(0,length(id)); names(cna_count) = id
cna = cna[which(is.element(cna$IDY,id)),]
for(i in 1:nrow(cna)){
	cur_id = cna$IDY[i]
	if(! is.element(cur_id,id)){
		next
	}
	cna_count[cur_id] = cna_count[cur_id] + 1
}


## merge mut and cna ##
total_count = mut_count + cna_count
up = 8
total_count[which(total_count>=up)] = up


## print out stats ##
print("#mut")
sum(mut_count)

print("#subjects with mut")
sum(mut_count >=1)

print("#subjects with mut(1,2,>=3)")
sum(mut_count ==1)
sum(mut_count ==2)
sum(mut_count >=3)

print("#cna")
sum(cna_count)

print("#subjects with cna")
sum(cna_count >=1)

print("#subjects with multiple cnas")
sum(cna_count >=2)

print("#subjects with at least one mut and/or CNA")
sum(mut_count >=1 | cna_count >= 1)

print("#subjects with both mut and cna")
sum(mut_count >=1 & cna_count >= 1)

print("#subjects with multiple alterations judged by mut or cna")
sum(mut_count >=2 | cna_count >= 2)

print("#subjects with multiple alterations judged by mut + cna")
sum(total_count >= 2)

a = sum(mut_count >=1 & cna_count >= 1)
b = sum(mut_count >=1) - sum(mut_count >=1 & cna_count >= 1)
c = sum(cna_count >=1) - sum(mut_count >=1 & cna_count >= 1)
d = sum(mut_count == 0 & cna_count == 0)
fisher.test(matrix(c(a,b,c,d),ncol=2))


## draw figure ##
pdf("hist_num_mut_cna_both.pdf",width=15,height=10)

plot(NULL, NULL,
xlim=c(-5,15),ylim=c(-1000,3000+1000),
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

	wid=0.5
	rect(i-wid/2,0,i+wid/2,none,col="gray",border="transparent")
	rect(i-wid/2,none,i+wid/2,none+mut_only,col=rgb(1,0,0,alpha=0.65),border="transparent")
	rect(i-wid/2,none+mut_only,i+wid/2,none+mut_only+cna_only,col=rgb(0,0,1,alpha=0.65),border="transparent")
	rect(i-wid/2,none+mut_only+cna_only,i+wid/2,none+mut_only+cna_only+both,col=rgb(1,0,1,alpha=0.65),border="transparent")
	text(i,none+mut_only+cna_only+both+100,none+mut_only+cna_only+both,font=2)

	text(i,-200,i,cex=1,font=2,adj=0.5)
}
text(4.5,-500,"Number of alterations",cex=1.3,font=2,adj=0.5)


segments(0,0,max(total_count)+1,0)
segments(0,0,0,3200)
lat_y1 = 1000*(0:3)
for(i in 1:length(lat_y1)){
	segments(0,lat_y1[i],-0.2,lat_y1[i])
	text(-0.4,lat_y1[i],lat_y1[i],adj=1,cex=1,font=2)
}
lat_y2 = setdiff(200*(1:16),lat_y1)
for(i in 1:length(lat_y2)){   
	segments(0,lat_y2[i],-0.1,lat_y2[i])
}
text(-2,1500,"Number of subjects",adj=0.5,cex=1.3,font=2,srt=90)


dx = 0  
dy = 1650
mag = 50

for(i in 5:max(total_count)){
	total = mag*sum(total_count  == i)
	none = mag*sum(total_count  == i & mut_count == 0 & cna_count == 0)
	mut_only = mag*sum(total_count  == i & mut_count >= 1 & cna_count == 0)
	cna_only = mag*sum(total_count  == i & mut_count == 0 & cna_count >= 1)
	both = mag*sum(total_count  == i & mut_count >= 1 & cna_count >= 1)

	wid = 0.5
	rect(i-wid/2+dx,0+dy,i+wid/2+dx,none+dy,col="gray",border="transparent")
	rect(i-wid/2+dx,none+dy,i+wid/2+dx,none+mut_only+dy,col=rgb(1,0,0,alpha=0.65),border="transparent")
	rect(i-wid/2+dx,none+mut_only+dy,i+wid/2+dx,none+mut_only+cna_only+dy,col=rgb(0,0,1,alpha=0.65),border="transparent")
	rect(i-wid/2+dx,none+mut_only+cna_only+dy,i+wid/2+dx,none+mut_only+cna_only+both+dy,col=rgb(1,0,1,alpha=0.65),border="transparent")

	text(i+dx,-200+dy,i,cex=1,font=2,adj=0.5)
}

dx = 5.5
dy = 1650

segments(-1+dx,0+dy,max(total_count)+1,0+dy)
segments(-1+dx,0+dy,-1+dx,40*mag+dy)
lat_y1 = 5*(0:8)
for(i in 1:length(lat_y1)){
	segments(-1+dx,lat_y1[i]*mag+dy,-1.2+dx,lat_y1[i]*mag+dy)
	text(-1.4+dx,lat_y1[i]*mag+dy,lat_y1[i],adj=1,cex=1,font=2)
}
lat_y2 = setdiff((0:40),lat_y1)
for(i in 1:length(lat_y2)){
	segments(-1+dx,lat_y2[i]*mag+dy,-1.1+dx,lat_y2[i]*mag+dy)
}


rect(10,200,10.5,300,col=rgb(1,0,0,alpha=0.65),border="transparent"); text(11,250,"Mutation alone",adj=0,font=2,cex=1)
rect(10,400,10.5,500,col=rgb(0,0,1,alpha=0.65),border="transparent"); text(11,450,"CNA alone",adj=0,font=2,cex=1)
rect(10,600,10.5,700,col=rgb(1,0,1,alpha=0.65),border="transparent"); text(11,650,"Mutation + CNA",adj=0,font=2,cex=1)


segments(4.5,450,8.5,450)#,col="gray")
segments(4.5,350,4.5,450)#,col="gray")
segments(8.5,350,8.5,450)#,col="gray")
segments(6.5,450,6.5,1100)#,col="gray")
polygon(c(6.2,6.8,6.5),c(1100,1100,1200),col="black")#col="gray")
#arrows(6.5,350,6.5,1050,length=0.2,angle=55)
dev.off()


## End of analysis
q()



