## Script for Supplementary Fig.1d ##

# import data
file = ~~ ## please access cna file on JGA
dat = read.table(file,header=T,sep="\t",quote="",comment.char="",stringsAsFactor=F)
dat = dat[which(dat$COPY_CHANGE != "unknown"),]

# extract data
size = dat$end - dat$start
cf = dat$CELL_FRAC
type = dat$COPY_CHANGE

# define color
cols = rep("gray",length(type))
cols[which(type=="neutral")] = "chartreuse3"
cols[which(type=="loss")] = "dodgerblue3"
cols[which(type=="gain")] = "firebrick3"

# plot
pdf("figure_cna_cf.pdf",width=7,height=7)

# plot cna data
plot(
size,cf,
log="xy",
pch=4,
cex=0.4,
xlab="",
ylab="",
xaxt="n",
yaxt="n",
bty="n",
xlim=c(0.001,500),
ylim=c(0.0005,1),
col=cols
)

# scale x (event size)
lat_x_small = c(0.01*(1:9),0.1*(1:9),1*(1:9),10*(1:9),100*(1:5))
lat_x_tall = c(0.05,0.1,0.5,1,5,10,50,100,500)
segments(0.01,0.001,500,0.001)
for(i in 1:length(lat_x_small)){
	segments(lat_x_small[i],0.001,lat_x_small[i],0.001*1.025)
}
for(i in 1:length(lat_x_tall)){
	segments(lat_x_tall[i],0.001,lat_x_tall[i],0.001*1.075)
}
segments(0.01,1,500,1)
for(i in 1:length(lat_x_small)){
	segments(lat_x_small[i],1,lat_x_small[i],1/1.025)
}
for(i in 1:length(lat_x_tall)){
	segments(lat_x_tall[i],1,lat_x_tall[i],1/1.075)
}

for(i in 1:length(lat_x_tall)){
	text(lat_x_tall[i],0.001/1.15,lat_x_tall[i],cex=0.6)
}
text(1.7,0.001/2,"Event size (Mb)",cex=0.8,adj=0.5)

# scale y (cell fraction)
lat_y_small = c(0.001*(1:9),0.01*(1:9),0.1*(1:9),1)
lat_y_tall = c(0.005,0.01,0.05,0.1,0.5,1)
segments(0.01,0.001,0.01,1)
for(i in 1:length(lat_y_small)){ 
	segments(0.01,lat_y_small[i],0.01*1.025,lat_y_small[i])
}
for(i in 1:length(lat_y_tall)){
	segments(0.01,lat_y_tall[i],0.01*1.075,lat_y_tall[i])
}
segments(500,0.001,500,1)
for(i in 1:length(lat_y_small)){
	segments(500,lat_y_small[i],500/1.025,lat_y_small[i])
}
for(i in 1:length(lat_y_tall)){
	segments(500,lat_y_tall[i],500/1.075,lat_y_tall[i]) 
}

for(i in 1:length(lat_y_tall)){
	text(0.01/1.15,lat_y_tall[i],lat_y_tall[i],cex=0.6,adj=1)
}
text(0.01/4,0.04,"Cell fraction",cex=0.8,adj=0.5,srt=90)

# legend
text(0.1,0.005*1.4,"CNN-LOH",cex=0.6,adj=0)
points(0.07,0.005*1.4,pch=4,cex=0.4,col="chartreuse3")
text(0.1,0.005,"Deletion",cex=0.6,adj=0)
points(0.07,0.005,pch=4,cex=0.4,col="dodgerblue3") 
text(0.1,0.005/1.4,"Duplication",cex=0.6,adj=0)
points(0.07,0.005/1.4,pch=4,cex=0.4,col="firebrick3") 

# end
dev.off()

