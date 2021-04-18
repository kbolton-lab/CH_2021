
writeForest = function(d,path.pdf){

pdf(path.pdf,width=10,height=10)

plot(
NULL,NULL,axes=F,
xlim=c(-5,5),ylim=c(0,60),
xlab="",ylab=""
)

segments(0,0,0,60,col="gray",lwd=1)
segments(-1,0,2,0,col="black",lwd=1)
lat = c((0.1)*(1:9),1:9,10*(1:9))
lat1 = c(0.1,1,10,100)
for(i in 1:length(lat)){
	segments(log10(lat[i]),0,log10(lat[i]),-0.2,lwd=1)
	segments(log10(lat1[i]),0,log10(lat1[i]),-0.4,lwd=1)
	text(log10(lat1[i]),-0.8,lat1[i],cex=0.5,font=2)
}

rx = 0.017*1.2; ry = 0.12*1.2; theta = seq(-pi,pi,length=100)
rx1 = rx2 = 0.035; ry1 = ry2 =0.07
for(i in 1:nrow(d)){
	text(-2,60-i,rownames(d)[i],adj=1,font=2,cex=0.4)
	segments(log10(d[i,"lower .95"]),60-i,log10(d[i,"upper .95"]),60-i,lwd=0.5)
	polygon(log10(d[i,"exp(coef)"])+(rx*cos(theta)),60-i+ry*sin(theta),col="violet",lwd=0.2)
}

dev.off()
}



writeForest2 = function(d,path.pdf){

pdf(path.pdf,width=10,height=10)

plot(
NULL,NULL,axes=F,
xlim=c(-2,2),ylim=c(0,60),
xlab="",ylab=""
)

segments(0,0,0,60,col="gray",lwd=1)
segments(-1,0,1,0,col="black",lwd=1)
lat1 = c(0.1,1,2,3,4,5,10)
lat = c(0.1*5,1+0.1*5,2+0.1*5,3+0.1*5,4+0.1*5)
for(i in 1:length(lat)){
	segments(log10(lat[i]),0,log10(lat[i]),-0.2,lwd=1)
	segments(log10(lat1[i]),0,log10(lat1[i]),-0.4,lwd=1)
	text(log10(lat1[i]),-0.8,lat1[i],cex=0.5,font=2)
}

rx = 0.017*1.2; ry = 0.12*1.2; theta = seq(-pi,pi,length=100)
rx1 = rx2 = 0.035; ry1 = ry2 =0.07
for(i in 1:nrow(d)){
	text(-2,60-i,rownames(d)[i],adj=1,font=2,cex=0.4)
	segments(log10(d[i,"lower .95"]),60-i,log10(d[i,"upper .95"]),60-i,lwd=0.5)
	polygon(log10(d[i,"exp(coef)"])+(rx*cos(theta)),60-i+ry*sin(theta),col="violet",lwd=0.2)
}

dev.off()
}




writeForest3 = function(d,d1,d2,path.pdf){

pdf(path.pdf,width=5,height=10)

plot(
NULL,NULL,axes=F,
xlim=c(-2,2),ylim=c(0,60),
xlab="",ylab=""
)

segments(0,0,0,60,col="gray",lwd=1)
segments(-1,0,2,0,col="black",lwd=1)
lat = c((0.1)*(1:9),1:9,10*(1:9))
lat1 = c(0.1,1,10,100)
for(i in 1:length(lat)){
	segments(log10(lat[i]),0,log10(lat[i]),-0.2,lwd=1)
	segments(log10(lat1[i]),0,log10(lat1[i]),-0.4,lwd=1)
	text(log10(lat1[i]),-0.8,lat1[i],cex=0.5,font=2)
}


rx = 0.017*1.2; ry = 0.12*1.2; theta = seq(-pi,pi,length=100)
#rx1 = rx2 = 0.035; ry1 = ry2 = 0.07
#rx1 = rx2 = 0.035; ry1 = ry2 = 0.07
for(i in 1:nrow(d)){
	size = log10(d[i,"n"])
	rx2 = 0.017*size; ry2 = 0.12*size; theta = seq(-pi,pi,length=100)
	polygon(2.05+(rx2*cos(theta)),60-i+ry2*sin(theta),col="violet",lwd=0.2)
	
	text(-2,60-i,rownames(d)[i],adj=1,font=2,cex=0.4)
	segments(log10(d[i,"lower_95"]),60-i,log10(d[i,"upper_95"]),60-i,lwd=0.5)
	polygon(log10(d[i,"exp_coef"])+(rx*cos(theta)),60-i+ry*sin(theta),col="violet",lwd=0.2)
	polygon(log10(d1[i,"exp_coef"])+(rx*cos(theta)),60-i+ry*sin(theta),col="red",lwd=0.2)
	polygon(log10(d2[i,"exp_coef"])+(rx*cos(theta)),60-i+ry*sin(theta),col="blue",lwd=0.2)
}

dev.off()
}


writeForestLegend = function(path.pdf){

pdf(path.pdf,width=5,height=10)
plot(NULL,NULL,axes=F,xlim=c(-2,2),ylim=c(0,60),xlab="",ylab="")

rx = 0.017*1.2; ry = 0.12*1.2; theta = seq(-pi,pi,length=100)
polygon(1+(rx*cos(theta)),60-10+ry*sin(theta),col="violet",lwd=0.2)
polygon(1+(rx*cos(theta)),60-11+ry*sin(theta),col="red",lwd=0.2)
polygon(1+(rx*cos(theta)),60-12+ry*sin(theta),col="blue",lwd=0.2)

for(i in 1:3){
	size = i
	rx2 = 0.017*size; ry2 = 0.12*size; theta = seq(-pi,pi,length=100)
	polygon(2.05+(rx2*cos(theta)),60-i+ry2*sin(theta),col="violet",lwd=0.2)
}

dev.off()
}




writeForest4 = function(d,d1,d2,path.pdf){

pdf(path.pdf,width=5,height=6)
par(mfrow=c(1,3))

plot(NULL,NULL,axes=F,xlim=c(-2,2),ylim=c(0,60),xlab="",ylab="")
segments(0,0,0,60,col="gray",lwd=1)
segments(-1,0,2,0,col="black",lwd=1)
lat = c((0.1)*(1:9),1:9,10*(1:9))
lat1 = c(0.1,1,10,100)
for(i in 1:length(lat)){
	segments(log10(lat[i]),0,log10(lat[i]),-0.2,lwd=1)
	segments(log10(lat1[i]),0,log10(lat1[i]),-0.4,lwd=1)
	text(log10(lat1[i]),-0.8,lat1[i],cex=0.5,font=2)
}
rx = 0.017*1.3*2.5; ry = 0.12*1.3; theta = seq(-pi,pi,length=100)
for(i in 1:nrow(d)){
#	size = log10(d[i,"n"])
#	rx2 = 0.017*size*5; ry2 = 0.12*size; theta = seq(-pi,pi,length=100)
#	polygon(2.05+(rx2*cos(theta)),60-i+ry2*sin(theta),col="violet",lwd=0.2)
	text(-1,60-i,d[i,"n"],adj=1,font=2,cex=0.4)
	text(-2,60-i,rownames(d)[i],adj=1,font=2,cex=0.4)
	segments(log10(d[i,"lower_95"]),60-i,log10(d[i,"upper_95"]),60-i,lwd=0.3)
	polygon(log10(d[i,"exp_coef"])+(rx*cos(theta)),60-i+ry*sin(theta),col=rgb(1,0,1,alpha=1),lwd=0.2)
}

plot(NULL,NULL,axes=F,xlim=c(-2,2),ylim=c(0,60),xlab="",ylab="")
segments(0,0,0,60,col="gray",lwd=1)
segments(-1,0,2,0,col="black",lwd=1)
lat = c((0.1)*(1:9),1:9,10*(1:9))
lat1 = c(0.1,1,10,100)
for(i in 1:length(lat)){
	segments(log10(lat[i]),0,log10(lat[i]),-0.2,lwd=1)
	segments(log10(lat1[i]),0,log10(lat1[i]),-0.4,lwd=1)
	text(log10(lat1[i]),-0.8,lat1[i],cex=0.5,font=2)
}
rx = 0.017*1.3*2.5; ry = 0.12*1.3; theta = seq(-pi,pi,length=100)
for(i in 1:nrow(d)){
#	size = log10(d1[i,"n"])
#	rx2 = 0.017*size; ry2 = 0.12*size; theta = seq(-pi,pi,length=100)
#	polygon(2.05+(rx2*cos(theta)),60-i+ry2*sin(theta),col="violet",lwd=0.2)
#	text(-2,60-i,rownames(d1)[i],adj=1,font=2,cex=0.4)
	segments(log10(d1[i,"lower_95"]),60-i,log10(d1[i,"upper_95"]),60-i,lwd=0.3)
	polygon(log10(d1[i,"exp_coef"])+(rx*cos(theta)),60-i+ry*sin(theta),col=rgb(1,0,0,alpha=1),lwd=0.2)
}

plot(NULL,NULL,axes=F,xlim=c(-2,2),ylim=c(0,60),xlab="",ylab="")
segments(0,0,0,60,col="gray",lwd=1)
segments(-1,0,2,0,col="black",lwd=1)
lat = c((0.1)*(1:9),1:9,10*(1:9))
lat1 = c(0.1,1,10,100)
for(i in 1:length(lat)){
	segments(log10(lat[i]),0,log10(lat[i]),-0.2,lwd=1)
	segments(log10(lat1[i]),0,log10(lat1[i]),-0.4,lwd=1)
	text(log10(lat1[i]),-0.8,lat1[i],cex=0.5,font=2)
}
rx = 0.017*1.3*2.5; ry = 0.12*1.3; theta = seq(-pi,pi,length=100)
for(i in 1:nrow(d)){
#	size = log10(d2[i,"n"])
#	rx2 = 0.017*size; ry2 = 0.12*size; theta = seq(-pi,pi,length=100)
#	polygon(2.05+(rx2*cos(theta)),60-i+ry2*sin(theta),col="violet",lwd=0.2)
#	text(-2,60-i,rownames(d2)[i],adj=1,font=2,cex=0.4)
	segments(log10(d2[i,"lower_95"]),60-i,log10(d2[i,"upper_95"]),60-i,lwd=0.3)
	polygon(log10(d2[i,"exp_coef"])+(rx*cos(theta)),60-i+ry*sin(theta),col=rgb(0,0,1,alpha=1),lwd=0.2)
}

dev.off()
}






writeForest5 = function(d,d1,d2,path.pdf){

pdf(path.pdf,width=5,height=6)
par(mfrow=c(1,3))

plot(NULL,NULL,axes=F,xlim=c(-2,2),ylim=c(0,60),xlab="",ylab="")
segments(0,0,0,60,col="gray",lwd=1)
segments(-2,0,2,0,col="black",lwd=1)
lat = c((0.01)*(1:9),(0.1)*(1:9),1:9,10*(1:9))
lat1 = c(0.01,0.1,1,10,100)
for(i in 1:length(lat)){
	segments(log10(lat[i]),0,log10(lat[i]),-0.2,lwd=1)
	segments(log10(lat1[i]),0,log10(lat1[i]),-0.4,lwd=1)
	text(log10(lat1[i]),-0.8,lat1[i],cex=0.5,font=2)
}
rx = 0.017*1.3*2.5; ry = 0.12*1.3; theta = seq(-pi,pi,length=100)
for(i in 1:nrow(d)){
#	size = log10(d[i,"n"])
#	rx2 = 0.017*size*5; ry2 = 0.12*size; theta = seq(-pi,pi,length=100)
#	polygon(2.05+(rx2*cos(theta)),60-i+ry2*sin(theta),col="violet",lwd=0.2)
	text(-1,60-i,d[i,"n"],adj=1,font=2,cex=0.4)
	text(-2,60-i,rownames(d)[i],adj=1,font=2,cex=0.4)
	segments(log10(d[i,"lower_95"]),60-i,log10(d[i,"upper_95"]),60-i,lwd=0.3)
	polygon(log10(d[i,"exp_coef"])+(rx*cos(theta)),60-i+ry*sin(theta),col=rgb(1,0,1,alpha=1),lwd=0.2)
}

plot(NULL,NULL,axes=F,xlim=c(-2,2),ylim=c(0,60),xlab="",ylab="")
segments(0,0,0,60,col="gray",lwd=1)
segments(-2,0,2,0,col="black",lwd=1)
lat = c((0.01)*(1:9),(0.1)*(1:9),1:9,10*(1:9))
lat1 = c(0.01,0.1,1,10,100)
for(i in 1:length(lat)){
	segments(log10(lat[i]),0,log10(lat[i]),-0.2,lwd=1)
	segments(log10(lat1[i]),0,log10(lat1[i]),-0.4,lwd=1)
	text(log10(lat1[i]),-0.8,lat1[i],cex=0.5,font=2)
}
rx = 0.017*1.3*2.5; ry = 0.12*1.3; theta = seq(-pi,pi,length=100)
for(i in 1:nrow(d)){
#	size = log10(d1[i,"n"])
#	rx2 = 0.017*size; ry2 = 0.12*size; theta = seq(-pi,pi,length=100)
#	polygon(2.05+(rx2*cos(theta)),60-i+ry2*sin(theta),col="violet",lwd=0.2)
#	text(-2,60-i,rownames(d1)[i],adj=1,font=2,cex=0.4)
	segments(log10(d1[i,"lower_95"]),60-i,log10(d1[i,"upper_95"]),60-i,lwd=0.3)
	polygon(log10(d1[i,"exp_coef"])+(rx*cos(theta)),60-i+ry*sin(theta),col=rgb(1,0,0,alpha=1),lwd=0.2)
}

plot(NULL,NULL,axes=F,xlim=c(-2,2),ylim=c(0,60),xlab="",ylab="")
segments(0,0,0,60,col="gray",lwd=1)
segments(-2,0,2,0,col="black",lwd=1)
lat = c((0.01)*(1:9),(0.1)*(1:9),1:9,10*(1:9))
lat1 = c(0.01,0.1,1,10,100)
for(i in 1:length(lat)){
	segments(log10(lat[i]),0,log10(lat[i]),-0.2,lwd=1)
	segments(log10(lat1[i]),0,log10(lat1[i]),-0.4,lwd=1)
	text(log10(lat1[i]),-0.8,lat1[i],cex=0.5,font=2)
}
rx = 0.017*1.3*2.5; ry = 0.12*1.3; theta = seq(-pi,pi,length=100)
for(i in 1:nrow(d)){
#	size = log10(d2[i,"n"])
#	rx2 = 0.017*size; ry2 = 0.12*size; theta = seq(-pi,pi,length=100)
#	polygon(2.05+(rx2*cos(theta)),60-i+ry2*sin(theta),col="violet",lwd=0.2)
#	text(-2,60-i,rownames(d2)[i],adj=1,font=2,cex=0.4)
	segments(log10(d2[i,"lower_95"]),60-i,log10(d2[i,"upper_95"]),60-i,lwd=0.3)
	polygon(log10(d2[i,"exp_coef"])+(rx*cos(theta)),60-i+ry*sin(theta),col=rgb(0,0,1,alpha=1),lwd=0.2)
}

dev.off()
}





