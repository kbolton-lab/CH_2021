
library(RColorBrewer)

show_di = function(d_lis,label,col_lines_i="red",ymax_df=NULL,lat1s=NULL,lat2s=NULL){
	col_polygon_i = adjustcolor(col_lines_i,alpha=0.3)
	
	xmax = ymax = 0
	for(i in 1:length(d_lis)){
	xmax = max(xmax,max(d_lis[[i]]$time))
	up = d_lis[[i]]$cumhaz + 1.96*d_lis[[i]]$std.err
	ymax = max(ymax,min(1,ceiling(max(up*20))/20))
	}
	
	if(!is.null(ymax_df)){
		ymax = ymax_df
	}
	
	plot(
	NULL,NULL,axes=F,
	xlim=c(-0.1*1,1),ylim=c(-0.5*ymax,ymax),
	xlab="",ylab=""
	)
	
	segments(0,0,0,ymax)
	segments(0,ymax,1,ymax)
	segments(0,0,1,0)
	segments(1,0,1,ymax)
#	text(0.5,-0.3*ymax,label_i,font=4)
	text(0.5,-0.3*ymax,label,font=2)
	
	for(i in 1:length(d_lis)){
		di = d_lis[[i]]
		
		up = di$cumhaz + 1.96*di$std.err
		lw = di$cumhaz - 1.96*di$std.err
		up[which(up>ymax)] = ymax
		lw[which(lw<0)] = 0
		
		diy = c(lw,rev(up))
		dix = c(di$time,rev(di$time))
		polygon(dix/xmax,diy,col=col_polygon_i[i],border="transparent")
		lines(sort(c(0,min(di$time/xmax),c(di$time/xmax,di$time/xmax)[-1])),sort(c(0,0,c(di$cumhaz,di$cumhaz)))[-length(di$cumhaz)*2],col=col_lines_i[i])
	}
	
	lat1 = setdiff((0:15)*10,(0:7)*20)
	lat2 = (0:7)*20
	segments(lat1/xmax,0,lat1/xmax,-ymax/200)
	segments(lat2/xmax,0,lat2/xmax,-ymax/100)
	text(lat2/xmax,-ymax/20,lat2,adj=0.5,cex=0.6)
	
	if(is.null(lat1s)){
		lat1 = setdiff((0:10)*ymax/10,(0:5)*2*ymax/10)
		lat2 = (0:5)*2*ymax/10
		segments(0,lat1,-0.005,lat1)
		segments(0,lat2,-0.01,lat2)
		text(-0.025,lat2,lat2,adj=1,cex=0.6)
	}else{
		lat1 = lat1s
		lat2 = lat2s
		segments(0,lat1,-0.005,lat1)
		segments(0,lat2,-0.01,lat2)
		text(-0.025,lat2,lat2,adj=1,cex=0.6)
	}

}




show_di_2 = function(d_lis,label,col_lines_i="red",ymax_df=NULL,lat1s=NULL,lat2s=NULL){
	col_polygon_i = adjustcolor(col_lines_i,alpha=0.3)
	
	xmax = ymax = 0
	for(i in 1:length(d_lis)){
	xmax = max(xmax,max(d_lis[[i]]$time))
	up = d_lis[[i]]$cumhaz + 1.96*d_lis[[i]]$std.err
	ymax = max(ymax,min(1,ceiling(max(up*20))/20))
	}
	
	if(!is.null(ymax_df)){
		ymax = ymax_df
	}
	ymax=1
	
	plot(
	NULL,NULL,axes=F,
	xlim=c(-0.1*1,1),ylim=c(-0.5*ymax,ymax),
	xlab="",ylab=""
	)
	
	segments(0,0,0,ymax)
	segments(0,ymax,1,ymax)
	segments(0,0,1,0)
	segments(1,0,1,ymax)
#	text(0.5,-0.3*ymax,label_i,font=4)
	text(0.5,-0.3*ymax,label,font=2)
	
	for(i in 1:length(d_lis)){
		di = d_lis[[i]]
		
		up = di$upper
		lw = di$lower
		up[which(up>ymax)] = ymax
		lw[which(lw<0)] = 0
		
		diy = c(lw,rev(up))
		dix = c(di$time,rev(di$time))
		polygon(dix/xmax,diy,col=col_polygon_i[i],border="transparent")
		lines(sort(c(0,min(di$time/xmax),c(di$time/xmax,di$time/xmax)[-1])),sort(c(di$surv[1],di$surv[1],c(di$surv,di$surv)),decreasing=T)[-(length(di$surv)+1)*2],col=col_lines_i[i])
	}
	
	lat1 = setdiff((0:15)*10,(0:7)*20)
	lat2 = (0:7)*20
	segments(lat1/xmax,0,lat1/xmax,-ymax/200)
	segments(lat2/xmax,0,lat2/xmax,-ymax/100)
	text(lat2/xmax,-ymax/20,lat2,adj=0.5,cex=0.6)
	
	if(is.null(lat1s)){
		lat1 = setdiff((0:10)*ymax/10,(0:5)*2*ymax/10)
		lat2 = (0:5)*2*ymax/10
		segments(0,lat1,-0.005,lat1)
		segments(0,lat2,-0.01,lat2)
		text(-0.025,lat2,lat2,adj=1,cex=0.6)
	}else{
		lat1 = lat1s
		lat2 = lat2s
		segments(0,lat1,-0.005,lat1)
		segments(0,lat2,-0.01,lat2)
		text(-0.025,lat2,lat2,adj=1,cex=0.6)
	}

}

