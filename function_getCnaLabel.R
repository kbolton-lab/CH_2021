
getCnaLabel = function(s){
	splited = unlist(strsplit(s,":"))
	type=NA
	prf = ifelse(splited[2]=="CNN-LOH","",ifelse(splited[2]=="Duplication","+",ifelse(splited[2]=="Deletion","del(","Unknown")))
	suf = ifelse(splited[2]=="CNN-LOH","UPD",ifelse(splited[2]=="Duplication","",ifelse(splited[2]=="Deletion",")","Unknown")))
	return(paste0(prf,splited[1],splited[3],suf))
}

