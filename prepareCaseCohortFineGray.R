

##############################
### case cohort, fine-gray ###
##############################

ccfg = function(cmp_data_org){

subcohort = cmp_data_org[which(cmp_data_org$is_subcohort==1),]
ccoh_ind = rep(NA,nrow(subcohort))
for(i in 1:nrow(subcohort)){
	if(subcohort$fail[i]=="Event"){
		ccoh_ind[i] = 1
	}else if(subcohort$fail[i]=="death"){
		ccoh_ind[i] = 3
	}else{
		ccoh_ind[i] = 0
	}
}
subcohort = cbind(subcohort,ccoh_ind)

case = cmp_data_org[which(cmp_data_org$is_subcohort==0 & cmp_data_org$fail=="Event"),]
ccoh_ind = rep(2,nrow(case))
case = cbind(case,ccoh_ind)

case_sub = rbind(case,subcohort)
an_entry = an_exit = an_ind = w = id = c()
analytic = as.data.frame(matrix(ncol=ncol(case_sub)))
colnames(analytic) = colnames(case_sub)
for(i in 1:nrow(case_sub)){
	if(case_sub$ccoh_ind[i]==1 | case_sub$ccoh_ind[i]==2){
		if(case_sub$ccoh_ind[i]==1){
			an_entry = c(an_entry,0)
			an_exit = c(an_exit,case_sub$ftime[i] - .0001)
			an_ind = c(an_ind,0)
			w = c(w,N_cohort/n_subcohort)
			id = c(id,rownames(case_sub)[i])
			analytic = rbind(analytic,case_sub[i,])
		}
		
		an_entry = c(an_entry,case_sub$ftime[i] - .0001)
		an_exit = c(an_exit,case_sub$ftime[i])
		an_ind = c(an_ind,1)
		w = c(w,1)
		id = c(id,rownames(case_sub)[i])
		analytic = rbind(analytic,case_sub[i,])
		
	}else if(case_sub$ccoh_ind[i]==0){
		an_entry = c(an_entry,0)
		an_exit = c(an_exit,case_sub$ftime[i])
		an_ind = c(an_ind,0)
		w = c(w,N_cohort/n_subcohort)
		id = c(id,rownames(case_sub)[i])
		analytic = rbind(analytic,case_sub[i,])
		
	}else if(case_sub$ccoh_ind[i]==3){
		an_entry = c(an_entry,0)
		an_exit = c(an_exit,case_sub$ftime[i])
		an_ind = c(an_ind,2)
		w = c(w,N_cohort/n_subcohort)
		id = c(id,rownames(case_sub)[i])
		analytic = rbind(analytic,case_sub[i,])
	}
}

an_ind = factor(an_ind,levels = c(0,1,2),labels = c("Censored","Event","Others"))

analytic2 = na.omit(analytic)
analytic2 = cbind(analytic2,an_entry,an_exit,an_ind,w,id)
analytic2$an_exit[which(analytic2$an_exit == 0)] = .0001
pdata <- finegray(Surv(an_entry,an_exit,an_ind)~.,weights=w,id=id,data=analytic2,etype="Event")

return(pdata)
}




#################
### fine-gray ###
#################

fg = function(cmp_data_org){
pdata <- finegray(Surv(ftime,fail)~.,data=cmp_data_org,etype="Event")
return(pdata)
}

