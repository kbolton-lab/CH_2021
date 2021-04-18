#################
### root path ###
#################

path = "/Users/saiki/Desktop/sync4/"


###########################
### import genomic data ###
###########################

source(paste0(path,"analysis/ImportGenomicData_HM_development.R"))
cmp_data_org = cmp_data_org[,setdiff(colnames(cmp_data_org),c("all_time_to_traf","traf_fail","status_HM","status_MYE","status_LYM","status_MDS","status_AML","status_NHL","status_MM","status_CVD","status_B","status_T","status_MPN","status_ALL","status_CML","status_CLL","all_disease"))]
cmp_data_org = na.omit(cmp_data_org)
path = "/Users/saiki/Desktop/sync4/"


######################
### select disease ###
######################

path_2 = paste0(path,"/analysis/fine_gray/HM/")


#################
### functions ###
#################

source(paste0(path,"/analysis/show_di.R"))
source(paste0(path,"/analysis/writeForest.R"))
source(paste0(path,"/analysis/getCnaLabel.R"))
check = function(x){length(unique(x$id))}
makeSf = function(d){ m = coxph(Surv(fgstart, fgstop, fgstatus) ~ 1, weight=fgwt, data=d); sf = survfit(m); return(sf) }


##############################
### case cohort, fine-gray ###
##############################

source("prepareCaseCohortFineGray.R")
pdata = ccfg(cmp_data_org)
pdata$fgstart = pdata$fgstart/30
pdata$fgstop = pdata$fgstop/30




##############################
### calculate hazard ratio ###
##############################

Names = c(
c("is_mut","is_mut_small","is_mut_large","is_mut_1","is_mut_2","is_mut_3","is_mut_n1","is_mut_n2","is_mut_n3"),
c("is_cna","is_cna_small","is_cna_large","is_cna_1","is_cna_2","is_cna_3","is_cna_n1","is_cna_n2","is_cna_n3"),
c("is_any","is_either","is_both","is_mut_only","is_cna_only"),
c("is_cis","is_dis"),
c("all_cis_tet2","all_cis_tp53","all_cis_dnmt3a","all_cis_jak2","all_cis_ezh2","all_cis_cbl","all_cis_runx1"),
c("x2pLOH","x4qLOH","x17pLOH"),c("is_tp53_ddpcr","is_jak2_ddpcr")
)
Labels = c(
c("SNV","SNV small","SNV large","VAF <5%","VAF 5-10%","VAF >10%","is_mut_n1","is_mut_n2","is_mut_n3"),
c("is_cna","is_cna_small","is_cna_large","is_cna_1","is_cna_2","is_cna_3","is_cna_n1","is_cna_n2","is_cna_n3"),
c("is_any","is_either","is_both","is_mut_only","is_cna_only"),
c("is_cis","is_dis"),
c("all_cis_tet2","all_cis_tp53","all_cis_dnmt3a","all_cis_jak2","all_cis_ezh2","all_cis_cbl","all_cis_runx1"),
c("x2pLOH","x4qLOH","x17pLOH"),c("is_tp53_ddpcr","is_jak2_ddpcr")
)

res_all = matrix(NA,ncol=11,nrow=length(Names))
rownames(res_all) = Names
for(i in 1:length(Names)){
cur_pdata = pdata[which(pdata[,Names[i]]==1 | pdata$is_none==1),] 
colnames(cur_pdata)[which(colnames(cur_pdata)==Names[i])] = "cur_alt"
n = length(unique(cur_pdata[which(cur_pdata$cur_alt==1),"id"]))
if(n<=5) next
m = coxph(Surv(fgstart,fgstop,fgstatus) ~ cur_alt+fac_age+is_male+is_OEE_1_0+is_OEE_1_2,weight=fgwt,data=cur_pdata)
res_all[i,] = c(summary(m)$coef[1,],summary(m)$conf.int[1,],n)
}

rownames(res_all) = Labels
colnames(res_all) = c(colnames(summary(m)$coef),colnames(summary(m)$conf.int),"n")
res_all


## cumulative incidence ##
## subjects without CH
pdata_none = pdata[which(pdata$is_none==1),]
m_none = coxph(Surv(fgstart, fgstop, fgstatus) ~ 1, weight=fgwt, data=pdata_none)
sf_none = survfit(m_none)


## subjects with SNV
pdata_mut_none = pdata[which(pdata$is_mut==0),]
pdata_mut_large = pdata[which(pdata$is_mut_large==1),]
pdata_mut_small= pdata[which(pdata$is_mut_small==1),]
sf_mut_none = survfit(coxph(Surv(fgstart, fgstop, fgstatus) ~ 1, weight=fgwt, data=pdata_mut_none))
sf_mut_large = survfit(coxph(Surv(fgstart, fgstop, fgstatus) ~ 1, weight=fgwt, data=pdata_mut_large))
sf_mut_small = survfit(coxph(Surv(fgstart, fgstop, fgstatus) ~ 1, weight=fgwt, data=pdata_mut_small))

pdf(paste0(path_2,"cum_inc_snv_small_large_fg.pdf"),width=4,height=5)
par(mfrow=c(1,1),mai=rep(0.2,4))
show_di(list(sf_mut_small,sf_mut_large,sf_mut_none),"SNV",c("pink","red","gray"),ymax_df=0.006)
dev.off()

pdata_mut_large_none = rbind(pdata_mut_large,pdata_none)
m_mut_large_none <- coxph(Surv(fgstart, fgstop, fgstatus) ~ is_mut + fac_age + is_male + is_OEE_1_0 +  is_OEE_1_2, weight=fgwt, data=pdata_mut_large_none)
pdata_mut_small_none = rbind(pdata_mut_small,pdata_none)
m_mut_small_none <- coxph(Surv(fgstart, fgstop, fgstatus) ~ is_mut + fac_age + is_male + is_OEE_1_0 +  is_OEE_1_2, weight=fgwt, data=pdata_mut_small_none)
comb_all = rbind(
c(summary(m_mut_large_none)$coef[1,],summary(m_mut_large_none)$conf.int[1,]),
c(summary(m_mut_small_none)$coef[1,],summary(m_mut_small_none)$conf.int[1,])
)
comb_all


## subjects with CNA
pdata_cna_none = pdata[which(pdata$is_cna==0),]
pdata_cna_large = pdata[which(pdata$is_cna_large==1),]
pdata_cna_small= pdata[which(pdata$is_cna_small==1),]
sf_cna_none = survfit(coxph(Surv(fgstart, fgstop, fgstatus) ~ 1, weight=fgwt, data=pdata_cna_none))
sf_cna_large = survfit(coxph(Surv(fgstart, fgstop, fgstatus) ~ 1, weight=fgwt, data=pdata_cna_large))
sf_cna_small = survfit(coxph(Surv(fgstart, fgstop, fgstatus) ~ 1, weight=fgwt, data=pdata_cna_small))

pdf(paste0(path_2,"cum_inc_cna_small_large_fg.pdf"),width=4,height=5)
par(mfrow=c(1,1),mai=rep(0.2,4))
show_di(list(sf_cna_small,sf_cna_large,sf_cna_none),"SNV",c("pink","red","gray"),ymax_df=0.003)
dev.off()

pdata_cna_large_none = rbind(pdata_cna_large,pdata_none)
m_cna_large_none = coxph(Surv(fgstart, fgstop, fgstatus) ~ is_cna + fac_age + is_male + is_OEE_1_0 +  is_OEE_1_2, weight=fgwt, data=pdata_cna_large_none)
pdata_cna_small_none = rbind(pdata_cna_small,pdata_none)
m_cna_small_none = coxph(Surv(fgstart, fgstop, fgstatus) ~ is_cna + fac_age + is_male + is_OEE_1_0 +  is_OEE_1_2, weight=fgwt, data=pdata_cna_small_none)
comb_all = rbind(
c(summary(m_cna_large_none)$coef[1,],summary(m_cna_large_none)$conf.int[1,]),
c(summary(m_cna_small_none)$coef[1,],summary(m_cna_small_none)$conf.int[1,])
)
comb_all


# subjects with combined SNV + CNA
pdata_both = pdata[which(pdata$is_mut==1 & pdata$is_cna==1),]
m_both = coxph(Surv(fgstart, fgstop, fgstatus) ~ 1, weight=fgwt, data=pdata_both)
sf_both = survfit(m_both)
pdata_either = pdata[which(xor(pdata$is_mut==1 , pdata$is_cna==1)),]
m_either = coxph(Surv(fgstart, fgstop, fgstatus) ~ 1, weight=fgwt, data=pdata_either)
sf_either = survfit(m_either)
pdata_mut_alone = pdata[which(pdata$is_mut==1 & pdata$is_cna==0),]
m_mut_alone = coxph(Surv(fgstart, fgstop, fgstatus) ~ 1, weight=fgwt, data=pdata_mut_alone)
sf_mut_alone = survfit(m_mut_alone)
pdata_cna_alone = pdata[which(pdata$is_mut==0 & pdata$is_cna==1),]
m_cna_alone = coxph(Surv(fgstart, fgstop, fgstatus) ~ 1, weight=fgwt, data=pdata_cna_alone)
sf_cna_alone = survfit(m_cna_alone)

pdf(paste0(path_2,"both_cum_inc_fg.pdf"),width=4,height=5)
par(mfrow=c(1,1),mai =rep(0.2,4))
show_di(list(sf_both,sf_either,sf_mut_alone,sf_cna_alone,sf_none),"Fine-Gray",c("purple","violet","red","skyblue","gray"),ymax_df=0.003)
dev.off()

pdata_both_mut = rbind(pdata_both,pdata_mut_alone)
m_both_mut <- coxph(Surv(fgstart, fgstop, fgstatus) ~ is_both + fac_age + is_male + is_OEE_1_0 +  is_OEE_1_2, weight=fgwt, data=pdata_both_mut)
pdata_both_cna = rbind(pdata_both,pdata_cna_alone)
m_both_cna <- coxph(Surv(fgstart, fgstop, fgstatus) ~ is_both + fac_age + is_male + is_OEE_1_0 +  is_OEE_1_2, weight=fgwt, data=pdata_both_cna)
pdata_both_either = rbind(pdata_both,pdata_either)
m_both_either <- coxph(Surv(fgstart, fgstop, fgstatus) ~ is_both + fac_age + is_male + is_OEE_1_0 +  is_OEE_1_2, weight=fgwt, data=pdata_both_either)
comb_all = rbind(
c(summary(m_both_mut)$coef[1,],summary(m_both_mut)$conf.int[1,]),
c(summary(m_both_cna)$coef[1,],summary(m_both_cna)$conf.int[1,]),
c(summary(m_both_either)$coef[1,],summary(m_both_either)$conf.int[1,])
)
comb_all


## individual SNV ##
pdf(paste0(path_2,"snv_cum_inc_fg_2.pdf"))
par(mfrow=c(4,4),mai =rep(0.2,4))
res_mut = matrix(ncol=11,nrow=ncol(all_mat_mut_is))
rownames(res_mut) = colnames(all_mat_mut_is)
colnames(res_mut) = c("coef","exp_coef","se_coef","z","p","exp_pcoef","exp_ncoef","lower_95","upper_95","n")

for(i in 1:ncol(all_mat_mut_is)){
cur_alt = colnames(all_mat_mut_is)[i]
pdata_cur = pdata[which(pdata[,cur_alt]==1 | pdata$is_none==1),c(cur_alt,"fgstart", "fgstop","fgstatus","fgwt","fac_age","is_male","is_OEE_1_0","is_OEE_1_2")]
n = length(unique(pdata[which(pdata[,cur_alt] == 1),"id"]))
if(n<5) next
#if(sum(pdata_cur[,cur_alt])<5) next
colnames(pdata_cur) = c("cur_alt","fgstart", "fgstop","fgstatus","fgwt","fac_age","is_male","is_OEE_1_0","is_OEE_1_2")
m_cur = coxph(Surv(fgstart,fgstop,fgstatus) ~ cur_alt+fac_age+is_male+is_OEE_1_0+is_OEE_1_2,weight=fgwt,data=pdata_cur)
res_mut[i,] = c(summary(m_cur)$coef["cur_alt",],summary(m_cur)$conf.int["cur_alt",],n)
pdata_sf = pdata[which(pdata[,cur_alt]==1),]
m_sf = coxph(Surv(fgstart, fgstop, fgstatus) ~ 1, weight=fgwt, data=pdata_sf)
sf_cur = survfit(m_sf)
}
dev.off()
res_mut


## individual CNA ##
cna_list = c()
pdf(paste0(path_2,"cna_cum_inc_fg_2.pdf"))
par(mfrow=c(4,4),mai =rep(0.2,4))
res_cna = matrix(ncol=10,nrow=ncol(all_mat_cna_is))
rownames(res_cna) = colnames(all_mat_cna_is)
colnames(res_cna) = c("coef","exp_coef","se_coef","z","p","exp_pcoef","exp_ncoef","lower_95","upper_95","n")
for(i in 1:ncol(all_mat_cna_is)){
cur_alt = gsub("^","X",gsub("-",".",gsub(":",".",colnames(all_mat_cna_is)[i])))
pdata_cur = pdata[which(pdata[,cur_alt]==1 | pdata$is_none==1),c(cur_alt,"fgstart", "fgstop","fgstatus","fgwt","fac_age","is_male","is_OEE_1_0","is_OEE_1_2")]
n = length(unique(pdata[which(pdata[,cur_alt] == 1),"id"]))
if(n<5) next
#if(sum(pdata_cur[,cur_alt])<5) next
colnames(pdata_cur) = c("cur_alt","fgstart", "fgstop","fgstatus","fgwt","fac_age","is_male","is_OEE_1_0","is_OEE_1_2")
m_cur = coxph(Surv(fgstart,fgstop,fgstatus) ~ cur_alt+fac_age+is_male+is_OEE_1_0+is_OEE_1_2,weight=fgwt,data=pdata_cur)
}
dev.off()
res_cna


## End of analysis ##
q()


