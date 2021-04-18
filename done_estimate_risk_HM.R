
## import genomic data ##
source(paste0(path,"ImportGenomicData_HM.R"))
cmp_data_org = cmp_data_org[,setdiff(colnames(cmp_data_org),c("all_time_to_traf","traf_fail","all_disease"))]
cmp_data_org = na.omit(cmp_data_org)


## select disease
cmp_data_org$fail = cmp_data_org$status_HM


## import blood cel counts ##
## data format: WBC Hb Ht PLT
cbc_dat = ~~~ 
cmp_data_org_2 = cbind(cmp_data_org,cbc_dat)
cmp_data_org_2 = na.omit(cmp_data_org_2)


## import functions ##
source("show_di.R")
source("writeForest.R")
source("getCnaLabel.R")
makeSf = function(d){ m = coxph(Surv(fgstart, fgstop, fgstatus) ~ 1, weight=fgwt, data=d); sf = survfit(m); return(sf) }
sf10 = function(sf){
print(1-sf$surv[max(which(sf$time<120))])
print((sf$cumhaz - 1.96*sf$std.err)[max(which(sf$time<120))])
print((sf$cumhaz + 1.96*sf$std.err)[max(which(sf$time<120))])
}

## case cohort, fine-gray ##
source("prepareCaseCohortFineGray.R")
pdata = ccfg(cmp_data_org)
pdata_2 = ccfg(cmp_data_org_2)


## cbc nomral and abnormal ##
Hb_l_thld = 11; Hb_h_m_thld = 16.5; Hb_h_f_thld = 16
Ht_h_thld = 49
WBC_l_thld = 3000; WBC_h_thld = 10000
Plt_l_thld = 10; Plt_h_thld = 45

pdata_ab = pdata_2[which(
(pdata_2$Ht >= Ht_h_thld) | 
((pdata_2$Hb >= Hb_h_m_thld & pdata_2$is_male==1) | (pdata_2$Hb >= Hb_h_f_thld & pdata_2$is_male==0)) | 
(pdata_2$WBC >= WBC_h_thld) | 
(pdata_2$PLT >= Plt_h_thld) |
(pdata_2$Hb <= Hb_l_thld) | 
(pdata_2$WBC <= WBC_l_thld) | 
(pdata_2$PLT<= Plt_l_thld)
),]

pdata_ab_multi = pdata_2[which(
(pdata_2$Hb <= Hb_l_thld) +
(pdata_2$WBC <= WBC_l_thld) + 
(pdata_2$PLT<= Plt_l_thld) >=2
),]

pdata_nm = pdata_2[which(
(pdata_2$Ht < Ht_h_thld) & 
((pdata_2$Hb < Hb_h_m_thld & pdata_2$is_male==1) | (pdata_2$Hb < Hb_h_f_thld & pdata_2$is_male==0)) & 
(pdata_2$WBC < WBC_h_thld) &
(pdata_2$PLT < Plt_h_thld) &
(pdata_2$Hb > Hb_l_thld) & 
(pdata_2$WBC > WBC_l_thld) & 
(pdata_2$PLT > Plt_l_thld)
),]

cmp_data_org_2_ab = cmp_data_org_2[which(
(cmp_data_org_2$Ht >= Ht_h_thld) | 
((cmp_data_org_2$Hb >= Hb_h_m_thld & cmp_data_org_2$is_male==1) | (cmp_data_org_2$Hb >= Hb_h_f_thld & cmp_data_org_2$is_male==0)) | 
(cmp_data_org_2$WBC >= WBC_h_thld) | 
(cmp_data_org_2$PLT >= Plt_h_thld) |
(cmp_data_org_2$Hb <= Hb_l_thld) | 
(cmp_data_org_2$WBC <= WBC_l_thld) | 
(cmp_data_org_2$PLT<= Plt_l_thld)),]

cmp_data_org_2_nm = cmp_data_org_2[which(
(cmp_data_org_2$Ht < Ht_h_thld) & 
((cmp_data_org_2$Hb < Hb_h_m_thld & cmp_data_org_2$is_male==1) | (cmp_data_org_2$Hb < Hb_h_f_thld & cmp_data_org_2$is_male==0)) & 
(cmp_data_org_2$WBC < WBC_h_thld) &
(cmp_data_org_2$PLT < Plt_h_thld) &
(cmp_data_org_2$Hb > Hb_l_thld) & 
(cmp_data_org_2$WBC > WBC_l_thld) & 
(cmp_data_org_2$PLT > Plt_l_thld)),]


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

res_all = matrix(NA,ncol=10,nrow=length(Names))
rownames(res_all) = Names
for(i in 1:length(Names)){
cur_pdata = pdata[which(pdata[,Names[i]]==1 | pdata$is_none==1),] 
colnames(cur_pdata)[which(colnames(cur_pdata)==Names[i])] = "cur_alt"
n = length(unique(cur_pdata[which(cur_pdata$cur_alt==1),"id"]))
if(n<=5) next
m = coxph(Surv(fgstart,fgstop,fgstatus) ~ cur_alt+fac_age+is_male+is_OEE_1_0+is_OEE_1_2,weight=fgwt,data=cur_pdata)
# m = coxph(Surv(fgstart,fgstop,fgstatus) ~ cur_alt,weight=fgwt,data=cur_pdata)
res_all[i,] = c(summary(m)$coef[1,],summary(m)$conf.int[1,],n)
}

rownames(res_all) = Labels
colnames(res_all) = c(colnames(summary(m)$coef),colnames(summary(m)$conf.int),"n")

## output data for Fig.5b ##
write.table(res_all,file=paste0(path_2,"hazard_ratio_all.txt"),col.name=T,row.name=T,sep="\t",quote=F)



############################
### cumulative incidence ###
############################

## no alterations ##
pdata_none = pdata[which(pdata$is_none==1),]
m_none = coxph(Surv(fgstart, fgstop, fgstatus) ~ 1, weight=fgwt, data=pdata_none)
sf_none = survfit(m_none)
cum10_noCH = 1-sf_none$surv[max(which(sf_none$time<120))]

## All CH ##
pdata_any = pdata[which(pdata$is_any==1),]
m_any = coxph(Surv(fgstart, fgstop, fgstatus) ~ 1, weight=fgwt, data=pdata_any)
sf_any = survfit(m_any)
pdf(paste0(path_2,"cum_inc_ch_fg.pdf"),width=4,height=5)
par(mfrow=c(1,1),mai=rep(0.2,4))
show_di(list(sf_any,sf_none),"All CH",c("red","gray"),ymax_df=0.015)
dev.off()
cum10_wtCH = 1-sf_any$surv[max(which(sf_any$time<120))]

## attributable mortality within 10 years ##
att10_CH = cum10_wtCH - cum10_noCH
att10_CH

## SNV CH ##
pdata_mut = pdata[which(pdata$is_mut==1),]
m_mut = coxph(Surv(fgstart, fgstop, fgstatus) ~ 1, weight=fgwt, data=pdata_mut)
sf_mut = survfit(m_mut)

## 10 year mortality in SNV+ subjects ##
CumMot10.SNV = 1-sf_mut$surv[max(which(sf_mut$time<120))]
CumMot10.SNV

## CNA CH ##
pdata_cna = pdata[which(pdata$is_cna==1),]
m_cna = coxph(Surv(fgstart, fgstop, fgstatus) ~ 1, weight=fgwt, data=pdata_cna)
sf_cna = survfit(m_cna)

## 10 year mortality in CNA+ subjects ##
CumMot10.CNA = 1-sf_mut$surv[max(which(sf_mut$time<120))]
CumMot10.CNA

## draw Fig.4a ##
pdf(paste0(path_2,"cum_inc_cna_fg.pdf"),width=4,height=5)
par(mfrow=c(1,1),mai=rep(0.2,4))
show_di(list(sf_any,sf_mut,sf_cna,sf_none),"",c("purple","red","skyblue","gray"),ymax_df=0.03)
dev.off()



## SNV + CNA ##
pdata_both = pdata[which(pdata$is_mut==1 & pdata$is_cna==1),]
m_both = coxph(Surv(fgstart, fgstop, fgstatus) ~ 1, weight=fgwt, data=pdata_both)
sf_both = survfit(m_both)
pdata_mut_alone = pdata[which(pdata$is_mut==1 & pdata$is_cna==0),]
m_mut_alone = coxph(Surv(fgstart, fgstop, fgstatus) ~ 1, weight=fgwt, data=pdata_mut_alone)
sf_mut_alone = survfit(m_mut_alone)
pdata_cna_alone = pdata[which(pdata$is_mut==0 & pdata$is_cna==1),]
m_cna_alone = coxph(Surv(fgstart, fgstop, fgstatus) ~ 1, weight=fgwt, data=pdata_cna_alone)
sf_cna_alone = survfit(m_cna_alone)

## draw Fig.4f ##
pdf(paste0(path_2,"both_cum_inc_fg.pdf"),width=4,height=5)
par(mfrow=c(1,1),mai =rep(0.2,4))
show_di(list(sf_both,sf_mut_alone,sf_cna_alone,sf_none),"Fine-Gray",c("purple","red","skyblue","gray"),ymax_df=0.025)
dev.off()



## Comparison: Both vs. SNV alone or Both vs. CNA alone ##
pdata_both = pdata[which(pdata$is_cna==1 & pdata$is_cna==1),]
pdata_mut_alone = pdata[which(pdata$is_mut==1 & pdata$is_cna==0),]
pdata_cna_alone = pdata[which(pdata$is_mut==0 & pdata$is_cna==1),]

pdata_both_mut = rbind(pdata_both,pdata_mut_alone)
m_both_mut <- coxph(Surv(fgstart, fgstop, fgstatus) ~ is_both + fac_age + is_male + is_OEE_1_0 +  is_OEE_1_2, weight=fgwt, data=pdata_both_mut)
pdata_both_cna = rbind(pdata_both,pdata_cna_alone)
m_both_cna <- coxph(Surv(fgstart, fgstop, fgstatus) ~ is_both + fac_age + is_male + is_OEE_1_0 +  is_OEE_1_2, weight=fgwt, data=pdata_both_cna)

## HR and 95%conf.int.: Both vs. SNV alone ##
c(summary(m_both_mut)$coef[1,],summary(m_both_mut)$conf.int[1,])
# HR=2.84, 95%CI:2.14-3.78

## HR and 95%conf.int.: Both vs. CNA alone ##
c(summary(m_both_cna)$coef[1,],summary(m_both_cna)$conf.int[1,])
# HR=2.64, 95%CI:1.94-3.60


## Cumilative mortality and hazard ratios for individual SNV ##
## draw Fig.S8 ##
pdf(paste0(path_2,"snv_cum_inc_fg_2.pdf"))
par(mfrow=c(4,4),mai =rep(0.2,4))
res_mut = matrix(ncol=10,nrow=ncol(all_mat_mut_is))
rownames(res_mut) = colnames(all_mat_mut_is)
colnames(res_mut) = c("coef","exp_coef","se_coef","z","p","exp_pcoef","exp_ncoef","lower_95","upper_95","n")

for(i in 1:ncol(all_mat_mut_is)){
cur_alt = colnames(all_mat_mut_is)[i]
pdata_cur = pdata[which(pdata[,cur_alt]==1 | pdata$is_none==1),c(cur_alt,"fgstart", "fgstop","fgstatus","fgwt","fac_age","is_male","is_OEE_1_0","is_OEE_1_2")]
n = length(unique(pdata[which(pdata[,cur_alt] == 1),"id"]))
if(n<5) next
colnames(pdata_cur) = c("cur_alt","fgstart", "fgstop","fgstatus","fgwt","fac_age","is_male","is_OEE_1_0","is_OEE_1_2")
m_cur = coxph(Surv(fgstart,fgstop,fgstatus) ~ cur_alt+fac_age+is_male+is_OEE_1_0+is_OEE_1_2,weight=fgwt,data=pdata_cur)
res_mut[i,] = c(summary(m_cur)$coef["cur_alt",],summary(m_cur)$conf.int["cur_alt",],n)
pdata_sf = pdata[which(pdata[,cur_alt]==1),]
m_sf = coxph(Surv(fgstart, fgstop, fgstatus) ~ 1, weight=fgwt, data=pdata_sf)
sf_cur = survfit(m_sf)
show_di(list(sf_cur,sf_none),cur_alt,c("red","gray"))
}
dev.off()

## output data for Fig.5c ##
write.table(res_mut,file=paste0(path_2,"data_mut_forest.txt"),col.name=T,row.name=T,sep="\t",quote=F)


## Cumilative mortality and hazard ratios for individual CNA ##
## draw Fig.S9 ##
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
if(!cvd) m_cur = coxph(Surv(fgstart,fgstop,fgstatus) ~ cur_alt+fac_age+is_male+is_OEE_1_0+is_OEE_1_2,weight=fgwt,data=pdata_cur)
if(cvd) m_cur = coxph(Surv(fgstart,fgstop,fgstatus) ~ cur_alt+fac_age+is_male+is_OEE_1_0+is_OEE_1_2+all_BMI+all_HT+all_DM+all_HL+all_Smoke+all_Drink,weight=fgwt,data=pdata_cur)
res_cna[i,] = c(summary(m_cur)$coef["cur_alt",],summary(m_cur)$conf.int["cur_alt",],n)
pdata_sf = pdata[which(pdata[,cur_alt]==1),]
m_sf = coxph(Surv(fgstart, fgstop, fgstatus) ~ 1, weight=fgwt, data=pdata_sf)
sf_cur = survfit(m_sf)
if(sum(sf_cur$n.event)==0) next
if(is.element(i,c(22,29,57))) next
if(length(grep("Unknown",cur_alt))==1) next
#show_di(list(sf_cur,sf_none),getCnaLabel(colnames(all_mat_cna_is)[i]),c("red","gray"),ymax_df=0.25)
show_di(list(sf_cur,sf_none),getCnaLabel(colnames(all_mat_cna_is)[i]),c("red","gray"))
cna_list = c(cna_list,rownames(res_cna)[i])
}
dev.off()

## output data for Fig.5d ##
write.table(res_cna,file=paste0(path_2,"data_cna_foest.txt"),col.name=T,row.name=T,sep="\t",quote=F)


## calculate expected 10 year mortality by different number of alterations ##
s0 <- survfit(coxph(Surv(fgstart, fgstop, fgstatus) ~ 1, weight=fgwt, data=pdata_1[which(pdata_1$all_n_alt==0),]))
s1 <- survfit(coxph(Surv(fgstart, fgstop, fgstatus) ~ 1, weight=fgwt, data=pdata_1[which(pdata_1$all_n_alt==1),]))
s2 <- survfit(coxph(Surv(fgstart, fgstop, fgstatus) ~ 1, weight=fgwt, data=pdata_1[which(pdata_1$all_n_alt==2),]))
s3 <- survfit(coxph(Surv(fgstart, fgstop, fgstatus) ~ 1, weight=fgwt, data=pdata_1[which(pdata_1$all_n_alt==3),]))
s4 <- survfit(coxph(Surv(fgstart, fgstop, fgstatus) ~ 1, weight=fgwt, data=pdata_1[which(pdata_1$all_n_alt>=4),]))
c0 = 1-s0$surv[max(which(s0$time<120))]
c1 = 1-s1$surv[max(which(s1$time<120))]
c2 = 1-s2$surv[max(which(s2$time<120))]
c3 = 1-s3$surv[max(which(s3$time<120))]
c4 = 1-s4$surv[max(which(s4$time<120))]

## expected 10 year mortality considering both snv and cna ##
c_all = (
c0*sum(cmp_data_org$all_n_alt == 0) +
c1*sum(cmp_data_org$all_n_alt == 1) +
c2*sum(cmp_data_org$all_n_alt == 2) +
c3*sum(cmp_data_org$all_n_alt == 3) +
c4*sum(cmp_data_org$all_n_alt >= 4)
)/length(cmp_data_org$all_n_alt)
c_all # 0.008022594

## expected 10 year mortality considering snv alone ##
c_mut = (
c0*sum(cmp_data_org$all_n_mut == 0) +
c1*sum(cmp_data_org$all_n_mut == 1) +
c2*sum(cmp_data_org$all_n_mut == 2) +
c3*sum(cmp_data_org$all_n_mut == 3) +
c4*sum(cmp_data_org$all_n_mut >= 4)
)/length(cmp_data_org$all_n_mut)
c_mut # 0.006661236


## expected 10 year mortality considering cna alone ##
c_cna = (
c0*sum(cmp_data_org$all_n_cna == 0) +
c1*sum(cmp_data_org$all_n_cna == 1) +
c2*sum(cmp_data_org$all_n_cna == 2) +
c3*sum(cmp_data_org$all_n_cna == 3) +
c4*sum(cmp_data_org$all_n_cna >= 4)
)/length(cmp_data_org$all_n_cna)
c_cna# 0.006092834


## Combined analysis increased the estimation of 10 year mortality by ##
c_all - c_mut # 0.001361358
c_all - c_cna # 0.00192976


## analysis stratified by numbers of SNVs ##
## stratified Log rank test ##
cmp_data_g2 = cmp_data_org[which(cmp_data_org$all_n_mut>=2),]; cmp_data_g2$age_strata = as.numeric(cmp_data_g2$fac_age>7)
cmp_data_m1 = cmp_data_org[which(cmp_data_org$all_n_mut==1),]; cmp_data_m1$age_strata = as.numeric(cmp_data_m1$fac_age>7)
cmp_data_m2 = cmp_data_org[which(cmp_data_org$all_n_mut==2),]; cmp_data_m2$age_strata = as.numeric(cmp_data_m2$fac_age>7)
cmp_data_m3 = cmp_data_org[which(cmp_data_org$all_n_mut==3),]; cmp_data_m3$age_strata = as.numeric(cmp_data_m3$fac_age>7)
cmp_data_m4 = cmp_data_org[which(cmp_data_org$all_n_mut>=4),]; cmp_data_m4$age_strata = as.numeric(cmp_data_m4$fac_age>7)
survdiff(Surv(ftime,fail=="Event") ~ is_cna+strata(age_strata)+strata(is_male),data=cmp_data_m1)
survdiff(Surv(ftime,fail=="Event") ~ is_cna+strata(age_strata)+strata(is_male),data=cmp_data_m2)
survdiff(Surv(ftime,fail=="Event") ~ is_cna+strata(age_strata)+strata(is_male),data=cmp_data_m3)
survdiff(Surv(ftime,fail=="Event") ~ is_cna+strata(age_strata)+strata(is_male),data=cmp_data_g2)

## draw ED Fig.8 c-e ##
pdata_g2_nca = pdata[which(pdata$all_n_mut>=2 & pdata$is_cna==0),]; pdata_g2_nca$age_strata = as.numeric(pdata_g2_nca$fac_age>7)
pdata_m1_nca = pdata[which(pdata$all_n_mut==1 & pdata$is_cna==0),]; pdata_m1_nca$age_strata = as.numeric(pdata_m1_nca$fac_age>7)
pdata_m2_nca = pdata[which(pdata$all_n_mut==2 & pdata$is_cna==0),]; pdata_m2_nca$age_strata = as.numeric(pdata_m2_nca$fac_age>7)
pdata_m3_nca = pdata[which(pdata$all_n_mut==3 & pdata$is_cna==0),]; pdata_m3_nca$age_strata = as.numeric(pdata_m3_nca$fac_age>7)
pdata_m4_nca = pdata[which(pdata$all_n_mut>=4 & pdata$is_cna==0),]; pdata_m4_nca$age_strata = as.numeric(pdata_m4_nca$fac_age>7)
pdata_g2_cna = pdata[which(pdata$all_n_mut>=2 & pdata$is_cna==1),]; pdata_g2_nca$age_strata = as.numeric(pdata_g2_nca$fac_age>7)
pdata_m1_cna = pdata[which(pdata$all_n_mut==1 & pdata$is_cna==1),]; pdata_m1_nca$age_strata = as.numeric(pdata_m1_nca$fac_age>7)
pdata_m2_cna = pdata[which(pdata$all_n_mut==2 & pdata$is_cna==1),]; pdata_m2_nca$age_strata = as.numeric(pdata_m2_nca$fac_age>7)
pdata_m3_cna = pdata[which(pdata$all_n_mut==3 & pdata$is_cna==1),]; pdata_m3_nca$age_strata = as.numeric(pdata_m3_nca$fac_age>7)
pdata_m4_cna = pdata[which(pdata$all_n_mut>=4 & pdata$is_cna==1),]; pdata_m4_nca$age_strata = as.numeric(pdata_m4_nca$fac_age>7)
pdf("effect_of_cna_stratified_by_SNV.pdf",height=5.5,width=16)
par(mfrow=c(1,4))
show_di(list(makeSf(pdata_m1_nca),makeSf(pdata_m1_cna)),"1 SNV",c("red","purple"),ymax_df=0.06)
show_di(list(makeSf(pdata_m2_nca),makeSf(pdata_m2_cna)),"2 SNV",c("red","purple"),ymax_df=0.06)
show_di(list(makeSf(pdata_m3_nca),makeSf(pdata_m3_cna)),"3 SNV",c("red","purple"),ymax_df=0.1)
show_di(list(makeSf(pdata_g2_nca),makeSf(pdata_g2_cna)),">=2 SNV",c("red","purple"),ymax_df=0.06)
dev.off()






################################
### bi-allele vs mono-allele ###
################################

pdata_biallele = pdata[which(pdata$TET2==1 | pdata$DNMT3A==1 | pdata$JAK2==1 | pdata$EZH2==1 | pdata$CBL==1 | pdata$TP53==1),]

pdata_biallele$is_biallelic = as.numeric(
(pdata_biallele$is_jak2_ddpcr == 1) | (pdata_biallele$is_tp53_ddpcr == 1) | 
(pdata_biallele$all_cis_dnmt3a==1) | (pdata_biallele$all_cis_tet2 == 1) |
(pdata_biallele$all_cis_ezh2==1) | (pdata_biallele$all_cis_cbl==1))

m <- coxph(Surv(fgstart, fgstop, fgstatus) ~ is_biallelic + fac_age + is_male + is_OEE_1_0 +  is_OEE_1_2, weight=fgwt, data=pdata_biallele)
out = summary(m)

pdata_mono = pdata_biallele[which(pdata_biallele$is_biallelic==0),]
pdata_bi = pdata_biallele[which(pdata_biallele$is_biallelic==1),]

m_mono <- coxph(Surv(fgstart, fgstop, fgstatus) ~ 1, weight=fgwt, data=pdata_mono)
sf_mono = survfit(m_mono)

m_bi <- coxph(Surv(fgstart, fgstop, fgstatus) ~ 1, weight=fgwt, data=pdata_bi)
sf_bi = survfit(m_bi)

pdf(paste0(path_2,"bi_mono_fg.pdf"),width=7,height=10)
show_di(list(sf_bi,sf_mono,sf_none),">=2 alterations",c("purple","violet","gray"),ymax_df=0.1)
dev.off()


## Comparison: probable bi-allelic vs probable mono-allelic ##
pdata_none = pdata[which(pdata$is_none==1),]
m_none = coxph(Surv(fgstart, fgstop, fgstatus) ~ 1, weight=fgwt, data=pdata_none)
sf_none = survfit(m_none)

pdata_b = pdata
pdata_b$is_biallelic = as.numeric( (pdata_b$all_cis_jak2==1) | 
(pdata_b$all_cis_tp53==1) | (pdata_b$all_cis_dnmt3a==1) | (pdata_b$all_cis_tet2 == 1) | 
(pdata_b$all_cis_ezh2==1) | (pdata_b$all_cis_cbl==1) | (pdata_b$all_cis_runx1==1))

pdata_bi_all = pdata_b[which(pdata_b$is_biallelic==1),]
pdata_bi_tp53 = pdata_b[which(pdata_b$all_cis_tp53==1),]
pdata_bi = pdata_b[which(pdata_b$is_biallelic==1 & pdata_b$all_cis_tp53==0),]
pdata_mono = pdata_b[which(pdata_b$is_biallelic==0 & pdata_b$is_both==1),]
pdata_snv = pdata[which(pdata$is_mut==1 & pdata$is_cna==0),]; pdata_snv$is_biallelic = 0

m_mono = coxph(Surv(fgstart, fgstop, fgstatus) ~ 1, weight=fgwt, data=pdata_mono); sf_mono = survfit(m_mono)
m_bi_all = coxph(Surv(fgstart, fgstop, fgstatus) ~ 1, weight=fgwt, data=pdata_bi_all); sf_bi_all = survfit(m_bi_all)
m_bi_tp53 = coxph(Surv(fgstart, fgstop, fgstatus) ~ 1, weight=fgwt, data=pdata_bi_tp53); sf_bi_tp53 = survfit(m_bi_tp53)
m_bi = coxph(Surv(fgstart, fgstop, fgstatus) ~ 1, weight=fgwt, data=pdata_bi); sf_bi = survfit(m_bi)
m_snv = coxph(Surv(fgstart, fgstop, fgstatus) ~ 1, weight=fgwt, data=pdata_snv); sf_snv = survfit(m_snv)

## draw ED Fig.6c ##
pdf(paste0(path_2,"bi_mono_fg.pdf"),width=7,height=10)
show_di(list(sf_bi_tp53,sf_bi_all,sf_bi,sf_mono,sf_snv,sf_none),">=2 alterations",c("darkgreen","olivedrab3","purple","violet","red","gray"),ymax_df=0.2)
dev.off()

## comparison betwen all biallelic and other combinations ##
pdata_bi_all_mono = rbind(pdata_bi_all,pdata_mono)
m_bi_all_mono = coxph(Surv(fgstart, fgstop, fgstatus) ~ is_biallelic+fac_age+is_male+is_OEE_1_0+is_OEE_1_2, weight=fgwt, data=pdata_bi_all_mono)
summary(m_bi_all_mono)

## comparison betwen biallelic (other than TP53) and SNV+CNA (mono) ##
pdata_bi_mono = rbind(pdata_bi,pdata_mono)
m_bi_mono = coxph(Surv(fgstart, fgstop, fgstatus) ~ is_biallelic + fac_age + is_male + is_OEE_1_0 +  is_OEE_1_2, weight=fgwt, data=pdata_bi_mono)
summary(m_bi_mono)

## comparison between SNV+CNA (mono) vs SNV alone ##
pdata_mono_snv = rbind(pdata_mono,pdata_snv)
m_both_mono_snv = coxph(Surv(fgstart, fgstop, fgstatus) ~ is_both + fac_age + is_male + is_OEE_1_0 +  is_OEE_1_2, weight=fgwt, data=pdata_mono_snv)
summary(m_both_mono_snv)

## comparison between biallelic (other than TP53) and SNV alone ##
pdata_bi_snv = rbind(pdata_bi,pdata_snv)
m_both_bi_snv = coxph(Surv(fgstart, fgstop, fgstatus) ~ is_biallelic + fac_age + is_male + is_OEE_1_0 +  is_OEE_1_2+all_n_alt+fac_size, weight=fgwt, data=pdata_bi_snv)
summary(m_both_bi_snv)



## multi-variate risk prediction model ##
Hb_l_thld = 11
WBC_l_thld = 3000
Plt_l_thld = 10

pdata_2$is_cytopenia = 0
pdata_2$is_cytopenia[which((pdata_2$Hb <= Hb_l_thld) | (pdata_2$WBC <= WBC_l_thld) | (pdata_2$PLT<= Plt_l_thld))] = 1

m_multi <- coxph(Surv(fgstart, fgstop, fgstatus) ~ 
all_n_alt +			# +1 alteration
fac_size +			# +10% cell fraction
is_both +			# coocurrence of SNV and CNA
is_cytopenia +		# presence of cytopenia
fac_age + 			# +10 years old
is_male + 			# male gender
is_OEE_1_0 + 		# examined by OEE v1.0
is_OEE_1_2,			# examined by OEE v1.2
weight=fgwt, data=pdata_2)

## output data for Fig.5a ##
comb_num_size_multi = cbind(summary(m_multi)$coef,summary(m_multi)$conf.int)
write.table(comb_num_size_multi,file=paste0(path_2,"comb_num_size_multi.txt"),col.name=T,row.name=T,sep="\t",quote=F)


## analysis stratified by number of alterations ##
## draw ED Fig.8f-h ##
pdf(paste0(path_2,"both_multi_inc.pdf"),width=15,height=4)
par(mfrow=c(1,3),mai=rep(0.2,4))

for(i in 2:4){
pdata_both = pdata[which(pdata$all_n_alt==i & pdata$is_both==1),]
sf_both = survfit( coxph(Surv(fgstart, fgstop, fgstatus) ~ 1, weight=fgwt, data=pdata_both))

pdata_either_multi = pdata[which((pdata$all_n_mut==i & pdata$all_n_cna==0) | (pdata$all_n_cna==i & pdata$all_n_mut==0)),]
pdata_mut_alone_multi = pdata[which(pdata$all_n_mut==i & pdata$all_n_cna==0),]
pdata_cna_alone_multi = pdata[which(pdata$all_n_cna==i & pdata$all_n_mut==0),]

sf_either_multi = survfit(coxph(Surv(fgstart, fgstop, fgstatus) ~ 1, weight=fgwt, data=pdata_either_multi))
sf_mut_alone_multi = survfit(coxph(Surv(fgstart, fgstop, fgstatus) ~ 1, weight=fgwt, data=pdata_mut_alone_multi))
sf_cna_alone_multi = survfit(coxph(Surv(fgstart, fgstop, fgstatus) ~ 1, weight=fgwt, data=pdata_cna_alone_multi))

pdata_none = pdata[which(pdata$is_none==1),]
m_cna_none <- coxph(Surv(fgstart, fgstop, fgstatus) ~ 1, weight=fgwt, data=pdata_none)
sf_none = survfit(m_none)

show_di(list(sf_both,sf_either_multi,sf_mut_alone_multi,sf_cna_alone_multi,sf_none),
paste0(i," alterations"),c("purple","violet","red","skyblue","gray"),ymax_df=0.12)

pdata_both_either_multi = rbind(pdata_both,pdata_either_multi)
m_both_either_multi <- coxph(Surv(fgstart, fgstop, fgstatus) ~ is_both + fac_age + is_male + is_OEE_1_0 +  is_OEE_1_2, weight=fgwt, data=pdata_both_either_multi)
pdata_both_mut_multi = rbind(pdata_both,pdata_mut_alone_multi)
m_both_mut_multi <- coxph(Surv(fgstart, fgstop, fgstatus) ~ is_both + fac_age + is_male + is_OEE_1_0 +  is_OEE_1_2, weight=fgwt, data=pdata_both_mut_multi)
pdata_both_cna_multi = rbind(pdata_both,pdata_cna_alone_multi)
m_both_cna_multi <- coxph(Surv(fgstart, fgstop, fgstatus) ~ is_both + fac_age + is_male + is_OEE_1_0 +  is_OEE_1_2, weight=fgwt, data=pdata_both_cna_multi)

## output results of comparison ##
print(rbind(comb_all,c(summary(m_both_either_multi)$coef[1,],summary(m_both_either_multi)$conf.int[1,]))

}
dev.off()



## stratification by both number of alterations and cell fractions ##
thld = 0.10*10
pdata_large_1 = pdata[which(pdata$all_n_alt==1 & pdata$fac_size>=thld),]
pdata_large_2 = pdata[which(pdata$all_n_alt>=2 & pdata$fac_size>=thld),]
pdata_large_3 = pdata[which(pdata$all_n_alt>=3 & pdata$fac_size>=thld),]
pdata_small_1 = pdata[which(pdata$all_n_alt==1 & pdata$fac_size<thld),]
pdata_small_2 = pdata[which(pdata$all_n_alt>=2 & pdata$fac_size<thld),]
pdata_small_3 = pdata[which(pdata$all_n_alt>=3 & pdata$fac_size<thld),]
pdata_none = pdata[which(pdata$is_none==1),]
sf_none = survfit(coxph(Surv(fgstart, fgstop, fgstatus) ~ 1, weight=fgwt, data=pdata_none))
sf_large_1 = survfit(coxph(Surv(fgstart, fgstop, fgstatus) ~ 1, weight=fgwt, data=pdata_large_1))
sf_large_2 = survfit(coxph(Surv(fgstart, fgstop, fgstatus) ~ 1, weight=fgwt, data=pdata_large_2))
sf_large_3 = survfit(coxph(Surv(fgstart, fgstop, fgstatus) ~ 1, weight=fgwt, data=pdata_large_3))
sf_small_1 = survfit(coxph(Surv(fgstart, fgstop, fgstatus) ~ 1, weight=fgwt, data=pdata_small_1))
sf_small_2 = survfit(coxph(Surv(fgstart, fgstop, fgstatus) ~ 1, weight=fgwt, data=pdata_small_2))
sf_small_3 = survfit(coxph(Surv(fgstart, fgstop, fgstatus) ~ 1, weight=fgwt, data=pdata_small_3))

## draw Fig.4d ##
pdf(paste0(path_2,"vaf_multi.pdf"),width=5,height=7)
show_di(list(sf_large_1,sf_large_2,sf_large_3,sf_small_1,sf_small_2,sf_small_3,sf_none),"",c("violet","violet","violet","red","red","red","gray"),ymax_df=0.10)
#show_di(list(sf_large_1,sf_large_2,sf_small_1,sf_small_2),"",c("violet","purple","pink","red"),ymax_df=0.06)
dev.off()


## ED Fig.6g ##
pdata_TP53 = pdata[which(pdata$TP53==1 | pdata$is_tp53_ddpcr==1),]
m <- coxph(Surv(fgstart, fgstop, fgstatus) ~ is_tp53_ddpcr+fac_age+is_male+is_OEE_1_0+is_OEE_1_2,weight=fgwt,data=pdata_TP53)
m_coef = cbind(summary(m)$coef,summary(m)$conf.int)
write.table(m_coef,file="tp53_17ploh.txt",col.name=T,row.name=T,sep="\t",quote=F)

pdata_TP53_bi = pdata[which(pdata$is_tp53_ddpcr==1),]
m_TP53_bi <- survfit(coxph(Surv(fgstart, fgstop, fgstatus) ~ 1, weight=fgwt, data=pdata_TP53_bi))
pdata_TP53_mono = pdata[which(pdata$TP53==1 & pdata$is_tp53_ddpcr==0),]
m_TP53_mono <- survfit(coxph(Surv(fgstart, fgstop, fgstatus) ~ 1, weight=fgwt, data=pdata_TP53_mono))

pdf(paste0(path_2,"tp53_17ploh.pdf"),width=5,height=7)
show_di(list(m_TP53_bi,m_TP53_mono),col=c("purple","red"),label="")
dev.off()



## Cumulative mortality by different number of alterations ##
pdata_none = pdata[which(pdata$is_none==1),]
m_none = coxph(Surv(fgstart, fgstop, fgstatus) ~ 1, weight=fgwt, data=pdata_none)
sf_none = survfit(m_none)

## draw Fig.4c ##
sf_list = list()
for(i in 1:4){
if(i<4){
pdata_i = pdata[which(pdata$all_n_alt==i),]
}else if(i==4){
pdata_i = pdata[which(pdata$all_n_alt>=i),]
}
m_i <- coxph(Surv(fgstart, fgstop, fgstatus) ~ 1, weight=fgwt, data=pdata_i)
sf_list[[i]] = survfit(m_i)
}
sf_list[[5]] = sf_none
pdf(paste0(path_2,"cum_inc_by_n_alt.pdf"),width=5,height=6.5)
show_di(sf_list,"",c("firebrick1","firebrick2","firebrick3","firebrick4","gray"),ymax_df=0.06)
dev.off()


## draw ED Fig.8a ##
sf_list = list()
for(i in 1:3){
if(i<3){
pdata_i = pdata[which(pdata$all_n_mut==i),]
}else if(i==3){
pdata_i = pdata[which(pdata$all_n_mut>=i),]
}
m_i <- coxph(Surv(fgstart, fgstop, fgstatus) ~ 1, weight=fgwt, data=pdata_i)
sf_list[[i]] = survfit(m_i)
}
sf_list[[4]] = sf_none
pdf(paste0(path_2,"cum_inc_by_n_mut.pdf"),width=5,height=6.5)
show_di(sf_list,"",c("firebrick2","firebrick3","firebrick4","gray"),ymax_df=0.06)
dev.off()


## draw ED Fig.8b ##
sf_list = list()
for(i in 1:3){
if(i<3){
pdata_i = pdata[which(pdata$all_n_cna==i),]
}else if(i==3){
pdata_i = pdata[which(pdata$all_n_cna>=i),]
}
m_i <- coxph(Surv(fgstart, fgstop, fgstatus) ~ 1, weight=fgwt, data=pdata_i)
sf_list[[i]] = survfit(m_i)
}
sf_list[[4]] = sf_none
pdf(paste0(path_2,"cum_inc_by_n_cna.pdf"),width=5,height=6.5)
show_di(sf_list,"",c("firebrick2","firebrick3","firebrick4","gray"),ymax_df=0.06)
dev.off()


## Comparison among subjects with different numbers of alterations ##
res_n_alt = matrix(nrow=0,ncol=9)
for(i in 1:4){
	if(i<4){
		pdata_i = pdata[which(pdata$is_none==1 | pdata$all_n_alt==i),]
	}else if(i==4){
		pdata_i = pdata[which(pdata$is_none==1 | pdata$all_n_alt>=i),]
	}
	m_i <- coxph(Surv(fgstart, fgstop, fgstatus) ~ is_any + fac_age + is_male + is_OEE_1_0 +  is_OEE_1_2, weight=fgwt, data=pdata_i)
	print(summary(m_i))
	res_n_alt = rbind(c(summary(m_i)$coef[1,],summary(m_i)$conf.int[1,]),res_n_alt)
}
for(i in 1:4){
	if(i<4){
		pdata_i = pdata[which(pdata$all_n_alt==i | pdata$all_n_alt==i+1),]
	}else if(i==4){
		pdata_i = pdata[which(pdata$all_n_alt==i | pdata$all_n_alt>=i+1),]
	}
	pdata_i$all_n_alt[which(pdata_i$all_n_alt==i)] = 0
	pdata_i$all_n_alt[which(pdata_i$all_n_alt>=i+1)] = 1
	m_i <- coxph(Surv(fgstart, fgstop, fgstatus) ~ all_n_alt + fac_age + is_male + is_OEE_1_0 +  is_OEE_1_2, weight=fgwt, data=pdata_i)
	print(summary(m_i))
	res_n_alt = rbind(c(summary(m_i)$coef[1,],summary(m_i)$conf.int[1,]),res_n_alt)
}
for(i in 1:3){
	if(i<3){
		pdata_i = pdata[which(pdata$all_n_alt==i | pdata$all_n_alt==i+2),]
	}else if(i==3){
		pdata_i = pdata[which(pdata$all_n_alt==i | pdata$all_n_alt>=i+2),]
	}
	pdata_i$all_n_alt[which(pdata_i$all_n_alt==i)] = 0
	pdata_i$all_n_alt[which(pdata_i$all_n_alt>=i+2)] = 1
	m_i <- coxph(Surv(fgstart, fgstop, fgstatus) ~ all_n_alt + fac_age + is_male + is_OEE_1_0 +  is_OEE_1_2, weight=fgwt, data=pdata_i)
	print(summary(m_i))
	res_n_alt = rbind(c(summary(m_i)$coef[1,],summary(m_i)$conf.int[1,]),res_n_alt)
}
for(i in 1:2){
	if(i<2){
		pdata_i = pdata[which(pdata$all_n_alt==i | pdata$all_n_alt==i+3),]
	}else if(i==2){
		pdata_i = pdata[which(pdata$all_n_alt==i | pdata$all_n_alt>=i+3),]
	}
	pdata_i$all_n_alt[which(pdata_i$all_n_alt==i)] = 0
	pdata_i$all_n_alt[which(pdata_i$all_n_alt>=i+3)] = 1
	m_i <- coxph(Surv(fgstart, fgstop, fgstatus) ~ all_n_alt + fac_age + is_male + is_OEE_1_0 +  is_OEE_1_2, weight=fgwt, data=pdata_i)
	print(summary(m_i))
	res_n_alt = rbind(c(summary(m_i)$coef[1,],summary(m_i)$conf.int[1,]),res_n_alt)
}
for(i in 1:1){
	pdata_i = pdata[which(pdata$all_n_alt==i | pdata$all_n_alt>=i+4),]
	pdata_i$all_n_alt[which(pdata_i$all_n_alt==i)] = 0
	pdata_i$all_n_alt[which(pdata_i$all_n_alt>=i+4)] = 1
	m_i <- coxph(Surv(fgstart, fgstop, fgstatus) ~ all_n_alt + fac_age + is_male + is_OEE_1_0 +  is_OEE_1_2, weight=fgwt, data=pdata_i)
	print(summary(m_i))
	res_n_alt = rbind(c(summary(m_i)$coef[1,],summary(m_i)$conf.int[1,]),res_n_alt)
}
res_n_alt


## Comparison among subjects with different numbers of SNVs ##
res_n_mut = matrix(nrow=0,ncol=9)
for(i in 1:3){
	if(i<3){
		pdata_i = pdata[which(pdata$is_none==1 | pdata$all_n_mut==i),]
	}else if(i==3){
		pdata_i = pdata[which(pdata$is_none==1 | pdata$all_n_mut>=i),]
	}
	m_i <- coxph(Surv(fgstart, fgstop, fgstatus) ~ is_mut + fac_age + is_male + is_OEE_1_0 +  is_OEE_1_2, weight=fgwt, data=pdata_i)
	print(summary(m_i))
	res_n_mut = rbind(c(summary(m_i)$coef[1,],summary(m_i)$conf.int[1,]),res_n_mut)
}
for(i in 1:3){
	if(i<3){
		pdata_i = pdata[which(pdata$all_n_mut==i | pdata$all_n_mut==i+1),]
	}else if(i==3){
		pdata_i = pdata[which(pdata$all_n_mut==i | pdata$all_n_mut>=i+1),]
	}
	pdata_i$all_n_mut[which(pdata_i$all_n_mut==i)] = 0
	pdata_i$all_n_mut[which(pdata_i$all_n_mut>=i+1)] = 1
	m_i <- coxph(Surv(fgstart, fgstop, fgstatus) ~ all_n_mut + fac_age + is_male + is_OEE_1_0 +  is_OEE_1_2, weight=fgwt, data=pdata_i)
	print(summary(m_i))
	res_n_mut = rbind(c(summary(m_i)$coef[1,],summary(m_i)$conf.int[1,]),res_n_mut)
}
for(i in 1:2){
	if(i<2){
		pdata_i = pdata[which(pdata$all_n_mut==i | pdata$all_n_mut==i+2),]
	}else if(i==2){
		pdata_i = pdata[which(pdata$all_n_mut==i | pdata$all_n_mut>=i+2),]
	}
	pdata_i$all_n_mut[which(pdata_i$all_n_mut==i)] = 0
	pdata_i$all_n_mut[which(pdata_i$all_n_mut>=i+2)] = 1
	m_i <- coxph(Surv(fgstart, fgstop, fgstatus) ~ all_n_mut + fac_age + is_male + is_OEE_1_0 +  is_OEE_1_2, weight=fgwt, data=pdata_i)
	print(summary(m_i))
	res_n_mut = rbind(c(summary(m_i)$coef[1,],summary(m_i)$conf.int[1,]),res_n_mut)
}
for(i in 1:1){
	pdata_i = pdata[which(pdata$all_n_mut==i | pdata$all_n_mut==i+3),]
	pdata_i$all_n_mut[which(pdata_i$all_n_mut==i)] = 0
	pdata_i$all_n_mut[which(pdata_i$all_n_mut>=i+3)] = 1
	m_i <- coxph(Surv(fgstart, fgstop, fgstatus) ~ all_n_mut + fac_age + is_male + is_OEE_1_0 +  is_OEE_1_2, weight=fgwt, data=pdata_i)
	print(summary(m_i))
	res_n_mut = rbind(c(summary(m_i)$coef[1,],summary(m_i)$conf.int[1,]),res_n_mut)
}
res_n_mut


## Comparison among subjects with different numbers of CNAs ##
res_n_cna = matrix(nrow=0,ncol=9)
for(i in 1:3){
	if(i<3){
		pdata_i = pdata[which(pdata$is_none==1 | pdata$all_n_cna==i),]
	}else if(i==3){
		pdata_i = pdata[which(pdata$is_none==1 | pdata$all_n_cna>=i),]
	}
	m_i <- coxph(Surv(fgstart, fgstop, fgstatus) ~ is_cna + fac_age + is_male + is_OEE_1_0 +  is_OEE_1_2, weight=fgwt, data=pdata_i)
	print(summary(m_i))
	res_n_cna = rbind(c(summary(m_i)$coef[1,],summary(m_i)$conf.int[1,]),res_n_cna)
}
for(i in 1:3){
	if(i<3){
		pdata_i = pdata[which(pdata$all_n_cna==i | pdata$all_n_cna==i+1),]
	}else if(i==3){
		pdata_i = pdata[which(pdata$all_n_cna==i | pdata$all_n_cna>=i+1),]
	}
	pdata_i$all_n_cna[which(pdata_i$all_n_cna==i)] = 0
	pdata_i$all_n_cna[which(pdata_i$all_n_cna>=i+1)] = 1
	m_i <- coxph(Surv(fgstart, fgstop, fgstatus) ~ all_n_cna + fac_age + is_male + is_OEE_1_0 +  is_OEE_1_2, weight=fgwt, data=pdata_i)
	print(summary(m_i))
	res_n_cna = rbind(c(summary(m_i)$coef[1,],summary(m_i)$conf.int[1,]),res_n_cna)
}
for(i in 1:2){
	if(i<2){
		pdata_i = pdata[which(pdata$all_n_cna==i | pdata$all_n_cna==i+2),]
	}else if(i==2){
		pdata_i = pdata[which(pdata$all_n_cna==i | pdata$all_n_cna>=i+2),]
	}
	pdata_i$all_n_cna[which(pdata_i$all_n_cna==i)] = 0
	pdata_i$all_n_cna[which(pdata_i$all_n_cna>=i+2)] = 1
	m_i <- coxph(Surv(fgstart, fgstop, fgstatus) ~ all_n_cna + fac_age + is_male + is_OEE_1_0 +  is_OEE_1_2, weight=fgwt, data=pdata_i)
	print(summary(m_i))
	res_n_cna = rbind(c(summary(m_i)$coef[1,],summary(m_i)$conf.int[1,]),res_n_cna)
}
for(i in 1:1){
	pdata_i = pdata[which(pdata$all_n_cna==i | pdata$all_n_cna>=i+3),]
	pdata_i$all_n_cna[which(pdata_i$all_n_cna==i)] = 0
	pdata_i$all_n_cna[which(pdata_i$all_n_cna>=i+3)] = 1
	m_i <- coxph(Surv(fgstart, fgstop, fgstatus) ~ all_n_cna + fac_age + is_male + is_OEE_1_0 +  is_OEE_1_2, weight=fgwt, data=pdata_i)
	print(summary(m_i))
	res_n_cna = rbind(c(summary(m_i)$coef[1,],summary(m_i)$conf.int[1,]),res_n_cna)
}
res_n_cna 


## Comparison: CH(+/-) and CBC(+/-) ##
pdata_ab_chp = pdata_ab[which(pdata_ab$is_any==1),]
pdata_ab_multi_chp = pdata_ab[which(pdata_ab_multi$is_any==1),]
pdata_nm_chp = pdata_nm[which(pdata_nm$is_any==1),]
pdata_ab_chn = pdata_ab[which(pdata_ab$is_any==0),]
pdata_nm_chn = pdata_nm[which(pdata_nm$is_any==0),]

pdata_ab_chp$is_cbc = 1
pdata_ab_chn$is_cbc = 1
pdata_nm_chp$is_cbc = 0
pdata_nm_chn$is_cbc = 0

sf_1 = survfit(coxph(Surv(fgstart, fgstop, fgstatus) ~ 1, weight=fgwt, data=pdata_nm_chn))
sf_2 = survfit(coxph(Surv(fgstart, fgstop, fgstatus) ~ 1, weight=fgwt, data=pdata_ab_chn))
sf_3 = survfit(coxph(Surv(fgstart, fgstop, fgstatus) ~ 1, weight=fgwt, data=pdata_nm_chp))
sf_4 = survfit(coxph(Surv(fgstart, fgstop, fgstatus) ~ 1, weight=fgwt, data=pdata_ab_chp))
sf_5 = survfit(coxph(Surv(fgstart, fgstop, fgstatus) ~ 1, weight=fgwt, data=pdata_ab_multi_chp))

## 10 year mortality of those without CH or CBC abnormality
sf10(sf_1)
## 10 year mortality of those with CH and multi-lineage cytopenia
sf10(sf_5)

## draw Fig.4e ##
pdf(paste0(path_2,"cum_inc.pdf"),width=5,height=6.5)
show_di(list(sf_1,sf_2,sf_3,sf_4),"",c("gray","pink","firebrick1","firebrick3"),ymax_df=0.025)
dev.off()


## Comparison: CH(+/-) and CBC(+/-) ##
## stratified Log rank test ##
## stratified by age and gender ##
cmp_data_org_2_ab$is_cbc = 1
cmp_data_org_2_ab$age_strata = as.numeric(cmp_data_org_2_ab$fac_age<7)　## < or >= median
cmp_data_org_2_ab$HM = as.numeric(cmp_data_org_2_ab$status_HM=="Event")
cmp_data_org_2_ab_chp = cmp_data_org_2_ab[which(cmp_data_org_2_ab$is_any==1),]
cmp_data_org_2_ab_chn = cmp_data_org_2_ab[which(cmp_data_org_2_ab$is_any==0),]

cmp_data_org_2_nm$is_cbc = 0
cmp_data_org_2_nm$age_strata = as.numeric(cmp_data_org_2_nm$fac_age<7)
cmp_data_org_2_nm$HM = as.numeric(cmp_data_org_2_nm$status_HM=="Event")
cmp_data_org_2_nm_chp = cmp_data_org_2_nm[which(cmp_data_org_2_nm$is_any==1),]
cmp_data_org_2_nm_chn = cmp_data_org_2_nm[which(cmp_data_org_2_nm$is_any==0),]

cmp_data_org_2_ab_chp_vs_nm_chn = rbind(cmp_data_org_2_ab_chp,cmp_data_org_2_nm_chn)
cmp_data_org_2_ab_chp_vs_ab_chn = rbind(cmp_data_org_2_ab_chp,cmp_data_org_2_ab_chn)
cmp_data_org_2_ab_chp_vs_nm_chp = rbind(cmp_data_org_2_ab_chp,cmp_data_org_2_nm_chp)
cmp_data_org_2_nm_chp_vs_nm_chn = rbind(cmp_data_org_2_nm_chp,cmp_data_org_2_nm_chn)
cmp_data_org_2_nm_chp_vs_ab_chn = rbind(cmp_data_org_2_nm_chp,cmp_data_org_2_ab_chn)
cmp_data_org_2_ab_chn_vs_nm_chn = rbind(cmp_data_org_2_ab_chn,cmp_data_org_2_nm_chn)

survdiff(Surv(ftime,HM) ~ is_cbc+strata(age_strata)+strata(is_male),data=cmp_data_org_2_ab_chp_vs_nm_chn)
survdiff(Surv(ftime,HM) ~ is_any+strata(age_strata)+strata(is_male),data=cmp_data_org_2_ab_chp_vs_ab_chn)
survdiff(Surv(ftime,HM) ~ is_cbc+strata(age_strata)+strata(is_male),data=cmp_data_org_2_ab_chp_vs_nm_chp)
survdiff(Surv(ftime,HM) ~ is_any+strata(age_strata)+strata(is_male),data=cmp_data_org_2_nm_chp_vs_nm_chn)
survdiff(Surv(ftime,HM) ~ is_any+strata(age_strata)+strata(is_male),data=cmp_data_org_2_nm_chp_vs_ab_chn)
survdiff(Surv(ftime,HM) ~ is_cbc+strata(age_strata)+strata(is_male),data=cmp_data_org_2_ab_chn_vs_nm_chn)


## End of analysis ##
q()
