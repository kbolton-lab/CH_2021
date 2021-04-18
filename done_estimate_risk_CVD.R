
## import genomic data ##

source("ImportGenomicData_CVD_OS.R")
cmp_data_org = cmp_data_org[,setdiff(colnames(cmp_data_org),c("all_time_to_traf","traf_fail"))]
cmp_data_org = na.omit(cmp_data_org)
id = rownames(cmp_data_org)
cmp_data_org = cbind(id,cmp_data_org)


## select disease ##

cmp_data_org$fail = cmp_data_org$status_CVD


## functions ##

source("show_di.R")
source("writeForest.R")
source("getCnaLabel.R")


## cohort study, fine-gray ##

source("prepareCaseCohortFineGray.R")
pdata = fg(cmp_data_org)
pdata = pdata[which(pdata$is_subcohort==1),]


## Cumulative mortality in no-CH cases ##
pdata_none = pdata[which(pdata$is_none==1),]
m_none = coxph(Surv(fgstart, fgstop, fgstatus) ~ 1, weight=fgwt, data=pdata_none)
sf_none = survfit(m_none)


## Effect of SNV ##
## draw Fig.6a ##
pdata_mut_none = pdata[which(pdata$is_mut==0 & pdata$is_mut_large==0),]
sf_mut_none = survfit(coxph(Surv(fgstart, fgstop, fgstatus) ~ 1, weight=fgwt, data=pdata_mut_none))
pdata_mut_small = pdata[which(pdata$is_mut==1 & pdata$is_mut_large==0),]
sf_mut_small = survfit(coxph(Surv(fgstart, fgstop, fgstatus) ~ 1, weight=fgwt, data=pdata_mut_small))
pdata_mut_large = pdata[which(pdata$is_mut==1 & pdata$is_mut_large==1),]
sf_mut_large = survfit(coxph(Surv(fgstart, fgstop, fgstatus) ~ 1, weight=fgwt, data=pdata_mut_large))
pdf("cum_inc_snv_small_large_fg.pdf",width=4,height=5)
par(mfrow=c(1,1),mai=rep(0.2,4))
show_di(list(sf_mut_large,sf_mut_small,sf_mut_none),"",c("red","pink","gray"))
dev.off()

## result of comparison ##
## SNV(small) vs no CH
## SNV(large) vs. no CH
pdata_mut_small = pdata[which(pdata$is_mut_large==0),]
m_mut_small = coxph(Surv(fgstart, fgstop, fgstatus) ~ is_mut+fac_age+is_male+is_OEE_1_0+is_OEE_1_2+all_BMI+all_HT+all_DM+all_HL+all_Smoke+all_Drink, weight=fgwt, data=pdata_mut_small)
summary(m_mut_small)
pdata_mut_large = pdata[which(pdata$is_mut_large==1 | pdata$is_mut==0),]
m_mut_large = coxph(Surv(fgstart, fgstop, fgstatus) ~ is_mut_large+fac_age+is_male+is_OEE_1_0+is_OEE_1_2+all_BMI+all_HT+all_DM+all_HL+all_Smoke+all_Drink, weight=fgwt, data=pdata_mut_large)
summary(m_mut_large)
m_res = rbind(
c(summary(m_mut_small)$coef["is_mut",],summary(m_mut_small)$conf.int["is_mut",]),
c(summary(m_mut_large)$coef["is_mut_large",],summary(m_mut_large)$conf.int["is_mut_large",]))
rownames(m_res) = c("VAF<5% vs no SNV","VAF>=5% vs no SNV")
m_res


## Effect of CNA ##
## draw Fig.6c ##
pdata_cna_none = pdata[which(pdata$is_cna==0 & pdata$is_cna_large==0),]
sf_cna_none = survfit(coxph(Surv(fgstart, fgstop, fgstatus) ~ 1, weight=fgwt, data=pdata_cna_none))
pdata_cna_small = pdata[which(pdata$is_cna==1 & pdata$is_cna_large==0),]
sf_cna_small = survfit(coxph(Surv(fgstart, fgstop, fgstatus) ~ 1, weight=fgwt, data=pdata_cna_small))
pdata_cna_large = pdata[which(pdata$is_cna==1 & pdata$is_cna_large==1),]
sf_cna_large = survfit(coxph(Surv(fgstart, fgstop, fgstatus) ~ 1, weight=fgwt, data=pdata_cna_large))
pdf("cum_inc_cna_small_large_fg.pdf",width=4,height=5)
par(mfrow=c(1,1),mai=rep(0.2,4))
show_di(list(sf_cna_large,sf_cna_small,sf_cna_none),"",c("red","pink","gray"))
dev.off()

## result of comparison ##
## CNA(small) vs no CH
## CNA(large) vs. no CH
pdata_cna_small = pdata[which(pdata$is_cna_large==0),]
m_cna_small = coxph(Surv(fgstart, fgstop, fgstatus) ~ is_cna+fac_age+is_male+is_OEE_1_0+is_OEE_1_2+all_BMI+all_HT+all_DM+all_HL+all_Smoke+all_Drink, weight=fgwt, data=pdata_cna_small)
summary(m_cna_small)
pdata_cna_large = pdata[which(pdata$is_cna_large==1 | pdata$is_cna==0),]
m_cna_large = coxph(Surv(fgstart, fgstop, fgstatus) ~ is_cna_large+fac_age+is_male+is_OEE_1_0+is_OEE_1_2+all_BMI+all_HT+all_DM+all_HL+all_Smoke+all_Drink, weight=fgwt, data=pdata_cna_large)
summary(m_cna_large)
m_res = rbind(
c(summary(m_cna_small)$coef["is_cna",],summary(m_cna_small)$conf.int["is_cna",]),
c(summary(m_cna_large)$coef["is_cna_large",],summary(m_cna_large)$conf.int["is_cna_large",]))
rownames(m_res) = c("VAF<5% vs no SNV","VAF>=5% vs no SNV")
m_res


## Effect of combined SNV & CNA ##
## draw Fig.6c ##
pdata_both = pdata[which(pdata$is_mut_large==1 & pdata$is_cna==1),]
sf_both = survfit(coxph(Surv(fgstart, fgstop, fgstatus) ~ 1, weight=fgwt, data=pdata_both))
pdata_mut_alone = pdata[which(pdata$is_mut_large==1 & pdata$is_cna==0),]
sf_mut_alone = survfit(coxph(Surv(fgstart, fgstop, fgstatus) ~ 1, weight=fgwt, data=pdata_mut_alone))
pdata_cna_alone = pdata[which(pdata$is_mut_large==0 & pdata$is_cna==1),]
sf_cna_alone = survfit(coxph(Surv(fgstart, fgstop, fgstatus) ~ 1, weight=fgwt, data=pdata_cna_alone))
pdata_none = pdata[which(pdata$is_none==1),]
sf_none = survfit(coxph(Surv(fgstart, fgstop, fgstatus) ~ 1, weight=fgwt, data=pdata_none))
pdf("both_cum_inc_fg_vaf_large.pdf",width=4,height=5)
par(mfrow=c(1,1),mai =rep(0.2,4))
show_di(list(sf_both,sf_mut_alone,sf_cna_alone,sf_none),"Fine-Gray",c("purple","red","skyblue","gray"),ymax_df=0.40)
dev.off()

## result of comparison ##
## Both vs SNV(>5%) alone
## Both vs. CNA alone
## Both vs. Either
pdata_both = pdata[which(pdata$is_mut_large==1 & pdata$is_cna==1),]
pdata_mut_alone = pdata[which(pdata$is_mut_large==1 & pdata$is_cna==0),]
pdata_both_mut_alone = rbind(pdata_both,pdata_mut_alone)
m_both_mut_alone = coxph(Surv(fgstart, fgstop, fgstatus) ~ is_cna+fac_age+is_male+is_OEE_1_0+is_OEE_1_2+all_BMI+all_HT+all_DM+all_HL+all_Smoke+all_Drink+fac_vaf, weight=fgwt, data=pdata_both_mut_alone)
pdata_both = pdata[which(pdata$is_mut_large==1 & pdata$is_cna==1),]
pdata_cna_alone = pdata[which(pdata$is_mut_large==0 & pdata$is_cna==1),]
pdata_both_cna_alone = rbind(pdata_both,pdata_cna_alone)
m_both_cna_alone = coxph(Surv(fgstart, fgstop, fgstatus) ~ is_mut_large+fac_age+is_male+is_OEE_1_0+is_OEE_1_2+all_BMI+all_HT+all_DM+all_HL+all_Smoke+all_Drink, weight=fgwt, data=pdata_both_cna_alone)
pdata_both = pdata[which(pdata$is_mut_large==1 & pdata$is_cna==1),]
pdata_mut_alone = pdata[which(pdata$is_mut_large==1 & pdata$is_cna==0),]
pdata_cna_alone = pdata[which(pdata$is_mut_large==0 & pdata$is_cna==1),]
pdata_both_either = rbind(pdata_both,pdata_mut_alone,pdata_cna_alone)
pdata_both_either$is_both[which(pdata_both_either$is_mut_large!=1 | pdata_both_either$is_cna!=1)] = 0
pdata_both_either$is_both[which(pdata_both_either$is_mut_large==1 & pdata_both_either$is_cna==1)] = 1
m_both_either = coxph(Surv(fgstart, fgstop, fgstatus) ~ is_both+fac_age+is_male+is_OEE_1_0+is_OEE_1_2+all_BMI+all_HT+all_DM+all_HL+all_Smoke+all_Drink, weight=fgwt, data=pdata_both_either)
m_res = rbind(
c(summary(m_both_mut_alone)$coef["is_cna",],summary(m_both_mut_alone)$conf.int["is_cna",]),
c(summary(m_both_cna_alone)$coef["is_mut_large",],summary(m_both_cna_alone)$conf.int["is_mut_large",]),
c(summary(m_both_either)$coef["is_both",],summary(m_both_either)$conf.int["is_both",])
)
rownames(m_res) = c("both vs mut (VAF>5%)","both vs cna","both vs either")
m_res


### prob. bi-allele vs prob. mono-allele ###

## define bi-allelic cases ##
pdata_2 = pdata
pdata_2$is_biallelic = as.numeric(
(pdata_2$all_cis_jak2==1) | 
(pdata_2$all_cis_tp53==1) | 
(pdata_2$all_cis_dnmt3a==1) | 
(pdata_2$all_cis_tet2 == 1) | 
(pdata_2$all_cis_ezh2==1) | 
(pdata_2$all_cis_cbl==1) | 
(pdata_2$all_cis_runx1==1))

## divide into categories ##
pdata_bi_all = pdata_2[which(pdata_2$is_biallelic==1),]
pdata_mono = pdata_2[which(pdata_2$is_biallelic==0 & pdata_2$is_both==1 & pdata_2$is_mut_large==1),]
pdata_snv = pdata[which(pdata$is_mut_large==1 & pdata$is_cna==0),]
pdata_snv$is_biallelic = 0

## Draw ED Fig.6d ##
pdata_none = pdata[which(pdata$is_none==1),]
m_none = coxph(Surv(fgstart, fgstop, fgstatus) ~ 1, weight=fgwt, data=pdata_none)
sf_none = survfit(m_none)
m_mono <- coxph(Surv(fgstart, fgstop, fgstatus) ~ 1, weight=fgwt, data=pdata_mono)
sf_mono = survfit(m_mono)
m_bi_all <- coxph(Surv(fgstart, fgstop, fgstatus) ~ 1, weight=fgwt, data=pdata_bi_all)
sf_bi_all = survfit(m_bi_all)
m_snv <- coxph(Surv(fgstart, fgstop, fgstatus) ~ 1, weight=fgwt, data=pdata_snv)
sf_snv = survfit(m_snv)
pdf("bi_mono_fg.pdf",width=7,height=10)
show_di(list(sf_bi_all,sf_mono,sf_snv,sf_none),">=2 alterations",c("purple","violet","red","gray"),ymax_df=0.5)
dev.off()

## Comparison among categories ##
pdata_bi_mono = rbind(pdata_bi,pdata_mono)
m_bi_mono <- coxph(Surv(fgstart, fgstop, fgstatus) ~ is_biallelic + fac_age + is_male + is_OEE_1_0 +  is_OEE_1_2, weight=fgwt, data=pdata_bi_mono)
summary(m_bi_mono)
pdata_mono_snv = rbind(pdata_mono,pdata_snv)
m_both_mono_snv <- coxph(Surv(fgstart, fgstop, fgstatus) ~ is_both + fac_age + is_male + is_OEE_1_0 +  is_OEE_1_2, weight=fgwt, data=pdata_mono_snv)
summary(m_both_mono_snv)
pdata_bi_snv = rbind(pdata_bi,pdata_snv)
m_both_bi_snv <- coxph(Surv(fgstart, fgstop, fgstatus) ~ is_biallelic + fac_age + is_male + is_OEE_1_0 +  is_OEE_1_2, weight=fgwt, data=pdata_bi_snv)
summary(m_both_bi_snv)


## both vs either, stratified by n alt ##

## Draw ED Fig.10d-f ##
pdf("both_multi_inc_mut_large.pdf",width=13,height=5)
par(mfrow=c(1,3),mai=rep(0.2,4))

pdata_mut = pdata[which(pdata$all_n_alt==1 & pdata$is_both==0 & pdata$is_mut_large==1),]
sf_mut = survfit( coxph(Surv(fgstart, fgstop, fgstatus) ~ 1, weight=fgwt, data=pdata_mut))
pdata_both = pdata[which(pdata$all_n_alt>=2 & pdata$is_both==1 & pdata$is_mut_large==1),]
sf_both = survfit( coxph(Surv(fgstart, fgstop, fgstatus) ~ 1, weight=fgwt, data=pdata_both))
pdata_mut_alone_multi = pdata[which(pdata$all_n_alt>=2 & pdata$all_n_cna==0 & pdata$is_mut_large==1),]
sf_mut_alone_multi = survfit(coxph(Surv(fgstart, fgstop, fgstatus) ~ 1, weight=fgwt, data=pdata_mut_alone_multi))
pdata_none = pdata[which(pdata$is_none==1),]
sf_none = survfit( coxph(Surv(fgstart, fgstop, fgstatus) ~ 1, weight=fgwt, data=pdata_none))

show_di(list(sf_both,sf_mut_alone_multi,sf_mut,sf_none),
paste0(">=2 alterations"),c("purple","red","orange","gray"),ymax_df=0.35)

pdata_both_mut_multi = rbind(pdata_both,pdata_mut_alone_multi)
m_both_mut_multi <- coxph(Surv(fgstart, fgstop, fgstatus) ~ is_both+fac_age+is_male+is_OEE_1_0+is_OEE_1_2+all_BMI+all_HT+all_DM+all_HL+all_Smoke+all_Drink,weight=fgwt,data=pdata_both_mut_multi)
pdata_both_mut = rbind(pdata_both,pdata_mut)
m_both_mut <- coxph(Surv(fgstart, fgstop, fgstatus) ~ is_both+fac_age+is_male+is_OEE_1_0+is_OEE_1_2+all_BMI+all_HT+all_DM+all_HL+all_Smoke+all_Drink,weight=fgwt,data=pdata_both_mut)

comb_all = rbind(
c(summary(m_both_mut_multi)$coef[1,],summary(m_both_mut_multi)$conf.int[1,]),
c(summary(m_both_mut)$coef[1,],summary(m_both_mut)$conf.int[1,])
)


for(i in 2:4){
pdata_mut = pdata[which(pdata$all_n_alt==1 & pdata$is_both==0 & pdata$is_mut_large==1),]
sf_mut = survfit( coxph(Surv(fgstart, fgstop, fgstatus) ~ 1, weight=fgwt, data=pdata_mut))
pdata_both = pdata[which(pdata$all_n_alt==i & pdata$is_both==1 & pdata$is_mut_large==1),]
sf_both = survfit( coxph(Surv(fgstart, fgstop, fgstatus) ~ 1, weight=fgwt, data=pdata_both))
pdata_mut_alone_multi = pdata[which(pdata$all_n_alt>=i & pdata$all_n_cna==0 & pdata$is_mut_large==1),]
sf_mut_alone_multi = survfit(coxph(Surv(fgstart, fgstop, fgstatus) ~ 1, weight=fgwt, data=pdata_mut_alone_multi))
pdata_none = pdata[which(pdata$is_none==1),]
sf_none = survfit( coxph(Surv(fgstart, fgstop, fgstatus) ~ 1, weight=fgwt, data=pdata_none))

show_di(list(sf_both,sf_mut_alone_multi,sf_mut,sf_none),
paste0(i," alterations"),c("purple","red","orange","gray"),ymax_df=0.35)

pdata_both_mut_multi = rbind(pdata_both,pdata_mut_alone_multi)
m_both_mut_multi <- coxph(Surv(fgstart, fgstop, fgstatus) ~ is_both+fac_age+is_male+is_OEE_1_0+is_OEE_1_2+all_BMI+all_HT+all_DM+all_HL+all_Smoke+all_Drink,weight=fgwt,data=pdata_both_mut_multi)
pdata_both_mut = rbind(pdata_both,pdata_mut)
m_both_mut <- coxph(Surv(fgstart, fgstop, fgstatus) ~ is_both+fac_age+is_male+is_OEE_1_0+is_OEE_1_2+all_BMI+all_HT+all_DM+all_HL+all_Smoke+all_Drink,weight=fgwt,data=pdata_both_mut)

comb_all = rbind(comb_all,
c(summary(m_both_mut_multi)$coef[1,],summary(m_both_mut_multi)$conf.int[1,]),
c(summary(m_both_mut)$coef[1,],summary(m_both_mut)$conf.int[1,])
)
}
dev.off()

## result of comparison ##
## Both vs. SNV(VAF>5%) alone
comb_all


## number of alterations ##

## draw ED Fig.10g ##
sf_list = list()
for(i in 1:3){
if(i<3){
pdata_i = pdata[which(pdata$all_n_alt==i),]
}else if(i==3){
pdata_i = pdata[which(pdata$all_n_alt>=i),]
}
m_i <- coxph(Surv(fgstart, fgstop, fgstatus) ~ 1, weight=fgwt, data=pdata_i)
sf_list[[i]] = survfit(m_i)
}
sf_list[[4]] = sf_none
pdf("cum_inc_by_n_alt.pdf",width=5,height=6.5)
show_di(sf_list,"",c("firebrick1","firebrick3","firebrick4","gray"),ymax_df=0.30)
dev.off()

## Cmparison among subjects with different numbers of alterations ##
res_n_alt = matrix(nrow=0,ncol=9)
for(i in 1:3){
	if(i<3){
		pdata_i = pdata[which(pdata$all_n_alt==0 | pdata$all_n_alt==i),]
	}else if(i==3){
		pdata_i = pdata[which(pdata$all_n_alt==0 | pdata$all_n_alt>=i),]
	}
	if(!cvd) m_i <- coxph(Surv(fgstart, fgstop, fgstatus) ~ is_any+fac_age+is_male+is_OEE_1_0+is_OEE_1_2,weight=fgwt,data=pdata_i)
	if(cvd) m_i <- coxph(Surv(fgstart, fgstop, fgstatus) ~ is_any+fac_age+is_male+is_OEE_1_0+is_OEE_1_2+all_BMI+all_HT+all_DM+all_HL+all_Smoke+all_Drink,weight=fgwt,data=pdata_i)
	
	print(summary(m_i))
	res_n_alt = rbind(c(summary(m_i)$coef[1,],summary(m_i)$conf.int[1,]),res_n_alt)
}
for(i in 0:2){
	if(i<2){
		pdata_i = pdata[which(pdata$all_n_alt==i | pdata$all_n_alt==i+1),]
	}else if(i==2){
		pdata_i = pdata[which(pdata$all_n_alt==i | pdata$all_n_alt>=i+1),]
	}
	pdata_i$all_n_alt[which(pdata_i$all_n_alt==i)] = 0
	pdata_i$all_n_alt[which(pdata_i$all_n_alt>=i+1)] = 1
	if(!cvd) m_i <- coxph(Surv(fgstart, fgstop, fgstatus) ~ all_n_alt+fac_age+is_male+is_OEE_1_0+is_OEE_1_2,weight=fgwt,data=pdata_i)
	if(cvd) m_i <- coxph(Surv(fgstart, fgstop, fgstatus) ~ all_n_alt+fac_age+is_male+is_OEE_1_0+is_OEE_1_2+all_BMI+all_HT+all_DM+all_HL+all_Smoke+all_Drink,weight=fgwt,data=pdata_i)

	print(summary(m_i))
	res_n_alt = rbind(c(summary(m_i)$coef[1,],summary(m_i)$conf.int[1,]),res_n_alt)
}
for(i in 0:1){
	if(i<1){
		pdata_i = pdata[which(pdata$all_n_alt==i | pdata$all_n_alt==i+2),]
	}else if(i==1){
		pdata_i = pdata[which(pdata$all_n_alt==i | pdata$all_n_alt>=i+2),]
	}
	pdata_i$all_n_alt[which(pdata_i$all_n_alt==i)] = 0
	pdata_i$all_n_alt[which(pdata_i$all_n_alt>=i+2)] = 1
	if(!cvd) m_i <- coxph(Surv(fgstart, fgstop, fgstatus) ~ all_n_alt+fac_age+is_male+is_OEE_1_0+is_OEE_1_2,weight=fgwt,data=pdata_i)
	if(cvd) m_i <- coxph(Surv(fgstart, fgstop, fgstatus) ~ all_n_alt+fac_age+is_male+is_OEE_1_0+is_OEE_1_2+all_BMI+all_HT+all_DM+all_HL+all_Smoke+all_Drink,weight=fgwt,data=pdata_i)

	print(summary(m_i))
	res_n_alt = rbind(c(summary(m_i)$coef[1,],summary(m_i)$conf.int[1,]),res_n_alt)
}
for(i in 0:0){
	if(i<0){
		pdata_i = pdata[which(pdata$all_n_alt==i | pdata$all_n_alt==i+3),]
	}else if(i==0){
		pdata_i = pdata[which(pdata$all_n_alt==i | pdata$all_n_alt>=i+3),]
	}
	pdata_i$all_n_alt[which(pdata_i$all_n_alt==i)] = 0
	pdata_i$all_n_alt[which(pdata_i$all_n_alt>=i+3)] = 1
	if(!cvd) m_i <- coxph(Surv(fgstart, fgstop, fgstatus) ~ all_n_alt+fac_age+is_male+is_OEE_1_0+is_OEE_1_2,weight=fgwt,data=pdata_i)
	if(cvd) m_i <- coxph(Surv(fgstart, fgstop, fgstatus) ~ all_n_alt+fac_age+is_male+is_OEE_1_0+is_OEE_1_2+all_BMI+all_HT+all_DM+all_HL+all_Smoke+all_Drink,weight=fgwt,data=pdata_i)
	print(summary(m_i))
	res_n_alt = rbind(c(summary(m_i)$coef[1,],summary(m_i)$conf.int[1,]),res_n_alt)
}

## result of comparison ##
res_n_alt


## analysis end ##
q()
