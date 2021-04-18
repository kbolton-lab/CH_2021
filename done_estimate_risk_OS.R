
## import genomic data ##
source("ImportGenomicData_CVD_OS.R")
cmp_data_org = cmp_data_org[,setdiff(colnames(cmp_data_org),c("all_time_to_traf","traf_fail"))]
id = rownames(cmp_data_org)
cmp_data_org = na.omit(cmp_data_org)


## select disease ##
cmp_data_org$fail = cmp_data_org$status_ALL


## functions ##
source("show_di.R")
source("writeForest.R")
source("getCnaLabel.R")
check = function(x){length(unique(x$id))}


## case cohort, fine-gray ##
source("prepareCaseCohortFineGray.R")
pdata = cmp_data_org
pdata = pdata[which(pdata$is_subcohort==1),]


## plot overall survivals ##
## calculation of hazard ratios ##

## OS in non-CH subjects ##
pdata_none = pdata[which(pdata$is_none==1),]
m_none = coxph(Surv(ftime,fail=="Event") ~ 1,data=pdata_none)
sf_none = survfit(m_none)


## Effeect of SNV on OS ##
## Draw ED Fig.10a ##
pdata_mut_none = pdata[which(pdata$is_mut==0 & pdata$is_mut_large==0),]
sf_mut_none = survfit(coxph(Surv(ftime,fail=="Event") ~ 1,data=pdata_mut_none))

pdata_mut_small = pdata[which(pdata$is_mut==1 & pdata$is_mut_large==0),]
sf_mut_small = survfit(coxph(Surv(ftime,fail=="Event") ~ 1,data=pdata_mut_small))

pdata_mut_large = pdata[which(pdata$is_mut==1 & pdata$is_mut_large==1),]
sf_mut_large = survfit(coxph(Surv(ftime,fail=="Event") ~ 1,data=pdata_mut_large))

pdf(paste0(path_2,"cum_inc_snv_small_large_fg.pdf"),width=4,height=5)
par(mfrow=c(1,1),mai=rep(0.2,4))
show_di_2(list(sf_mut_large,sf_mut_small,sf_mut_none),"SNV",c("red","pink","gray"))
dev.off()

## result of comparison ##
## mut small vs. none
## mut large vs. none
pdata_mut_small = pdata[which(pdata$is_mut_large==0),]
m_mut_small = coxph(Surv(ftime,fail=="Event") ~ is_mut+fac_age+is_male+is_OEE_1_0+is_OEE_1_2+all_BMI+all_HT+all_DM+all_HL+all_Smoke+all_Drink,data=pdata_mut_small)
pdata_mut_large = pdata[which(pdata$is_mut_large==1 | pdata$is_mut==0),]
m_mut_large = coxph(Surv(ftime,fail=="Event") ~ is_mut_large+fac_age+is_male+is_OEE_1_0+is_OEE_1_2+all_BMI+all_HT+all_DM+all_HL+all_Smoke+all_Drink,data=pdata_mut_large)
m_res = rbind(
c(summary(m_mut_small)$coef["is_mut",],summary(m_mut_small)$conf.int["is_mut",]),
c(summary(m_mut_large)$coef["is_mut_large",],summary(m_mut_large)$conf.int["is_mut_large",]))
rownames(m_res) = c("VAF<5% vs no SNV","VAF>=5% vs no SNV")
m_res


## Effect of CNA on OS ##
## Draw ED Fig.10b ##
pdata_cna_none = pdata[which(pdata$is_cna==0 & pdata$is_cna_large==0),]
sf_cna_none = survfit(coxph(Surv(ftime,fail=="Event") ~ 1,data=pdata_cna_none))
pdata_cna_small = pdata[which(pdata$is_cna==1 & pdata$is_cna_large==0),]
sf_cna_small = survfit(coxph(Surv(ftime,fail=="Event") ~ 1,data=pdata_cna_small))
pdata_cna_large = pdata[which(pdata$is_cna==1 & pdata$is_cna_large==1),]
sf_cna_large = survfit(coxph(Surv(ftime,fail=="Event") ~ 1,data=pdata_cna_large))
pdf(paste0(path_2,"cum_inc_cna_small_large_fg.pdf"),width=4,height=5)
par(mfrow=c(1,1),mai=rep(0.2,4))
show_di_2(list(sf_cna_large,sf_cna_small,sf_cna_none),"CNA",c("red","pink","gray"),ymax_df=0.25)
dev.off()

## result of comparison ##
## mut small vs. none
## mut large vs. none
pdata_cna_small = pdata[which(pdata$is_cna_large==0),]
m_cna_small = coxph(Surv(ftime,fail=="Event") ~ is_cna+fac_age+is_male+is_OEE_1_0+is_OEE_1_2+all_BMI+all_HT+all_DM+all_HL+all_Smoke+all_Drink,data=pdata_cna_small)
pdata_cna_large = pdata[which(pdata$is_cna_large==1 | pdata$is_cna==0),]
m_cna_large = coxph(Surv(ftime,fail=="Event") ~ is_cna_large+fac_age+is_male+is_OEE_1_0+is_OEE_1_2+all_BMI+all_HT+all_DM+all_HL+all_Smoke+all_Drink,data=pdata_cna_large)
m_res = rbind(
c(summary(m_cna_small)$coef["is_cna",],summary(m_cna_small)$conf.int["is_cna",]),
c(summary(m_cna_large)$coef["is_cna_large",],summary(m_cna_large)$conf.int["is_cna_large",]))
rownames(m_res) = c("Cell fraction<5% vs no CNA","Cell fraction>=5% vs no CNA")
m_res


## Combined effect of SNV and CNA on OS ##
## Draw ED Fig.10c ##
pdata_both = pdata[which(pdata$is_mut_large==1 & pdata$is_cna==1),]
m_both = coxph(Surv(ftime,fail=="Event") ~ 1,data=pdata_both)
sf_both = survfit(m_both)
pdata_mut_alone = pdata[which(pdata$is_mut_large==1 & pdata$is_cna==0),]
m_mut_alone = coxph(Surv(ftime,fail=="Event") ~ 1,data=pdata_mut_alone)
sf_mut_alone = survfit(m_mut_alone)
pdata_cna_alone = pdata[which(pdata$is_mut_large==0 & pdata$is_cna==1),]
m_cna_alone = coxph(Surv(ftime,fail=="Event") ~ 1,data=pdata_cna_alone)
sf_cna_alone = survfit(m_cna_alone)
pdata_none = pdata[which(pdata$is_none==1),]
sf_none = survfit(coxph(Surv(ftime,fail=="Event") ~ 1,data=pdata_none))
pdf(paste0(path_2,"both_cum_inc_fg.pdf"),width=4,height=5)
par(mfrow=c(1,1),mai =rep(0.2,4))
show_di_2(list(sf_both,sf_mut_alone,sf_cna_alone,sf_none),"Fine-Gray",c("purple","red","skyblue","gray"))
dev.off()

## result of comparison ##
## Both vs. SNV alone
## Both vs. CNA alone
## Both vs. Either
pdata_both_mut_alone = rbind(pdata_both,pdata_mut_alone)
m_both_mut_alone = coxph(Surv(ftime,fail=="Event") ~ is_both+fac_age+is_male+is_OEE_1_0+is_OEE_1_2+all_BMI+all_HT+all_DM+all_HL+all_Smoke+all_Drink,data=pdata_both_mut_alone)
pdata_both_cna_alone = rbind(pdata_both,pdata_cna_alone)
m_both_cna_alone = coxph(Surv(ftime,fail=="Event") ~ is_both+fac_age+is_male+is_OEE_1_0+is_OEE_1_2+all_BMI+all_HT+all_DM+all_HL+all_Smoke+all_Drink,data=pdata_both_cna_alone)
pdata_both_either = rbind(pdata_both,pdata_mut_alone,pdata_cna_alone)
m_both_either = coxph(Surv(ftime,fail=="Event") ~ is_both+fac_age+is_male+is_OEE_1_0+is_OEE_1_2+all_BMI+all_HT+all_DM+all_HL+all_Smoke+all_Drink,data=pdata_both_either)
m_res = rbind(
c(summary(m_both_mut_alone)$coef["is_both",],summary(m_both_mut_alone)$conf.int["is_both",]),
c(summary(m_both_cna_alone)$coef["is_both",],summary(m_both_cna_alone)$conf.int["is_both",]),
c(summary(m_both_either)$coef["is_both",],summary(m_both_either)$conf.int["is_both",]))
rownames(m_res) = c("both vs mut (VAF>5%)","both vs cna","both vs either")
m_res



## prob. bi-allele vs prob. mono-allele ##
## determine prob. bi-allelic cases ##
pdata_2 = pdata
pdata_2$is_biallelic = as.numeric(
(pdata_2$all_cis_jak2==1) | 
(pdata_2$all_cis_tp53==1) | 
(pdata_2$all_cis_dnmt3a==1) | 
(pdata_2$all_cis_tet2 == 1) | 
(pdata_2$all_cis_ezh2==1) | 
(pdata_2$all_cis_cbl==1) | 
(pdata_2$all_cis_runx1==1))

## stratify subjects ##
pdata_bi = pdata_2[which(pdata_2$is_biallelic==1),]
pdata_mono = pdata_2[which(pdata_2$is_biallelic==0 & pdata_2$is_both==1),]
pdata_snv = pdata[which(pdata$is_mut==1 & pdata$is_cna==0),]
pdata_snv$is_biallelic = 0

## Draw ED Fig.6e ##
pdata_none = pdata_2[which(pdata_2$is_none==1),]
m_none = coxph(Surv(ftime,fail=="Event") ~ 1,data=pdata_none); sf_none = survfit(m_none)
m_mono <- coxph(Surv(ftime,fail=="Event") ~ 1,data=pdata_mono); sf_mono = survfit(m_mono)
m_bi <- coxph(Surv(ftime,fail=="Event") ~ 1,data=pdata_bi); sf_bi = survfit(m_bi)
m_snv <- coxph(Surv(ftime,fail=="Event") ~ 1,data=pdata_snv); sf_snv = survfit(m_snv)
pdf(paste0(path_2,"bi_mono_fg.pdf"),width=7,height=10)
show_di_2(list(sf_bi_all,sf_mono,sf_snv,sf_none),">=2 alterations",c("purple","violet","red","gray"),ymax_df=1)
dev.off()

## perform comparison ##
pdata_bi_mono = rbind(pdata_bi,pdata_mono)
m_bi_mono <- coxph(Surv(ftime,fail=="Event") ~ is_biallelic+fac_age+is_male+is_OEE_1_0+is_OEE_1_2+all_BMI+all_HT+all_DM+all_HL+all_Smoke+all_Drink,data=pdata_bi_mono)
pdata_mono_snv = rbind(pdata_mono,pdata_snv)
m_both_mono_snv <- coxph(Surv(ftime,fail=="Event") ~ is_both+fac_age+is_male+is_OEE_1_0+is_OEE_1_2+all_BMI+all_HT+all_DM+all_HL+all_Smoke+all_Drink,data=pdata_mono_snv)
pdata_bi_snv = rbind(pdata_bi,pdata_snv)
m_both_bi_snv <- coxph(Surv(ftime,fail=="Event") ~ is_biallelic+fac_age+is_male+is_OEE_1_0+is_OEE_1_2+all_BMI+all_HT+all_DM+all_HL+all_Smoke+all_Drink,data=pdata_bi_snv)
pdata_mono_snv = rbind(pdata_mono,pdata_snv)
m_both_mono_snv <- coxph(Surv(ftime,fail=="Event") ~ is_both+fac_age+is_male+is_OEE_1_0+is_OEE_1_2+all_BMI+all_HT+all_DM+all_HL+all_Smoke+all_Drink,data=pdata_mono_snv)

## Result of comparison ##
## bi vs. mono
## both(bi) vs. snv alone
## both(mono) vs. snv alone
m_res = rbind(
c(summary(m_bi_mono)$coef["is_both",],summary(m_bi_mono)$conf.int["is_both",]),
c(summary(m_both_mono_snv)$coef["is_both",],summary(m_both_mono_snv)$conf.int["is_both",]),
c(summary(m_both_bi_snv)$coef["is_both",],summary(m_both_bi_snv)$conf.int["is_both",]))
rownames(m_res) = c("bi vs. mono","both(mono) vs. snv alone","both(bi) vs. snv alone")
m_res


## End of analysis ##
q()

