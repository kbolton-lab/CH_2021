
## import genomic data ##
source(paste0(path,"ImportGenomicData_HM_OR.R"))
cmp_data_org = cmp_data_org[,setdiff(colnames(cmp_data_org),c("all_time_to_traf","traf_fail"))]
cmp_data_org = na.omit(cmp_data_org)


## select disease ##
cmp_data_org$HM = cmp_data_org$status_LYM


## functions ##
source("show_di.R")
source("writeForest.R")
source("getCnaLabel.R")
makeSf = function(d){ m = coxph(Surv(fgstart, fgstop, fgstatus) ~ 1, weight=fgwt, data=d); sf = survfit(m); return(sf) }


## case cohort, fine-gray ##
pdata = cmp_data_org


## calculate hazard ratio ##
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

res_all = matrix(NA,ncol=8,nrow=length(Names))
rownames(res_all) = Names
for(i in 1:length(Names)){
cur_pdata = pdata[which(pdata[,Names[i]]==1 | pdata$is_none==1),] 
colnames(cur_pdata)[which(colnames(cur_pdata)==Names[i])] = "cur_alt"
n = sum(cur_pdata$cur_alt==1)
if(n<=5) next
m = glm(HM ~ cur_alt+fac_age+is_male+is_OEE_1_0+is_OEE_1_2+is_cancer,family=binomial,data=cur_pdata)
or = exp(summary(m)$coef[2,1])
lw = exp(summary(m)$coef[2,1] - 1.96*summary(m)$coef[2,2])
up = exp(summary(m)$coef[2,1] + 1.96*summary(m)$coef[2,2])
res_all[i,] = c(summary(m)$coef[2,],or,lw,up,n)
}
rownames(res_all) = Labels
colnames(res_all) = c(colnames(summary(m)$coef),"exp(coef)","lower .95","upper .95","n")


## individual SNV ##
res_mut = matrix(ncol=8,nrow=ncol(all_mat_mut_is))
rownames(res_mut) = colnames(all_mat_mut_is)
colnames(res_mut) = c("coef","se_coef","z","p","exp(coef)","lower .95","upper .95","n")

for(i in 1:ncol(all_mat_mut_is)){
cur_alt = colnames(all_mat_mut_is)[i]
cur_pdata = pdata[which(pdata[,cur_alt]==1 | pdata$is_none==1),c("HM",cur_alt,"fac_age","is_male","is_OEE_1_0","is_OEE_1_2","is_cancer")]
colnames(cur_pdata)[which(colnames(cur_pdata)==cur_alt)] = "cur_alt"
n = sum(pdata[,cur_alt] == 1)
if(n<5) next
colnames(cur_pdata) = c("HM","cur_alt","fac_age","is_male","is_OEE_1_0","is_OEE_1_2","is_cancer")
m = glm(HM ~ cur_alt+fac_age+is_male+is_OEE_1_0+is_OEE_1_2+is_cancer,family=binomial,data=cur_pdata)
or = exp(summary(m)$coef[2,1])
lw = exp(summary(m)$coef[2,1] - 1.96*summary(m)$coef[2,2])
up = exp(summary(m)$coef[2,1] + 1.96*summary(m)$coef[2,2])
res_mut[i,] = c(summary(m)$coef[2,],or,lw,up,n)
}


## individual CNA ##
res_cna = matrix(ncol=8,nrow=ncol(all_mat_cna_is))
rownames(res_cna) = colnames(all_mat_cna_is)
colnames(res_cna) = c("coef","se_coef","z","p","exp(coef)","lower .95","upper .95","n")

for(i in 1:ncol(all_mat_cna_is)){
cur_alt = gsub("^","X",gsub("-",".",gsub(":",".",colnames(all_mat_cna_is)[i])))
cur_pdata = pdata[which(pdata[,cur_alt]==1 | pdata$is_none==1),c("HM",cur_alt,"fac_age","is_male","is_OEE_1_0","is_OEE_1_2","is_cancer")]
colnames(cur_pdata)[which(colnames(cur_pdata)==cur_alt)] = "cur_alt"
n = sum(pdata[,cur_alt] == 1)
#if(n<5) next
colnames(cur_pdata) = c("HM","cur_alt","fac_age","is_male","is_OEE_1_0","is_OEE_1_2","is_cancer")
m = glm(HM ~ cur_alt+fac_age+is_male+is_OEE_1_0+is_OEE_1_2+is_cancer,family=binomial,data=cur_pdata)
or = exp(summary(m)$coef[2,1])
lw = exp(summary(m)$coef[2,1] - 1.96*summary(m)$coef[2,2])
up = exp(summary(m)$coef[2,1] + 1.96*summary(m)$coef[2,2])
res_cna[i,] = c(summary(m)$coef[2,],or,lw,up,n)
}


## sort by effect size ##
rownames(res_cna) = sapply(rownames(res_cna),getCnaLabel)
res_mut = res_mut[order(res_mut[,"exp(coef)"],decreasing=T),]
res_cna = res_cna[order(res_cna[,"exp(coef)"],decreasing=T),]


## ED Fig.9a ##
## result of analysis
res_all
res_mut
res_cna


## End of analysis
q()




