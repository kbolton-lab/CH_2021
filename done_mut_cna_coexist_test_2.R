
## draw Supplementary Fig.6f ##

genes = c("DNMT3A","TET2","JAK2","TP53","NRAS","EZH2","CBL")
chr = c(2,4,9,17,1,7,11)
starts = c(25455830,106034562,4949510,7563450,115247085,148504464,119076986)
ends = c(25565459,106234240,5163918,7592173,115259515,148581441,119178859)
type = c("LOH","LOH","UPD","LOH","UPD","LOH","LOH")
gene_reg = data.frame(genes,chr,starts,ends,type)

# CH
mut_path = "MatMed_2021_mut_table.txt"
cna_path = "MatMed_2021_cna_table.csv"

mut = read.table(mut_path,header=T,stringsAsFactor=F,sep="\t",quote="",colClass="character")
mut = mut[which(!is.na(mut$Chr)),]
mut$Gene.refGene[which(mut$Gene.refGene=="TET2;TET2")] = "TET2"
mut$Gene.refGene[which(mut$Gene.refGene=="U2AF1;U2AF1L5 ")] = "U2AF1"

cna = read.csv(cna_path,header=T,stringsAsFactor=F,quote="",comment.char="")
cna = cna[which(!is.na(cna$chr)),]

id = as.character(1:11234)

starts = starts/(10^6)
ends = ends/(10^6)
type = c("LOH","LOH","neutral","LOH","neutral","LOH","LOH")


res_coex_all = matrix(ncol=2,nrow=0)
n_coex = rep(NA,4)
d_mut = d_cna = list()
h_mut = h_cna = list()

xisq_list = xisq_null_list = list()
for(g in 1:4){
res = matrix(ncol=2,nrow=length(id))
rownames(res) = id

for(i in 1:length(id)){
	cur_id = id[i]
	cur_mut = mut[which(mut$id == cur_id & mut$Gene.refGene == genes[g]),]
	cur_cna = cna[which(cna$id == cur_id & cna$COPY_CHANGE == "neutral" & cna$chr==chr[g] & (as.numeric(cna$start)-20-ends[g])*(as.numeric(cna$end)+20-starts[g])<=0 ),]
	
	if(nrow(cur_mut) != 0){
		res[i,1] = max(as.numeric(cur_mut$misRate))
	}else{
		res[i,1] = 0
	}
	if(nrow(cur_cna) != 0){
		res[i,2] = max(as.numeric(cur_cna$CELL_FRAC))
	}else{
		res[i,2] = 0
	}
}

tmp = res[which(res[,1] != 0 & res[,2] != 0),]
res_coex_all = rbind(res_coex_all,tmp)
n_coex[g] = nrow(tmp)

h_mut[[g]] = res[which(res[,1]!=0),1]
h_cna[[g]] = res[which(res[,2]!=0),1]
d_mut[[g]] = density(res[which(res[,1]!=0),1],bw=.1,from=0,to=0.5)
d_cna[[g]] = density(res[which(res[,2]!=0),2],bw=.1,from=0,to=0.5)

}



i=4
m = d_mut[[i]]$y
c = d_cna[[i]]$y

h_mut[[i]][which(h_mut[[i]]>0.5)] = 0.50
h_cna[[i]][which(h_cna[[i]]>0.5)] = 0.50
h1 = hist(h_mut[[i]],breaks=seq(0,0.5,length=25))$density #counts
b1 = hist(h_mut[[i]],breaks=seq(0,0.5,length=25))$breaks
h2 = hist(h_cna[[i]],breaks=seq(0,0.5,length=25))$density #counts
b2 = hist(h_cna[[i]],breaks=seq(0,0.5,length=25))$breaks
dev.off()


gene_fac = c(
rep(1,n_coex[1]),
rep(2,n_coex[2]),
rep(3,n_coex[3]),
rep(4,n_coex[4])
)

length(h_mut[[1]])
length(h_cna[[1]])

xi = yi = rep(NA,55)
#xi = yi = rep(NA,n_coex[1])
for(i in 1:length(xi)){
xi[i] = sample(1:512,1,prob=d_mut[[1]]$y)/(512*2) #sample(h_mut[[1]],1)  #
yi[i] = sample(1:512,1,prob=d_cna[[1]]$y)/(512*2) #sample(h_cna[[1]],1)  #
}

act0 = sum(res_coex_all[,1] >= res_coex_all[,2])

p = rep(NA,50)

#ks = c(30,32,34,36,38,40)
#ks = c(1,5,10,15,20,25,30,35,40)
#ks = c(30,33,38)
ks = 20:50
totals = list()
for(k in ks){

act = act0 - k

N = 50000
total = rep(0,N)
for(i in 1:N){

cur_exl = sample(1:sum(n_coex),k,replace=F)
gene_fac_tmp = gene_fac[-cur_exl]

for(j in 1:length(gene_fac_tmp)){
	c = gene_fac_tmp[j]
	xi = sample(h_mut[[c]],1)  #sample(1:512,1,prob=d_mut[[c]]$y)
	yi = sample(h_cna[[c]],1)  #sample(1:512,1,prob=d_cna[[c]]$y)
	total[i] = total[i] + as.numeric(xi >= yi)
}
}
total = total + k
totals[[k]] = total

p[k] = sum(as.numeric(total >= act0))/N
print(k)
print(p[k])
}



pdf("p_plot.pdf")
plot(NULL,NULL,
axes=F,
xlim=c(-5,55),ylim=c(-0.5,5.5),
xlab="",ylab=""
)


#ps = read.table("simulation_p_2.txt",header=F,sep="\t",quote="")
#p = ps[,2]
p[which(is.na(p))] = min(na.omit(p))

p[which(p<10^(-4))] = 10^(-4)
q[which(q<10^(-4))] = 10^(-4)

segments(0,0,0,4.3)
segments(0,0,53,0)

latx = (0:10)*5
for(i in 1:length(latx)){
	segments(latx[i],0,latx[i],-0.1)
	text(latx[i],-0.2,sum(n_coex)-latx[i],cex=.6)
}
latx = 0:50
for(i in 1:length(latx)){
	segments(latx[i],0,latx[i],-0.05)
}
laty = 1:4
for(i in 1:length(laty)){
	segments(0,laty[i],-0.8,laty[i])
	text(-1.5,laty[i],laty[i],cex=.7)
}

for(i in 1:length(p)){
segments(i,-log10(p[i]),i+1,-log10(p[i+1]))
}

segments(0,-log10(0.05),50,-log10(0.05),lty="dashed",col=2)


dev.off()



