######################################################################
#### Consolidation of CRISPR-Cas annotations and gene predictions ####
######################################################################

#### Load packages ####
library(ggpubr)
library(ggplot)
library(ggboxplot)
library(dplyr)
library(plyr)

#### 0 - Technical comparison of detection tools ####

casc <- read.csv("~/Thesis/4Data/CRISPR_comparison/0technical/technical_casc.tsv", header=FALSE)
cd <- read.csv("~/Thesis/4Data/CRISPR_comparison/0technical/technical_cd.tsv", header=FALSE)
ccf <- read.csv("~/Thesis/4Data/CRISPR_comparison/0technical/technical_ccf.tsv", header=FALSE)
names(casc)=c("sample","runtime","sd", "crispr","cas")
names(cd)=c("sample","runtime","sd", "crispr","cas","IE","IF","IC")
names(ccf)=c("sample","runtime","sd", "crispr","cas","IE","IF","IC","other/FP")

shapiro.test(casc$runtime) #p<0.05
shapiro.test(cd$runtime) #p>0.05
shapiro.test(ccf$runtime) #p<0.05

wilcox.test(casc$runtime, cd$runtime) # p = 0.000442
wilcox.test(casc$runtime, ccf$runtime) # p = 0.0004671
wilcox.test(cd$runtime, ccf$runtime) # p = 0.0003539

crispr=cbind(casc[,c(1,4)],cd[,4],ccf[,4])
names(crispr)=c('Sample', 'CASC','CD', 'CCF')
View(crispr) 
cas=cbind(casc[,c(1,5)],cd[,5],ccf[,5])
names(cas)=c('Sample', 'CASC','CD', 'CCF')
View(Cas)


#### 1 - Hybrid vs draft genomes - UZL set ####

uzl_lr <- read.csv("~/Thesis/4Data/FinalData/CCF/47_lr_ccf_results.csv", header=FALSE)
uzl_sr <- read.csv("~/Thesis/4Data/FinalData/CCF/47_sr_ccf_results.csv", header=FALSE)
names(uzl_lr)=c("Sample","crispr", "spacers","cas","IE","IF","IC","other/FP", "type_other", "Cas_genes")
names(uzl_sr)=c("Sample","crispr", "spacers","cas","IE","IF","IC","other/FP", "type_other", "Cas_genes")

#one strain with type U = type IV CRISPR-Cas
uzl_u = which(uzl_lr$type_other == "['U']") 
uzl_lr$U = ifelse(which(uzl_lr$type_other == "['U']"), 1, 0)
uzl_sr$U = ifelse(which(uzl_sr$type_other == "['U']"), 1, 0)
uzl_lr$cas[15] = uzl_lr$cas[15]+1
uzl_sr$cas[15] = uzl_sr$cas[15]+1

uzl_lr$cc = ifelse(uzl_lr$crispr > 0 & uzl_lr$cas > 0, 1, 0)
uzl_sr$cc = ifelse(uzl_sr$crispr > 0 & uzl_sr$cas > 0, 1, 0)

sum(uzl_lr$cc) #30 /47 = 63%
sum(uzl_sr$cc) #30
sum(uzl_lr$crispr > 0) #33 /47 = 70%
sum(uzl_sr$crispr > 0) #33
sum(uzl_lr$cas > 0) #31 /47 = 65%
sum(uzl_sr$cas > 0) #31

which(uzl_lr$cas != uzl_sr$cas) 
which(uzl_lr$spacers != uzl_sr$spacers) 
which(uzl_lr$crispr != uzl_sr$crispr) 

######### 2 - Exploring and curating other collections ##########
#### 2.1 504 UMC strains (draft) ####
umc_sr <- read.csv("~/Thesis/4Data/FinalData/CCF/505_ccf_results.csv", header=FALSE)
names(umc_sr)=c("Sample","crispr", "spacers","cas","IE","IF","IC","other/FP", "type_other", "Cas_genes")
umc_sr=umc_sr[-c(101),] #remove contaminated strain

#other types (except for IIIU)? -> U (2), IA (1), ID (10), IIIA (1)
summary(umc_sr$type_other)

#are they valid?
which(umc_sr$type_other == "['U']") 
umc_sr$Cas_genes[c(40,489)]
# 1) again 4 typeU genes: csf4, csf3, csf2, csf1 (strain also has IF system, of which one gene likely is on other contig)
# 2) only one gene: csf2 -> probably not valid
umc_sr[40, 'U'] = 1
umc_sr[-40, 'U'] = 0
umc_sr$cas[40] = umc_sr$cas[40]+1

which(umc_sr$type_other == "['IA', 'IIIU']") 
umc_sr$Cas_genes[100] #only one gene: casRa (strain contains IE system)

umc_sr_ID=which(umc_sr$type_other %in% c("['ID']", "['IIIU', 'ID']", "['ID', 'IIIU']"))
umc_sr$Cas_genes[umc_sr_ID] #each time only cas3_typeID, alone on contig as well
#so idem as many cas3_typeI false positives

#check for irregularities
umc_check=which(umc_sr$cas > 1) #28 strains
umc_sr_sub=subset(umc_sr, cas>1, select=c('IE', 'IF', 'IC', 'Cas_genes'))
#8 strains have IF split over contigs
umc_sr$IF[umc_check[c(3,4,12,16,22,23,24,26)]] = umc_sr$IF[umc_check[c(3,4,12,16,22,23,24,26)]] -1
umc_sr$cas[umc_check[c(3,4,12,16,22,23,24,26)]] = umc_sr$cas[umc_check[c(3,4,12,16,22,23,24,26)]] -1
#1 strain has IE split over contigs
umc_sr$IE[umc_check[27]] = umc_sr$IE[umc_check[27]] -1
umc_sr$cas[umc_check[27]] = umc_sr$cas[umc_check[27]] -1

umc_sr$cc = ifelse(umc_sr$crispr > 0 & umc_sr$cas > 0, 1, 0)

sum(umc_sr$cc == 1) #285 /504 = 57%
sum(umc_sr$crispr > 0) #349 /504 = 69%
sum(umc_sr$cas > 0) #286 /504 = 57%


#### 2.2 JPP collection (hybrid) ####
jpp_lr <- read.csv("~/Thesis/4Data/FinalData/CCF/19_hist_lr_ccf_results.csv", header=FALSE)
names(jpp_lr)=c("Sample","crispr", "spacers","cas","IE","IF","IC","other/FP", "type_other", "Cas_genes")

###other types (except for IIIU)? -> U (1)
summary(jpp_lr$type_other)
jpp_lr$Cas_genes[which(jpp_lr$type_other == "['U']")] #csf2, csf1, csf4, csf3 
jpp_lr[17, 'U'] = 1
jpp_lr[-17, 'U'] = 0
jpp_lr$cas[17] = jpp_lr$cas[17]+1

jpp_lr$cc = ifelse(jpp_lr$crispr > 0 & jpp_lr$cas > 0, 1, 0)

sum(jpp_lr$cc == 1) #12 /19 = 63%
sum(jpp_lr$crispr > 0) #17 /19 = 89%
sum(jpp_lr$cas > 0) #13 /19 = 68%

#irregularities?
which(jpp_lr$cas > 1) # 3 strains
hist_lr_sub=subset(jpp_lr, cas>1, select=c('IE', 'IF', 'IC', 'Cas_genes'))
#all are okay


#### 2.3 NCBI genomes ####
ncbi <- read.csv("~/Thesis/4Data/FinalData/CCF/4812_ccf_results.csv", header=FALSE)
names(ncbi)=c("sample","crispr", "spacers","cas","IE","IF","IC","other/FP", "type_other", "Cas_genes")

#load other meta-data
ncbi_useful_data <- read.delim("~/Thesis/4Data/SizePlot/ncbi_useful_data.txt")
ncbi_full = merge(ncbi_useful_data, ncbi,by.x=c("RefSeq.FTP"), by.y=c("sample"))
names(ncbi_full)[1] = "Sample"
names(ncbi_full)[17] = "Cas_genes"

###other types (except for IIIU)? -> U (64), IA (42), ID (242), IIIA (2), IIIB (3)
summary(ncbi$type_other)

#are they valid?
ncbi_U=which(ncbi_full$type_other %in% c("['U']","['ID', 'U', 'IIIU']","['U', 'ID', 'IIIU']","['U', 'ID']","['U', 'IIIU', 'ID']",
		"['U', 'IIIU']","['IIIU', 'U']","['ID', 'IIIU', 'U']","['IIIU', 'ID', 'U']")) 
check_U=ncbi_full[ncbi_U, c("Sample","crispr","cas","type_other","other/FP", "Replicons", "Cas_genes")]

#irregularities for: 4 (IV-B), 38+57+59 (only csf2), 53 (extra csf3), 61 (split over contigs) (of ncbi_U)
ncbi_U_valid = ncbi_U[-c(38,57,59)]
ncbi_full[ncbi_U_valid, 'U'] = 1
ncbi_full[-ncbi_U_valid, 'U'] = 0
ncbi_full$cas[ncbi_U_valid] = ncbi_full$cas[ncbi_U_valid]+1
#remark: amount of loci are not counted so it's just the presence/absence (<> other types)

#other irregularities
ncbi_IA=which(ncbi$type_other %in% c("['IA', 'IIIU']","['IA']","['IIIU', 'IA']")) 
ncbi$Cas_genes[ncbi_IA] #only casRa
ncbi_ID=which(ncbi$type_other %in% c("['ID']", "['IIIU', 'ID']", "['ID', 'IIIU']","['ID', 'IIIU', 'U']","['ID', 'U', 'IIIU']",
		"['IIIU', 'ID', 'U']","['U', 'ID', 'IIIU']","['U', 'ID']","['U', 'IIIU', 'ID']"))
ncbi$Cas_genes[ncbi_ID] #only cas3
ncbi_IIIA=which(ncbi$type_other %in% c("['IIIA']","['IIIU', 'IIIA']")) #481 and 2199
ncbi$Cas_genes[c(481,2199)] #only one gene: csm2
ncbi_IIIB=which(ncbi$type_other %in% c("['IIIB']","['IIIB', 'IIIU']","['IIIU', 'IIIB']")) #977, 984 and 1839
ncbi$Cas_genes[c(977,984,1839)] #cmr1 gene in 2 cases, cmr6 in another

ncbi_full$cc = ifelse(ncbi_full$crispr > 0 & ncbi_full$cas > 0, 1, 0)
sum(ncbi_full$cc == 1) #2220 /4812 = 46%
sum(ncbi_full$crispr > 0) #2834 /4812 = 59%
sum(ncbi_full$cas > 0) #2323 /4812 = 48%

#not checked for other irregularities due to draft genomes


#### 2.4 Comparison draft vs complete within NCBI collection ####
#make variable assembly: complete/chromosome level vs. scaffold/contig level
ncbi_full$Assembly=as.factor(ifelse(ncbi_full$Level %in% c('Chromosome','Complete Genome'), 'Complete', 'Draft'))
ncbi_complete=which(ncbi_full$Assembly == 'Complete')

wilcox.test(ncbi_full$crispr[ncbi_complete], ncbi_full$crispr[-ncbi_complete]) #0.8836
wilcox.test(ncbi_full$cas[ncbi_complete], ncbi_full$cas[-ncbi_complete]) #0.8631
wilcox.test(ncbi_full$spacers[ncbi_complete], ncbi_full$spacers[-ncbi_complete]) #0.8884

tab=table(ncbi_full$Assembly, ncbi_full$cc)
chisq.test(tab) #p = 0.454


######### 3 - Full dataset: 5382 genomes + add Size/GC/CDS #########
all_lr <- read.csv("~/Thesis/4Data/FinalData/CCF/66_lr_ccf_results.csv", header=FALSE)
names(all_lr)=c("Sample","crispr", "spacers","cas","IE","IF","IC","other/FP", "type_other", "Cas_genes")
#load workspace from previous script (for size and GC)
ALL_lr = merge(UZL_lr, JPP_lr, all = TRUE) 

all_lr$Strain = all_lr$Sample
all_lr$Assembly = 'Complete'
all_lr = merge(all_lr, ALL_lr[,c(1,17,18)], by='Sample')
all_lr$Total.length = all_lr$Total.length/1000000
names(all_lr)[c(15,16)]=c('GC','Size')

umc_sr$Assembly = 'Draft'
umc_sr$Strain = umc_sr$Sample
#add size and GC
umc_sr= merge(umc_sr, UMC_sr[,c(1,16,17)], by = 'Sample')
umc_sr$Total.length = umc_sr$Total.length/1000000
names(umc_sr)[c(15,16)]=c('Size','GC')

names(ncbi_full)[c(4,5)]=c('Size','GC')

#add number of predicted genes from Prokka results
prokka_47 <- read.delim("~/Thesis/4Data/FinalData/Prokka_reports/47_multiqc_prokka.txt")
prokka_33_505_19 <- read.delim("~/Thesis/4Data/FinalData/Prokka_reports/33+505+19_multiqc_prokka.txt")
names(prokka_47)[7]='Proteins'
names(prokka_33_505_19)[5]='Proteins'
umc_sr = merge(umc_sr, prokka_33_505_19[,c(1,5)], by='Sample', all.x=TRUE)
all_lr = merge(all_lr, prokka_47[,c(1,7)], by='Sample', all.x=TRUE)

x=which(is.na(all_lr$Proteins))
all_lr$Proteins[x] = prokka_33_505_19[c(531:549), c('Proteins')]
ncbi_full$Proteins[which(ncbi_full$Proteins == '-')]=NA

ncbi_full$Proteins = as.numeric(levels(ncbi_full$Proteins))[ncbi_full$Proteins]
all=merge(ncbi_full[,-c(6:7,15:17)], all_lr[,-c(8:10)], by=c('Sample','Strain','Size','GC','Proteins','IE','IF', 'IC', 'crispr','cas','U','cc','spacers','Assembly'), all=TRUE)
all=merge(all, umc_sr[,-c(8:10)], by=c('Sample','Strain','Size','GC','Proteins','IE','IF', 'IC', 'crispr','cas','U','cc','spacers','Assembly'), all=TRUE)
all$cc = as.factor(all$cc)

#### 3.1 Assess CRISPR-Cas in full dataset ####
sum(cc == 1) #2547 /5382 = 47%
all_cc = which(cc == 1)

all$corc=as.factor(ifelse(all$cas > 0, 1, ifelse(all$crispr > 0, 1 ,0)))
all_corc = which(all$corc == 1) #3339 /5382 = 62%

all_crispr = which(all$crispr > 0) #3233 /5382 = 60%
all_cas = which(all$cas > 0) #2653 /5383 = 49%

which(crispr==0 & cas>0) #106 /5382 = 2%
which(crispr>0 & cas==0) #686 /5382 = 13%

summary(spacers) #min 0 max 102
mean(spacers[all_cc]) #41.95
mean(spacers[all_crispr]) #35.56

summary(crispr) #min 0 max 8
mean(crispr[all_cc]) #2.79
mean(crispr[all_crispr]) #2.42

###CRISPR-Cas (sub)types
f=which(IE == 0 & IF>0 & IC == 0 & U == 0 ) 
e=which(IE > 0 & IF == 0 & IC == 0 & U == 0 ) 
c=which(IE == 0 & IF == 0 & IC > 0 & U == 0 ) 
u=which(IE == 0 & IF == 0 & IC == 0 & U > 0 ) 

sum(IE[all_cc] > 0) #529 /2653 = 20%
sum(IF[all_cc] > 0) #2046 /2653 = 77%
sum(IC[all_cc] > 0) #154 /2653 = 6%
sum(U > 0) #64 /2547 = 3%
sum(U[all_cc] > 0) #26 /2547 = 1% #most strains have no array!

#Combinations (% relative to amount of strains with cas)
ef=which(IE>0 & IF>0) #142 /2653 = 5%
fc=which(IC>0 & IF>0) #49 /2653 = 2%
ec=which(IE>0 & IC>0) #4 /2653 = 0.2%
uf=which(U>0 & IF>0) #15
ue=which(U>0 & IE>0) #1

#Define type variables: for partial systems (corc) and complete systems (cc)
all$ctype_corc=ifelse(corc == 0, "None", "Only CRISPR")
all$ctype_corc[f]="IF"
all$ctype_corc[e]="IE"
all$ctype_corc[c]="IC"
all$ctype_corc[u]="IV"
all$ctype_corc[c(ef,fc, ec, uf, ue)]="combination"

all$ctype_corc_order <- factor(all$ctype_corc, levels = c('None','Only CRISPR', 'IE', 'IF', 'IC', 'combination', 'IV'), labels = c('None', 'Only CRISPR','I-E', 'I-F','I-C','Combination', 'IV')

all$ctype=all$ctype_corc
all$ctype[-all_cc]="None" #only consider type of those with CC
all$ctype_order <- factor(all$ctype, levels = c('None', 'IE', 'IF', 'IC', 'combination', 'IV'), labels = c('None','I-E','I-F','I-C','Combination','IV')

detach(genomes)
attach(all)

###some visuals
#pie plot 
summary(as.factor(all$ctype_corc))
x <-  c(109,404,209,1883,48)
labels <-  c('I-C','I-E','combination','I-F','IV')
piepercent<- round(100*x/sum(x), 1)
pie(x, labels = paste(labels, piepercent, '%'), col = brewer.pal(6,'Pastel1'), cex=0.8)

#box plot CRISPR
give.n <- function(x){
  return(c(y = 107, label = length(x)))
  #return(c(y = median(x)*1.05, label = length(x))) 
}
ggboxplot(all[which(all$crispr>0),], x="ctype_corc_order", y="spacers", color="black", fill = 'ctype_corc_order', alpha=0.4, palette = 'jco') + 
  #geom_text(aes(y = max(Size),label = paste0("n = ",length(Size)),vjust = 0))+
  stat_summary(fun.data = give.n,  geom = "text", fun.y = median) +
  theme(legend.position = "none")+ylab("Number of spacers") + xlab('')+
  font("ylab", size=14)+scale_x_discrete(labels= c('Only CRISPR','I-E','I-F','I-C','IV','combination'))

#bar plot CC
set = factor(c('RefSeq','JPP','UMCU','UZL'), levels=c('RefSeq','JPP','UMCU','UZL'))
set = rep(set, each=4)
type = factor(c('CRISPR-Cas', 'Only CRISPR', 'Only Cas', 'None'), levels=c('None','Only Cas','Only CRISPR', 'CRISPR-Cas'))
type = rep(type, 4)
prop = c( 2220/4812, 614/4812, 103/4812, 1875/4812, 17/24,5/24,1/24,1/24, 285/504, 64/504, 1/504, 154/504, 30/47, 3/47, 1/47, 13/47)*100
prop = round(prop, 1)
df = data.frame(set, type, prop)
df_sum = ddply(df, 'set', transform, label_y = cumsum(prop))

ggplot(data = df_sum, aes(x = set, y = prop, fill=type)) + 
  geom_bar(stat="identity", alpha=0.55, width = 0.7, col='black')+
  scale_fill_brewer(palette="RdYlBu")+
  coord_flip()+
  geom_text(aes(y=label_y-prop/2, label=prop), color="black",size=4)+
  theme_classic()+
  ylab("") + xlab("")+
  ylim(0,100)+font("xlab",size=14)+
  theme(legend.position='top',legend.title = element_blank(), axis.text.y=element_text(size = 14), axis.line = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank())



#### 3.2 CRISPR-Cas vs. genome size ####

shapiro.test(sample(Size, 5000)) #test takes max 5000 values
# p < 2.2e-16

#Complete systems (CRISPR and Cas)
wilcox.test(Size~cc) #p < 2.2e-16
wilcox.test(Size~cc, alternative='greater') #idem
mena(all[all_cc, c('Size')]); mean(all[-all_cc, c('Size')])
#6.557128 vs. #6.688394

#Cas alone 
wilcox.test(Size[all_cas],Size[-all_cas]) # p < 2.2e-16
mean(Size[all_cas]);mean(Size[-all_cas]) #6.563951 vs 6.686859

#Crispr alone 
mean(Size[all_crispr]);mean(Size[-all_crispr]) # 6.563714 vs 6.720388
wilcox.test(Size[all_crispr],Size[-all_crispr]) # p < 2.2e-16

#Partial systems (CRISPR or Cas)
wilcox.test(Size[all_corc], Size[-all_corc]) # p < 2.2e-16
mean(Size[all_corc]);mean(Size[-all_corc]) #6.568927 vs 6.719997

#Compare types
compare_means(Size~ctype,data = all,ref.group = "None",method ="wilcox.test")
compare_means(Size~ctype_corc,data = all,ref.group = "None",method ="wilcox.test")

ggboxplot(all, x = "ctype_order", y = "Size", color = "ctype", 
         palette = "jco", fill='ctype', alpha=0.3, xlab="CC type", ylab = "Genome size (Mbp)")+
		stat_compare_means(label = "p.signif", method = "wilcox.test", ref.group = "None") 

ggboxplot(all, x = "ctype_corc_order", y = "Size", color = "ctype_corc", 
         palette = "jco", fill='ctype_corc', alpha=0.3, xlab="CC type", ylab = "Genome size (Mbp)")+
		stat_compare_means(label = "p.signif", method = "wilcox.test", ref.group = "None") 


##Genome size distributions
#1) both CRISPR and cas
ggplot(all, aes(x=Size, color=cc, fill=cc))+
  geom_histogram(aes(y=..density..), alpha=0.5, position="identity")+
  geom_density(alpha=0.4)+
  ylab("Density") + xlab("Genome size (Mbp)") 
  
#without type IV: 
ggplot(all[-c(u,ue,uf),], aes(x=Size, color=cc, fill=cc))+
  geom_histogram(aes(y=..density..), alpha=0.5, position="identity")+
  geom_density(alpha=0.4)+
  ylab("Density") + xlab("Genome size (Mbp)") 

#2) CRISPR or Cas 
ggplot(all, aes(x=Size, color=corc, fill=corc))+
  geom_histogram(aes(y=..density..), alpha=0.5, position="identity")+
  geom_density(alpha=0.4)+
  ylab("Density") + xlab("Genome size (Mbp)") 
  
#3) per type 
#active systems
ggplot(all, aes(x=Size, color=ctype, fill=ctype))+
  geom_histogram(aes(y=..density..), alpha=0.5, position="identity")+
  geom_density(alpha=0.4)+ ylab("Density") + xlab("Genome size (Mbp)") +
  facet_grid(rows=vars(ctype_order))
  
#include partial systems
ggplot(all, aes(x=Size, color=ctype_corc, fill=ctype_corc))+
  geom_histogram(aes(y=..density..), alpha=0.5, position="identity")+
  geom_density(alpha=0.4)+
  facet_grid(rows=vars(ctype_corc_order))+
  ylab("Density") + xlab("Genome size (Mbp)") 


#### 3.3 CRISPR-Cas vs. origin (clinical?) ####

# load data from IPCD site
IPCD_clin = read_excel("Thesis/4Data/Population structure/IPCD_ids.xlsx")
a = merge(all, IPCD_clin, by="Strain", all.x = T)
a$Host[c(4827:4892, 4893:5396)]="Human" #strains from JPP, UZL and UMCU collection are clinical as well
a$Host[which(is.na(a$Host))]='background'
summary(as.factor(a$Host))

# CRISPR-Cas?
a$temp = ifelse(a$cc==1, 'CRISPR and Cas', ifelse(a$cas > 0, 'Only Cas', ifelse(a$crispr >0, 'Only CRISPR', 'None')))
round(prop.table(table(a$Host, a$temp), margin = 1), 2)
round(prop.table(table(a$Host, a$cc), margin = 1), 2)

chisq.test(table(a$Host[which(a$Host %in% c('background','Human'))], a$cc[which(a$Host %in% c('background','Human'))]))
#p = 7.054e-05
chisq.test(table(a$Host[which(a$Host %in% c('background','Human'))], a$temp[which(a$Host %in% c('background','Human'))]))
#p = 8.704e-08

  
#### 3.4 CRISPR-Cas vs. number of proteins ####
summary(all$Proteins) #many outliers... only consider valid ones 

#exclude outliers:
q1 = quantile(Proteins, na.rm = T)[2];q3 = quantile(Proteins,na.rm = T)[4]
outlierup = which(Proteins > (q3+1.5*IQR(Proteins, na.rm = T))) # 5, 1267, 1269, 1280
outlierdown = which(Proteins < (q1-1.5*IQR(Proteins, na.rm = T))) #1688, 4268

# -> 1267, 1269, 1280 are close to max of others, doesn't seem wrong
#others are really too different, although they seem ok in terms of other data...: just exclude them here
outliers=c(5, 1688, 4268) 

summary(Proteins[-outliers]) #mean = 6145

shapiro.test(sample(Proteins[-outliers],5000)) # p < 2.2e-16
wilcox.test(Proteins[-outliers]~cc[-outliers]) # p < 2.2e-16
mean(Proteins[which(cc[-outliers]==1)], na.rm = T) #6084.799
mean(Proteins[which(cc[-outliers]==0)], na.rm = T) #6198.952

compare_means(Proteins~ctype_corc, data=all[-outliers,], ref.group = "None", method="wilcox.test")

ggboxplot(all[-outliers,], x = "ctype_corc_order", y = "Proteins", color = "ctype_corc", 
          palette = "jco", fill='ctype_corc', alpha=0.3, xlab="CC type", ylab = "Number of proteins")+
		stat_compare_means(label = "p.signif", method = "wilcox.test", ref.group = "None") 

ggplot(all[-outliers,], aes(x=Proteins, color=corc, fill=corc))+
  geom_histogram(aes(y=..density..), alpha=0.5, position="identity")+
  geom_density(alpha=0.4)+
  ylab("Density") + xlab("Number of genes") 


#### 3.5 CRISPR-Cas vs. GC content ####
ggplot(all, aes(x=GC, color=corc, fill=corc))+
  geom_histogram(aes(y=..density..), alpha=0.5, position="identity", binwidth = 0.1)+
  geom_density(alpha=0.4)+
  ylab("Density") + xlab("GC content (%)") 


######### 4 - Add population structure data #########
#Difference in size could be due to phylogenetic group

### 4.1 load data from comparative analysis and Freschi et al.(2019)
popstr603 <- read_excel("C:/Users/emmav/OneDrive/Documenten/Thesis/4Data/Population structure/603_phyl_groups.xlsx",col_types = c("text", "text"))
Freschi_SGI <- read.delim("~/Thesis/4Data/Population structure/Freschi_info.txt")
names_together <- read.table("~/Thesis/4Data/Population structure/4812_names_together.txt", colClasses = c("factor","factor"), header = T)

#Freschi data needs converting of IDs via names from IPCD:
IPCD_ids <- read_excel("~/Thesis/4Data/Population structure/IPCD_ids.xlsx", sheet = "filtered", col_types = c("text", "text", "text"))
IPCD_ids$strain = paste('ipcd', IPCD_ids$`Isolate ID`, sep='')
names(IPCD_ids)[2]='Strain'
x=merge(IPCD_ids[,c(2,4)], Freschi_SGI, by='strain', all.y=T)
x=merge(x, names_together, by='strain', all.x = T)
x=merge(x, all[,c(1,2)], by='Strain', all.x = T)

x$Sample=coalesce(x$Sample.x, x$Sample.y)
x=x[complete.cases(x$Sample),c(1,3,6)]
#weirdly 15 samples are duplicated... -> remove manually
n_occur <- data.frame(table(x$Sample))
x[which(x$Sample %in% n_occur[which(n_occur$Freq > 1),1]),]
x=x[-c(325,602,603,608,622,699,700,701,703,704,705,707,710,738),]
x=x[-c(326,328),]

all= merge(all, popstr603, by = 'Strain', all.x = TRUE)
all = merge(all, x[,c(2,3)], by='Sample', all.x=T)
all$Group = coalesce(all$Group, as.character(all$group))
all=all[,-c(22)]
all$Group = as.factor(all$Group)
detach(all);attach(all)

summary(Group)  #imbalanced! 
#1448 group 1, 350 group 2, 17 group 3, 13 group 4, 18 group 5

### 4.2 Phylogenetic group vs. genome size
compare_means(Size~Group, data=all)

give.n <- function(x){
  return(c(y = 7.8, label = length(x)))
  #return(c(y = median(x)*1.05, label = length(x))) 
}

ggboxplot(all[complete.cases(Group),], x="Group", y="Size", color="black", fill = 'Group', alpha=0.4, palette='jco') + 
  stat_summary(fun.data = give.n,  geom = "text", fun.y = median) +
  theme(legend.position = "none")+ylab("Genome size (Mbp)") +
  stat_compare_means(label = "p.format", method = "wilcox.test", comparisons = list(c('1','2'))) 

#Genome size histogram for groups 1 and 2
ggplot(all[which(Group %in%c(1,2)),], aes(x=Size))+#, color=corc, fill=corc))+
  geom_histogram(aes(y=..density..),color='blue', fill='blue', alpha=0.5, position="identity", binwidth = 0.1)+
  geom_density(alpha=0.4)+
  geom_vline(data = ddply(all[which(Group %in%c(1,2)),], c("Group"), summarize, avg = mean(Size)), aes(xintercept=avg),colour="red", linetype="dashed", size=1) +
  facet_grid(rows=vars(Group))+
  ylab("Density") + xlab("Genome size (Mbp)") +
  font("xlab", size=14)+ font("ylab", size=14)+ theme_bw()

### 4.3 Phylogenetic group vs. CRISPR-Cas occurrences
tab1=table(Group, ctype_corc)
round(prop.table(tab1, margin=1)*100, digits=1)
tab2=table(Group[which(all$Group %in% c(1,2) & all$ctype_corc != 'None')], ctype_corc[which(all$Group %in% c(1,2) & all$ctype_corc != 'None')])
round(prop.table(tab2, margin=1)*100, digits=1)
tab3=table(Group, corc)
chisq.test(matrix(c(454,994,136,214), nrow = 2, byrow = T)) #0.008804
round(prop.table(tab3, margin=1)*100, digits=1)

ggplot(all[which(Group %in%c(1,2)),], aes(x=Size), color=corc, fill=corc))+
  geom_histogram(aes(y=..density..),color='blue', fill='blue', alpha=0.5, position="identity", binwidth = 0.1)+
  geom_density(alpha=0.4)+
  facet_grid(rows=vars(Group))+
  ylab("Density") + xlab("Genome size (Mbp)") 

ggplot(all[which(all$Group %in% c(1,2) & all$crispr > 0),], aes(x=Group, y=spacers, fill=Group)) + 
  geom_violin(trim=F, alpha=0.2)+
  geom_jitter(aes(shape=as.factor(ctype_corc_order), color=as.factor(ctype_corc_order)), alpha=1.5, position=position_jitter(0.2), cex = 1.1)+
  theme_bw()+ylab('Spacer count')+font("xlab", size=14)+ font("ylab", size=14)+
  theme(legend.position = 'top')
