##################################################
#### Consolidation of anti-CRISPR annotations ####
##################################################

#### 0 - Load packages and data ####
library(ggplot2)
library(ggboxplot)
library(ggpubr)
library(gridExtra)

#get workspace from CRISPR-Cas.R ('all' dataframe in particular)

#load AcrFinder results
acr_all <- read.csv("~/Thesis/4Data/FinalData/Acr/all_acr_results.csv", header=FALSE)
names(acr_all)=c('Sample','acrIE', 'acrIF', 'Acr_genes_hom','nb_hits_hom', 'nb_other_hom', 'gbaLC', 'gbaMC_IE', 'gbaMC_IF', 'gbaMC_IC', 'gbaHC_IE', 'gbaHC_IF', 'gbaHC_IC', 'phage')
acr_all=acr_all[-c(4913),] #remove non-Pa strain

all=merge(all, acr_all, by=c('Sample'), all = TRUE)

#detach(all)
attach(all)

#### 1 - Homology statistics ####
summary(acr_all)[,c(2,3,5,6)]
#nb_hits_hom (and nb_other_hom): max 24 (16) in PaLo287 -> correct? 8 times same locus with one IE and 2 IF
summary(acr_all[,4])
#most abundant locus is combination of AcrIE1 and AcrIF12, second and third also combination of IE and IF

H_hit=which(nb_hits_hom > 0) #2251 /5382 = 42 %
sum(nb_hits_hom[all_cc] > 0) #989 /2251 = 44%
                              #989 /5382 = 19%
sum(nb_hits_hom[-all_cc] > 0) #1262 /2251 =56%
                              #1262 /5382 = 23%
#more hits in strains without active systems

sum(nb_hits_hom[all_corc] >0) #1348 /5382 = 25%
sum(nb_hits_hom[-all_corc] >0) #903 /5382 = 17%
#more hits in strains with partial systems

sum(acrIF > 0 & acrIE > 0) #1842

#Other genes in the loci (~Aca)
sum(nb_other_hom > 0) #1831 / 2251

#pie plot of genes per type
x <-  c(1931,2694,231)
labels <-  c('AcrIE','AcrIF','AcrIE4-IF7') #+1 AcrVA5...
piepercent<- round(100*x/sum(x), 1)
pie(x, labels = paste(labels,piepercent, '%'), col = brewer.pal(3, 'Pastel1'))


#### 2 - GBA statistics ####
summary(acr_all[,c(7:13)])
# no HC I-C hits

LC=which(gbaLC != 0) #772 /5382 = 14%
MC=which(gbaMC_IE > 0 | gbaMC_IF > 0 | gbaMC_IC > 0) #737 /5382 = 14%
HC=which(gbaHC_IE > 0 | gbaHC_IF > 0) #139 /5382 = 3%

which(gbaMC_IE > 0) #131 /5382 = 2%
which(gbaHC_IE > 0) #96 /5382 = 2%
which(gbaMC_IF > 0 ) #580 /5382 = 11%
which(gbaHC_IF > 0) #43 /5382 = 0.8%
which(gbaMC_IC > 0) #47 /5382 = 0.9%

#all GBA are per definition in strains with crispr arrays

all$GBA=0
all$GBA[c(LC,MC, HC)]=1
all$GBA_woLC=0
all$GBA_woLC[c(MC, HC)]=1
detach(all);attach(all)

which(cc == 0 & GBA == 1) #4
which(cc == 1 & GBA == 1) #1596

#pie plot of amount of hits per confidence level
x <-  c(2920,2266,320)
labels <-  c('Low','Medium','High')
piepercent<- round(100*x/sum(x), 1)
gba=pie(x, labels = paste(labels, piepercent, '%'), col = brewer.pal(3, 'Pastel2'))
title(main = 'Confidence levels of GBA hits' )
legend("topright", labels, cex = 0.8, fill = rainbow(length(x)))

#summarizing barplot
type=c('CRISPR-Cas','CRISPR','Cas','None')
type=factor(type, levels=c('None','CRISPR','Cas','CRISPR-Cas'))
type=rep(type,2)
Strategy=as.factor(c(rep('Homology',4),rep('GBA',4)))
prop=round(c(989/2547, 314/686, 45/106,903/2043,1596/2547,4/686,0,0)*100,2)
df=data.frame(type, prop, Strategy)
ggplot(data = df, aes(x = type, y = prop, fill=Strategy)) + 
  geom_bar(stat="identity", position=position_dodge(), alpha=0.8)+
  #scale_fill_brewer(palette="RdYlBu")+
  geom_text(aes(y=(prop+2.5), label=prop), color="black",size=3.5, position=position_dodge(0.8))+
  theme_classic()+
  ylim(0,70) + ylab("Proportion of strains with hits (%)") + xlab("")+
  font("xlab",size=14)+theme(legend.position = 'top', legend.title = element_blank())


#### 3 - Draft vs. complete assemblies ####

## in 47 UZL strains
acr_UZL_sr <- read.csv("~/Thesis/4Data/Anticrispr/Results/acr_47_sr.csv", header=FALSE)
names(acr_UZL_sr)=c('Sample','acrIE', 'acrIF', 'Acr_genes_hom','nb_hits_hom', 'nb_other_hom', 'gbaLC', 'gbaMC_IE', 'gbaMC_IF', 'gbaMC_IC', 'gbaHC_IE', 'gbaHC_IF', 'gbaHC_IC', 'phage')
acr_UZL_lr=all[-c(1:4812,4821,4823,4825,4827,4832:4835,4838,4839,4843,4845,4854,4846,4848,4852,4854,4860,4865,4867,4879:5382), c(1,22:34)]

#inspection: homology same, GBA is mixed up in confidence levels (which makes sense: contigs make that distances to STS cannot be assessed)
#phage hits same

#within RefSeq genomes?
x=table(all$GBA[c(1:4812)], all$Assembly[c(1:4812)])
chisq.test(x) #0.3179
wilcox.test(all$nb_hits_hom ~ all$Assembly) # p = 0.04881 
x=table((all$nb_hits_hom[c(1:4812)]>0), all$Assembly[c(1:4812)])
chisq.test(x) #0.05293


#### 4 - Combine Acr and CRISPR-Cas to obtain corrected annotations  ####

all$acr = ifelse(GBA == 1, 1, ifelse(nb_hits_hom > 0, 1, 0))
#without LC GBA hits (after inspection of size histogram)
all$acr_woLC = 0
all$acr_woLC[c(MC,HC,H_hit)]=1
detach(all);attach(all)

## 1) add homology hits to complete systems (CC)
IE_adj=which(ctype=='IE' & acrIE > 0) 
IF_adj=which(ctype=='IF' & acrIF > 0) 

which(acrIE[ec] > 0) #0
a=which(acrIE[ef] > 0) 
b=which(acrIF[fc] > 0) 
c=which(acrIF[ef] > 0) 
ac=unique(c(a,c)) 
combo_adj=unique(c(fc[b],ef[ac])) 

#define corrected CC as cc_adj
all$cc_adj=cc
all$cc_adj[IE_adj]=0
all$cc_adj[IF_adj]=0
all$cc_adj[combo_adj]=0
#define corrected type annotations as ctype_adj
all$ctype_adj=ctype
all$ctype_adj[IE_adj]="Inactivated"
all$ctype_adj[IF_adj]="Inactivated"
all$ctype_adj[combo_adj] = "Inactivated" 

detach(all);attach(all)

## 2) also add GBA hits (medium and high confidence) to CC
IE_adj_2=which(ctype_adj=='IE' & (gbaHC_IE > 0 | gbaMC_IE > 0)) 
IF_adj_2=which(ctype_adj=='IF' & (gbaHC_IF > 0 | gbaMC_IF > 0)) 
IC_adj_2=which(ctype_adj=='IC' & gbaMC_IC > 0) 

which(gbaHC_IE[ec] > 0 | gbaMC_IE[ec] > 0) 
d=which(gbaHC_IE[ef] > 0 | gbaMC_IE[ef] > 0) 
e=which(gbaHC_IF[fc] > 0 | gbaMC_IF[fc] > 0) 
f=which(gbaHC_IF[ef] > 0 | gbaMC_IF[ef] > 0) 
df=unique(c(d,f)) 
g=which(gbaMC_IC[fc]>0) 
eg=unique(c(e,g)) 
h=which(gbaMC_IC[ec]>0)
combo_adj2=unique(c(fc[eg],ef[df], ec[h])) 

#define cc_adj2
all$cc_adj2=cc_adj
all$cc_adj2[IE_adj_2]=0
all$cc_adj2[IF_adj_2]=0
all$cc_adj2[IC_adj_2]=0
all$cc_adj2[combo_adj2]=0
#define ctype_adj2
all$ctype_adj2=ctype_adj
all$ctype_adj2[IE_adj_2]="Inactivated"
all$ctype_adj2[IF_adj_2]="Inactivated"
all$ctype_adj2[IC_adj_2]="Inactivated"
all$ctype_adj2[combo_adj2]="Inactivated"

all$ctype_adj2_order <- factor(all$ctype_adj2, levels = c('None', 'IE', 'IF', 'IC', 'combination', 'IV', 'Inactivated'))

## 3) consider low confidence GBA hits
rest=which(cc_adj2 == 1 & gbaLC > 1) #603 strains
hist(Size)
hist(Size[rest]) 
# -> don't include

## 4) now corrected CORC (partial systems included) for homology and GBA (without LC)
IE_adj3=which(ctype_corc=='IE' & (acrIE > 0 | gbaHC_IE > 0 | gbaMC_IE > 0)) 
IF_adj3=which(ctype_corc=='IF' & (acrIF > 0 | gbaHC_IF > 0 | gbaMC_IF > 0)) 
IC_adj3=which(ctype_corc=='IC' & gbaMC_IC > 0) 

#!! only CRISPR strains: only possible to do naïve correction
CR_adj3 = which(ctype_corc == 'Only CRISPR' & acr_woLC == 1) #314 = almost half of all only CRISPR strains

combo_adj3=unique(c(combo_adj2, combo_adj)) 

#define corrected CORC as corc_adj
all$corc_adj=corc
all$corc_adj[IE_adj3]=0
all$corc_adj[IF_adj3]=0
all$corc_adj[IC_adj3]=0
all$corc_adj[combo_adj3]=0
all$corc_adj[CR_adj3]=0
#define corrected CORC types as ctype_corc_adj
all$ctype_corc_adj=ctype_corc
all$ctype_corc_adj[IE_adj3]="Inactivated"
all$ctype_corc_adj[IF_adj3]="Inactivated"
all$ctype_corc_adj[IC_adj3]="Inactivated"
all$ctype_corc_adj[CR_adj3]="Inactivated"
all$ctype_corc_adj[combo_adj3] = "Inactivated"
all$ctype_corc_adj[which(ctype_corc=="None")]="None"

all$ctype_corc_adj_order <- factor(ctype_corc_adj, levels = c('None','Inactivated','Only CRISPR', 'IE', 'IF', 'IC','IV', 'combination'), labels = c('None','Inactivated','Only CRISPR','I-E','I-F','I-C','IV','Combination'))

detach(all);attach(all)

#remove type IV strains (large genomes due to plasmid)
all_woIV = all
all_woIV = all_woIV[-c(which(all_woIV$U > 0)),] #64 strains removed
#detach(all);attach(all_woIV)


#### 5 - Revisit genome size (distributions) ####

wilcox.test(Size[which(acr==1)], Size[which(acr==0)]) #<2.2e-16
mean(Size[which(acr==1)]) - mean(Size[which(acr==0)]) #0.088 Mbp
wilcox.test(Size[which(acr_woLC==1)], Size[which(acr_woLC==0)]) #<2.2e-16
mean(Size[which(acr_woLC==1)]) - mean(Size[which(acr_woLC==0)]) #0.157 Mbp

#1) cc_adj (only homology)
corr1=ggplot(all, aes(x=Size, color=as.factor(cc_adj), fill=as.factor(cc_adj)))+
  geom_histogram(aes(y=..density..), alpha=0.5, position="identity")+
  geom_density(alpha=0.4)+
  ylab("Density") + xlab("Genome size (Mbp)") 

#peak in smaller genomes is again bigger...?
#what is in this peak?
#check separately those that had cc in first place:
peak=which(corc==1 & cc_adj==0 & Size < 6.6) #601 strains
summary(as.factor(ctype_corc[peak])) #500 IF (83%), 63 IE, 4 IC and 34 (5.6%) combination

#2) cc_adj2 and ctype_adj2 (GBA as well)
corr2=ggplot(all, aes(x=Size, color=as.factor(cc_adj2), fill=as.factor(cc_adj2)))+
   geom_histogram(aes(y=..density..), alpha=0.5, position="identity")+
   geom_density(alpha=0.4)+
   ylab("Density") + xlab("Genome size (Mbp)") 

grid.arrange(corr1, corr2,ncol=2)

pertype_nice=ggplot(all, aes(x=Size, color=as.factor(ctype_adj2_order), fill=as.factor(ctype_adj2_order)))+
   geom_histogram(aes(y=..density..), alpha=0.5, position="identity")+
   geom_density(alpha=0.4)+
   facet_grid(rows=vars(ctype_order))+
   ylab("Density") + xlab("Genome size (Mbp)") 

compare_means(Size~ctype,data = all,ref.group = "None",method ="wilcox.test") #no acr
compare_means(Size~ctype_adj,data = all,ref.group = "None",method ="wilcox.test") #only hom
compare_means(Size~ctype_adj2,data = all,ref.group = "None",method ="wilcox.test") # also GBA (no LC)
#p-values get more significant with the corrections


#3) corc_adj (hom & gba)
ggplot(all, aes(x=Size, color=as.factor(ctype_corc_adj_order), fill=as.factor(ctype_corc_adj_order)))+
  geom_histogram(aes(y=..density..), alpha=0.5, position="identity")+
  geom_density(alpha=0.4)+
  facet_grid(rows=vars(ctype_corc_order))+
    ylab("Density") + xlab("Genome size (Mbp)") + theme_bw()
    
pertype_corc= ggplot(all, aes(x=Size, color=as.factor(ctype_corc_adj), fill=as.factor(ctype_corc_adj)))+
  geom_histogram(aes(y=..density..), alpha=0.5, position="identity")+
  geom_density(alpha=0.4)+
  facet_grid(rows=vars(ctype_corc_order))+
    ylab("Density") + xlab("Genome size (Mbp)") 
  
#remarkedly: type I-F does not show difference between inactivated and active systems!
#see also boxplot:

give.n <- function(x){
  return(c(y = 7.7, label = length(x)))
  #return(c(y = median(x)*1.05, label = length(x))) 
}
boxplot_final=ggboxplot(all, x = "ctype_corc_adj_order", y = "Size", color = "black", 
    palette = "jco", fill='ctype_corc', alpha=0.3, xlab="CRISPR-Cas type", ylab = "Genome size (Mbp)")+
	stat_compare_means(label = "p.format", method = "wilcox.test", comparisons =list(c('Inactivated','IF')))+
	stat_summary(fun.data = give.n,  geom = "text", fun.y = median) +
	theme(legend.position = "none")+  font("xlab",size=14)+font("ylab",size=14)+
	scale_x_discrete(labels= c('None','Inactivated','Only CRISPR','I-E','I-F','I-C','IV','combination'))

#however: there's still a significant difference 
mean(Size[which(ctype_corc == 'IF' & ctype_corc_adj == 'Inactivated')]); mean(Size[which(ctype_corc_adj == 'IF')])
#6.56 and 6.50
wilcox.test(Size[which(ctype_corc_adj == 'IF' & ctype_corc_adj == 'Inactivated')], Size[which(ctype_corc_adj == 'IF')])
#p = 8.363e-07 the difference is very significant still


#### 6 - Revisit phylogenetic groups ####

tab=table(Group, ctype_corc_adj)
round(prop.table(tab2, margin=1)*100, digits=1)

ggplot(all[which(Group %in%c(1,2)),], aes(x=Size, color=corc_adj, fill=corc_adj))+
  geom_histogram(aes(y=..density..), alpha=0.5, position="identity")+
  geom_density(alpha=0.4)+
  facet_grid(rows=vars(Group))+
  ylab("Density") + xlab("Genome size (Mbp)") +
  font("xlab", size=14)+ font("ylab", size=14)
#only look at group 1 and 2 because too few data for the others
