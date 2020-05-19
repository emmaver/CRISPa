########################################################################
#### Filter NCBI genomes, size distribution and assembly assessment ####
########################################################################

#### Load packages ####
library(ggplot2)
library(ggboxplot)
library(ggpubr)


#### Filter NCBI genomes for outlying GC% and size ####
setwd("C:/Users/emmav/OneDrive/Documenten")

#load genomes_list downloaded from NCBI
genomes <- read.delim("~/Thesis/4Data/NCBI/genomes_list_10-02-20_adapted.txt")

attach(genomes)

#GC filtering
summary(GC.) #especially low GC outliers
q1 = quantile(GC.)[2]; q3 = quantile(GC.)[4]
highgc = which(GC. > (q3+1.5*IQR(GC.))) 
lowgc = which(GC. < (q1-1.5*IQR(GC.))) 

#size filtering
summary(Size..Mb.) #both high and low size outliers
q1 = quantile(Size..Mb.)[2]; q3 = quantile(Size..Mb.)[4]
big = which(Size..Mb. > (q3+1.5*IQR(Size..Mb.))) 
small = which(Size..Mb. < (q1-1.5*IQR(Size..Mb.))) 

#resulting tsv
nopass = unique(append(big, c(small, lowgc, highgc)))
genomes[nopass, c('GC.', 'Size..Mb.','Level','RefSeq.FTP')]

pass_genomes = genomes[-nopass,]
noref=which(pass_genomes$RefSeq.FTP=="-")
ref_genomes = pass_genomes[-noref,]
write.table(ref_genomes[,"RefSeq.FTP"], file='pass_genomes_refseq.txt', col.names = F, row.names = F, quote = F)


#### Genome size distribution plot of all P. aeruginosa genomes ####

#NCBI collection
ncbi=ref_genomes[-which(rownames(ref_genomes) %in% c('265','1085')),] 
#filter out GCF_900167195.1_IMG-taxon_2568526011_annotated_assembly and GCF_900168095.1_IMG-taxon_2606217733_annotated_assembly
ncbi_df=subset(ncbi, select=c('RefSeq.FTP','Strain','Level','Size..Mb.','GC.','Replicons','Genes','Proteins'))
write.table(ncbi_df, "ncbi_useful_data.txt", row.names = F, col.names = T, sep="\t") # to use in further analyses

#UZL collection (hybrid assemblies)
UZL_lr=read.delim("~/Thesis/4Data/SizePlot/UZL_multiqc_quast_lr.txt")
UZL=UZL_lr$Total.length/1000000

#UMC collection (draft assemblies)
UMC_sr=read.delim("~/Thesis/4Data/SizePlot/UMC_multiqc_quast_sr.txt")
UMC=UMC_sr$Total.length/1000000

#extreme values for length -> contamination (Kraken analysis)
UMC=UMC[-c(102,109,3,291,302,465)]
UMC_sr = UMC_sr[-c(102,109,3,291,302,465),]

#JPP collection (hybrid assemblies)
JPP_lr=read.delim("~/Thesis/4Data/SizePlot/JPP_multiqc_quast_lr.txt")
JPP=JPP_lr$Total.length/1000000

total=append(ncbi$Size..Mb., c(UZL,UMC,JPP))
mean(total) #6.626273

#make plot of genome size 
df=data.frame(total)
ggplot(df, aes(x=total))+
  geom_histogram(aes(y=..density..),colour='blue', fill='blue', alpha=0.5, binwidth = 0.04)+
  geom_vline(aes(xintercept=mean(total)),colour="red", linetype="dashed", size=1)+
  annotate("text", x=6.85, y=1.5, label= "Average = 6.63 Mbp", colour="red",fontface=2 ) +
  geom_density(colour="black")+
  ylab("Density") + xlab("Genome size (Mbp)") + #ggtitle("Genome length distribution of 5382 P. aeruginosa strains") +
  theme(plot.title = element_text(hjust = 0.5))+theme_bw()+
  font("xlab", size=14)+ font("ylab", size=14)

#make plot of GC %
total_gc=append(ncbi$GC.,c(UMC_sr$GC...., UZL_lr$GC...., JPP_lr$GC....))
summary(total_gc) #avg 66.33 
dfgc=data.frame(total_gc)
ggplot(dfgc, aes(x=total_gc))+
  geom_histogram(aes(y=..density..),colour='blue', fill='blue', alpha=0.5, binwidth = 0.04)+
  geom_vline(aes(xintercept=mean(total_gc)),colour="red", linetype="dashed", size=1)+
  annotate("text", x=66.20, y=2, label= "Average = 66.33 %", colour="red",fontface=2 ) +
  geom_density(colour="black")+
  ylab("Density") + xlab("GC content (%)") + ggtitle("GC content of 5382 P. aeruginosa strains") +
  theme(plot.title = element_text(hjust = 0.5))+theme_bw()+
  font("xlab", size=14)+ font("ylab", size=14)

#### Hybrid vs. short-read assembly ####

x=c(rep('Hybrid', 47), rep('Short-read', 504))
y=c(UZL_lr$N50, UMC_sr$N50)
z=c(UZL_lr$X..contigs, UMC_sr$X..contigs)
w=c(UZL_lr$Total.length, UMC_sr$Total.length)
df=as.data.frame(cbind(x,y,z,w))
df$y = as.numeric(y)/1000000
df$z= as.numeric(z)
df$w = as.numeric(w)/1000000

shapiro.test(y) #very sign
shapiro.test(z) # idem
shapiro.test(w) #idem
y.test=compare_means(y ~ x, data = df,  method = "wilcox.test") #mean: 5807745 vs. 97126
z.test=compare_means(z ~ x, data = df,  method = "wilcox.test") #mean: 4.227273 vs. 309.630952
w.test=compare_means(w ~ x, data = df,  method = "wilcox.test") #mean: 6567680 vs 6496180

one = ggboxplot(df, 'x', 'y', ggtheme = theme_pubr(), xlab = F, ylab = 'N50 (Mbp)', bxp.errorbar = T, add = "mean", add.params = list(color="red"))+
  font("ylab", face="bold", size=14)+
  font("x.text", size=14)+
  stat_pvalue_manual(mutate(y.test, y.position = 8.2), label = "p = {p.adj}")
two = ggboxplot(df, 'x', 'z', ggtheme = theme_pubr(), xlab = F, ylab = 'Number of contigs', bxp.errorbar = T, add = "mean", add.params = list(color="red"))+
  font("ylab", face="bold", size=14)+
  font("x.text", size=14)+
  stat_pvalue_manual(mutate(z.test, y.position = 6050), label = "p = {p.adj}")
three = ggboxplot(df, 'x', 'w', ggtheme = theme_pubr(), xlab = F, ylab = 'Genome size (Mbp)', bxp.errorbar = T, add = "mean", add.params = list(color="red"))+
  font("ylab", face="bold", size=14)+
  font("x.text", size=14)+
  stat_pvalue_manual(mutate(w.test, y.position = 7.8), label = "p = {p.adj}")
grid.arrange(one,two,three,ncol=3)


