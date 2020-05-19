######################################################################
#### Consolidation of infectivity data and corresponding analyses ####
######################################################################

##### 0 - Load packages and data ####
library(readxl)
library(ggfortify)
library(corrplot)
library(caret)
library(tree)

#get workspace from Acr.R ('all'/'all_woIV' dataframes in particular)

##load UMC infectivity data 
ID_conversion <- read_excel("C:/Users/emmav/OneDrive/Documenten/Thesis/4Data/Phage_infectivity/UMC_ordernr-seqfileID_conversion.xls")
ID_conversion = ID_conversion[-c(496),] #remove contaminated strain
PA_IMI_8_phages <- read_excel("~/Thesis/4Data/Phage infectivity/UMC_original.xlsx", sheet = "master_sheet (2)", range = "A1:I453")
names(PA_IMI_8_phages)[1]='ID'

x = merge(PA_IMI_8_phages, ID_conversion, by.x = 'ID', by.y = 'phageresistnr')
notsequenced=which(is.na(x$seqnr_corrected))

UMC_infect = x[-notsequenced, c(11, 2:9)]
names(UMC_infect)[c(1,2,4,6:8)]=c('Strain','ph141','LKD16','LUZ7','LUZ19','LUZ24') #ph141 = 14-1
UMC_inf = UMC_infect

#convert PFU/ml to 0-1 encoding
for (i in 2:9){
  UMC_inf[,i] = ifelse(UMC_inf[,i]>0, 1, 0)
}
#2 of the outliers have to be removed:
UMC_inf = UMC_inf[-c(149,300),]

##load LoGt infectivity data
LoGT_inf <- read_excel("~/Thesis/4Data/Phage infectivity/HostRange_LoGT.xlsx", range = "A1:X72")
names(LoGT_inf)[c(1,4,15)] = c('Strain', 'phiKMV', 'ph141')

##load SGI infectivity data
SGI_inf <- read_excel("~/Thesis/4Data/Phage infectivity/Infectivity_all.xlsx", sheet = "SGI")
names(SGI_inf)[c(1,2)] = c('Strain','ph141')
#3 strains not in full dataset: 
SGI_inf = SGI_inf[-c(79,85,253),]

##combine all 3 sets
x = merge(Saar_inf, SGI_inf, all=T)
x = merge(x, UMC_inf, all=T)
all_complete = merge(x, all, by='Strain')

all_complete_woIV = all_complete
all_complete_woIV = all_complete_woIV[-c(which(all_complete_woIV$U > 0)),] #8 removed
attach(all_complete_woIV)



#### 1 - Exploratory analysis: PCA of LoGT set alone ####

#without IV: 66 strains in total +3 incomplete cases, 23 phages 
#add other variables for visualizations
x_ext = merge(Saar_inf, all_woIV, by = c('Strain'))
x_ext = x_ext[complete.cases(x_ext[,c(2:24)]),]

x.pca = princomp(x_ext[,c(2:24)]) #centered but not scaled, on covariance matrix
x.pca
summary(x.pca) #first two: 46 % of variance
screeplot(x.pca, type = 'lines')

##visualization
autoplot(x.pca, data = x_ext, colour = 'Group', shape = 'Group', size = 3, frame = TRUE)+theme_bw()
autoplot(x.pca, data = x_ext, colour = 'Size', size = 3)
x_ext$SizeGroup = ifelse(x_ext$Size > 6.7, 'big', 'small')
autoplot(x.pca, data = x_ext, colour = 'SizeGroup',shape='SizeGroup', size = 3)

summary(as.factor(x_ext$ctype)) #33 I-F, 22 None and 11 of others (small set of strains)

x=data.frame(cbind(x.pca$scores, x_ext))
x$f = factor(x$ctype == 'IF', labels = c('Other types','Type I-F'))
ggplot(data = x, aes(x = Comp.1, y = Comp.2)) + 
  geom_point(aes(color = ctype_corc_adj_order, fill = ctype_corc_adj_order, shape=ctype_corc_adj_order, size = Size))+
  facet_grid(rows= vars(x$f) )+
  theme_bw()+scale_shape_manual(values = c(22, 21, 23, 10, 25, 8, 12))+
  ylab("Principal componenent 2 (10.85%)") + xlab("Principal component 1 (35.34%)")+
  theme(legend.position = 'top', legend.title = element_blank(), strip.text.y = element_text(size = 14), legend.text = element_text(size = 12))


#### 2 - Correlation analysis ####

#use complete dataset: 5 phages tested on 736 strains
all_complete_woIV$TIF = ifelse(all_complete_woIV$ctype=="IF", TRUE, FALSE)

#function to compute significance matrix
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ..., method = 'spearman')
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

#LUZ19

luz19 = all_complete_woIV[complete.cases(LUZ19),c("LUZ19","Size","Proteins","cas","spacers","cc_adj2","corc_adj","ctype_adj2","ctype_corc_adj","TIF","nb_hits_hom","nb_other_hom","GBA_woLC","ctype_corc")] 
luz19$cc_adj2 = as.numeric(luz19$cc_adj2); luz19$corc_adj = as.numeric(luz19$corc_adj)
p.mat = cor.mtest(luz19[,-c(3,8,9,10,14)])
cmat = cor(luz19[,-c(3,8,9,10,14)], method = 'spearman')
colnames(cmat) = c('LUZ19', 'Genome size','Cas','Spacers','Corrected CC','Corrected CorC','Homology-Acr','Homology-Aca','GBA')
rownames(cmat) = c('LUZ19', 'Genome size','Cas','Spacers','Corrected CC','Corrected CorC','Homology-Acr','Homology-Aca','GBA')
corrplot(cmat,p.mat = p.mat, addCoef.col = "black",tl.cex = 1.7, pch.cex = 7, number.cex = 1.5, number.font = 1, cl.cex = 1.5,
         method = 'color', type = 'lower', tl.col = "black",tl.srt = 45, diag = T, outline = T) 

#14-1, LBL3, LUZ7 and LKD16

ph141 = all_complete_woIV[complete.cases(ph141),c("ph141","Size","Proteins","cas","spacers","cc_adj2","corc_adj","ctype_adj2","ctype_corc_adj","TIF","nb_hits_hom","nb_other_hom","GBA_woLC")] 
ph141$cc_adj2 = as.numeric(ph141$cc_adj2); ph141$corc_adj = as.numeric(ph141$corc_adj)
c1 = cor(ph141[,-c(3,8,9,10)], method = 'spearman')[,1]
p1 = cor.mtest(ph141[,-c(3,8,9,10)])[,1]

lbl3 = all_complete_woIV[complete.cases(LBL3),c("LBL3","Size","Proteins","cas","spacers","cc_adj2","corc_adj","ctype_adj2","ctype_corc_adj","TIF","nb_hits_hom","nb_other_hom","GBA_woLC")] 
lbl3$cc_adj2 = as.numeric(lbl3$cc_adj2); lbl3$corc_adj = as.numeric(lbl3$corc_adj)
c2 = cor(lbl3[,-c(3,8,9,10)], method = 'spearman')[,1]
p2 = cor.mtest(lbl3[,-c(3,8,9,10)])[,1]

luz7 = all_complete_woIV[complete.cases(LUZ7),c("LUZ7","Size","Proteins","cas","spacers","cc_adj2","corc_adj","ctype_adj2","ctype_corc_adj","TIF","nb_hits_hom","nb_other_hom","GBA_woLC")] 
luz7$cc_adj2 = as.numeric(luz7$cc_adj2); luz7$corc_adj = as.numeric(luz7$corc_adj)
c3 = cor(luz7[,-c(3,8,9,10)], method = 'spearman')[,1] 
p3 = cor.mtest(luz7[,-c(3,8,9,10)])[,1]

lkd16 = all_complete_woIV[complete.cases(LKD16),c("LKD16","Size","Proteins","cas","spacers","cc_adj2","corc_adj","ctype_adj2","ctype_corc_adj","TIF","nb_hits_hom","nb_other_hom","GBA_woLC")] 
lkd16$cc_adj2 = as.numeric(lkd16$cc_adj2); lkd16$corc_adj = as.numeric(lkd16$corc_adj)
c4 = cor(lkd16[,-c(3,8,9,10)], method = 'spearman')[,1]
p4 = cor.mtest(lkd16[,-c(3,8,9,10)])[,1]

c=rbind(c1,c2,c3,c4)[,-1]
p = rbind(p1,p2,p3,p4)[,-1]
colnames(c) = c('Genome size','Cas','Spacers','Corrected CC','Corrected CorC','Homology-Acr','Homology-Aca','GBA')
rownames(c) = c('14-1','LBL3','LUZ7','LKD16')
corrplot(t(c),p.mat = t(p), addCoef.col = "black", tl.cex = 1.7, pch.cex = 8, number.cex = 1.5, number.font = 1, cl.cex = 1.5,
         method = 'color', tl.col = "black",tl.srt = 45, outline = T) 




#### 3 - Modelling the infectivity of LUZ19 ####

luz19$GBA_woLC = as.factor(luz19$GBA_woLC); luz19$cc_adj2 = as.factor(luz19$cc_adj2); luz19$corc_adj = as.factor(luz19$corc_adj); luz19$ctype_adj2 = as.factor(luz19$ctype_adj2); 
luz19$ctype_corc_adj = as.factor(luz19$ctype_corc_adj);luz19$ctype_corc = as.factor(luz19$ctype_corc); luz19$LUZ19= as.factor(luz19$LUZ19)

#make train and test set
set.seed(100)
trainset = createDataPartition(luz19$LUZ19, p=0.9, list = F)
trainData = luz19[trainset, ] 
testData = luz19[-trainset, ]
down_train = downSample(x = trainData[,-c(1)], y = as.factor(trainData$LUZ19), yname='LUZ19')

#exploratory, on full dataset
#ctype_corc is included in luz19 dataframe to assess type I-F effect but should not be included in the models
summary(tree(LUZ19~.-ctype_corc, luz19)) 
summary(tree(LUZ19~.-Proteins-ctype_corc, luz19))
summary(tree(LUZ19~.-Size-Proteins-ctype_corc, luz19)) 
summary(tree(LUZ19~.-Size-Proteins-ctype_corc-ctype_corc_adj-ctype_adj2-cc_adj2-corc_adj, luz19))

#now use training and test sets
tree1 = tree(LUZ19~.-Proteins-ctype_corc, down_train)
plot(tree1);text(tree1, pretty = 0) 
cv.tree1 = cv.tree(tree1, FUN = prune.misclass)
plot(cv.tree1) 
prune.tree1 = prune.misclass(tree1, best = 5)
plot(prune.tree1);text(prune.tree1, pretty=0)
tree.pred = predict(prune.tree1, testData, type="class")
table(tree.pred, testData$LUZ19)  


#### 4 - LUZ19 infectivity rates per type ####

summary(luz19$ctype_corc_adj)
f=summary(luz19$LUZ19[which(luz19$ctype_corc_adj == 'IF')])
e=summary(luz19$LUZ19[which(luz19$ctype_corc_adj == 'IE')])
c=summary(luz19$LUZ19[which(luz19$ctype_corc_adj == 'IC')]) # too few observations
co=summary(luz19$LUZ19[which(luz19$ctype_corc_adj == 'combination')])
no=summary(luz19$LUZ19[which(luz19$ctype_corc_adj == 'None')])
oc=summary(luz19$LUZ19[which(luz19$ctype_corc_adj == 'Only CRISPR')])
i=summary(luz19$LUZ19[which(luz19$ctype_corc_adj == 'Inactivated')])

chisq.test(rbind(no, i)) #p = 0.1431
chisq.test(rbind(f, e)) # p = 0.69

chisq.test(rbind(e+f+oc, no+i)) # p = 4.684e-06

#type I-F effect?
f2=summary(luz19$LUZ19[which(luz19$ctype_corc == 'IF' & luz19$ctype_corc_adj=="Inactivated")])
i2=summary(luz19$LUZ19[which(luz19$ctype_corc_adj == 'Inactivated' & luz19$ctype_corc != 'IF')])
chisq.test(rbind(f2, i2)) #not different p = 0.3345
