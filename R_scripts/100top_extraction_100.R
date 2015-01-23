### Estract Results Of MARINA Shadow and Synergy

library(ssmarina)
library(limma)

load("p100_results_MARINA.RData")
load("Rdatas/dset.RData")

design = matrix(rep(0,1760), nrow=880)
colnames(design) = c("enf","san")

design[enfMuest,1]=1
design[sanMuest,2]=1

cont.matrix = makeContrasts("enf - san", levels=design)

toclean <- ls()


###############################################################################################
###############################################################################################
##  100 Masters regulators

ShA_top100mrs_100 <- as.vector(summary(ShA_mrshadow_100,mrs=100)$MARINA.results$Regulon)

hugo_top100mrs_100 <- as.vector(summary(hugo_mrshadow_100,mrs=100)$MARINA.results$Regulon)

vaq_top100mrs_100 <- as.vector(summary(vaq_mrshadow_100,mrs=100)$MARINA.results$Regulon)

ani_top100mrs_100 <- as.vector(summary(ani_mrshadow_100,mrs=100)$MARINA.results$Regulon)

all_top100mrs_100 <- unique(c(ani_top100mrs_100,vaq_top100mrs_100,hugo_top100mrs_100,ShA_top100mrs_100))


exp_all_top100mrs_100 <- dset[all_top100mrs_100,]

all_top100mrs_100fit  = lmFit(exp_all_top100mrs_100, design)
all_top100mrs_100fit2 = contrasts.fit(all_top100mrs_100fit, cont.matrix)
all_top100mrs_100fit2 = eBayes(all_top100mrs_100fit2)

all_top100mrs_100topTable <- topTable(all_top100mrs_100fit2, coef=1, adjust = "fdr", n = length(all_top100mrs_100))

write.table(all_top100mrs_100topTable,file="./lists/p100/all_top100mrs_100_LfCh.txt",quote=FALSE,sep="\t")
###############################################################################################

exp_ShA_top100mrs_100 <- dset[ShA_top100mrs_100,]

ShA_top100mrs_100fit  = lmFit(exp_ShA_top100mrs_100, design)
ShA_top100mrs_100fit2 = contrasts.fit(ShA_top100mrs_100fit, cont.matrix)
ShA_top100mrs_100fit2 = eBayes(ShA_top100mrs_100fit2)

ShA_top100mrs_100topTable <- topTable(ShA_top100mrs_100fit2, coef=1, adjust = "fdr", n = length(ShA_top100mrs_100))

write.table(ShA_top100mrs_100topTable,file="./lists/p100/ShA_top100mrs_100_LfCh.txt",quote=FALSE,sep="\t")
###############################################################################################
exp_hugo_top100mrs_100 <- dset[hugo_top100mrs_100,]

hugo_top100mrs_100fit  = lmFit(exp_hugo_top100mrs_100, design)
hugo_top100mrs_100fit2 = contrasts.fit(hugo_top100mrs_100fit, cont.matrix)
hugo_top100mrs_100fit2 = eBayes(hugo_top100mrs_100fit2)

hugo_top100mrs_100topTable <- topTable(hugo_top100mrs_100fit2, coef=1, adjust = "fdr", n = length(hugo_top100mrs_100))

write.table(hugo_top100mrs_100topTable,file="./lists/p100/hugo_top100mrs_100_LfCh.txt",quote=FALSE,sep="\t")
###############################################################################################
exp_ani_top100mrs_100 <- dset[ani_top100mrs_100,]

ani_top100mrs_100fit  = lmFit(exp_ani_top100mrs_100, design)
ani_top100mrs_100fit2 = contrasts.fit(ani_top100mrs_100fit, cont.matrix)
ani_top100mrs_100fit2 = eBayes(ani_top100mrs_100fit2)

ani_top100mrs_100topTable <- topTable(ani_top100mrs_100fit2, coef=1, adjust = "fdr", n = length(ani_top100mrs_100))

write.table(ani_top100mrs_100topTable,file="./lists/p100/ani_top100mrs_100_LfCh.txt",quote=FALSE,sep="\t")
###############################################################################################
exp_vaq_top100mrs_100 <- dset[vaq_top100mrs_100,]

vaq_top100mrs_100fit  = lmFit(exp_vaq_top100mrs_100, design)
vaq_top100mrs_100fit2 = contrasts.fit(vaq_top100mrs_100fit, cont.matrix)
vaq_top100mrs_100fit2 = eBayes(vaq_top100mrs_100fit2)

vaq_top100mrs_100topTable <- topTable(vaq_top100mrs_100fit2, coef=1, adjust = "fdr", n = length(vaq_top100mrs_100))

write.table(vaq_top100mrs_100topTable,file="./lists/p100/vaq_top100mrs_100_LfCh.txt",quote=FALSE,sep="\t")
###############################################################################################

write.table(ShA_top100mrs_100,file="./lists/p100/p100_ShA_top100mrs.txt",quote=FALSE,sep="\t")
write.table(vaq_top100mrs_100,file="./lists/p100/p100_vaq_top100mrs.txt",quote=FALSE,sep="\t")
write.table(hugo_top100mrs_100,file="./lists/p100/p100_hugo_top100mrs.txt",quote=FALSE,sep="\t")
write.table(ani_top100mrs_100,file="./lists/p100/p100_ani_top100mrs.txt",quote=FALSE,sep="\t")
write.table(all_top100mrs_100,file="./lists/p100/p100_all_top100.txt",quote=FALSE,sep="\t")






summary(ShA_mrshadow_100,mrs=100)$MARINA.results[,1:5]
summary(hugo_mrshadow_100,mrs=100)$MARINA.results[,1:5]
summary(vaq_mrshadow_100,mrs=100)$MARINA.results[,1:5]
summary(ani_mrshadow_100,mrs=100)$MARINA.results[,1:5]

