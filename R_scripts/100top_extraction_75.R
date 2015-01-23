### Estract Results Of MARINA Shadow and Synergy

library(ssmarina)
library(limma)

load("p75_results_MARINA.RData")
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

ShA_top100mrs_75 <- as.vector(summary(ShA_mrshadow_75,mrs=100)$MARINA.results$Regulon)

hugo_top100mrs_75 <- as.vector(summary(hugo_mrshadow_75,mrs=100)$MARINA.results$Regulon)

vaq_top100mrs_75 <- as.vector(summary(vaq_mrshadow_75,mrs=100)$MARINA.results$Regulon)

ani_top100mrs_75 <- as.vector(summary(ani_mrshadow_75,mrs=100)$MARINA.results$Regulon)

all_top100mrs_75 <- unique(c(ani_top100mrs_75,vaq_top100mrs_75,hugo_top100mrs_75,ShA_top100mrs_75))


exp_all_top100mrs_75 <- dset[all_top100mrs_75,]

all_top100mrs_75fit  = lmFit(exp_all_top100mrs_75, design)
all_top100mrs_75fit2 = contrasts.fit(all_top100mrs_75fit, cont.matrix)
all_top100mrs_75fit2 = eBayes(all_top100mrs_75fit2)

all_top100mrs_75topTable <- topTable(all_top100mrs_75fit2, coef=1, adjust = "fdr", n = length(all_top100mrs_75))

write.table(all_top100mrs_75topTable,file="./lists/p75/all_top100mrs_75_LfCh.txt",quote=FALSE,sep="\t")
###############################################################################################

exp_ShA_top100mrs_75 <- dset[ShA_top100mrs_75,]

ShA_top100mrs_75fit  = lmFit(exp_ShA_top100mrs_75, design)
ShA_top100mrs_75fit2 = contrasts.fit(ShA_top100mrs_75fit, cont.matrix)
ShA_top100mrs_75fit2 = eBayes(ShA_top100mrs_75fit2)

ShA_top100mrs_75topTable <- topTable(ShA_top100mrs_75fit2, coef=1, adjust = "fdr", n = length(ShA_top100mrs_75))

write.table(ShA_top100mrs_75topTable,file="./lists/p75/ShA_top100mrs_75_LfCh.txt",quote=FALSE,sep="\t")
###############################################################################################
exp_hugo_top100mrs_75 <- dset[hugo_top100mrs_75,]

hugo_top100mrs_75fit  = lmFit(exp_hugo_top100mrs_75, design)
hugo_top100mrs_75fit2 = contrasts.fit(hugo_top100mrs_75fit, cont.matrix)
hugo_top100mrs_75fit2 = eBayes(hugo_top100mrs_75fit2)

hugo_top100mrs_75topTable <- topTable(hugo_top100mrs_75fit2, coef=1, adjust = "fdr", n = length(hugo_top100mrs_75))

write.table(hugo_top100mrs_75topTable,file="./lists/p75/hugo_top100mrs_75_LfCh.txt",quote=FALSE,sep="\t")
###############################################################################################
exp_ani_top100mrs_75 <- dset[ani_top100mrs_75,]

ani_top100mrs_75fit  = lmFit(exp_ani_top100mrs_75, design)
ani_top100mrs_75fit2 = contrasts.fit(ani_top100mrs_75fit, cont.matrix)
ani_top100mrs_75fit2 = eBayes(ani_top100mrs_75fit2)

ani_top100mrs_75topTable <- topTable(ani_top100mrs_75fit2, coef=1, adjust = "fdr", n = length(ani_top100mrs_75))

write.table(ani_top100mrs_75topTable,file="./lists/p75/ani_top100mrs_75_LfCh.txt",quote=FALSE,sep="\t")
###############################################################################################
exp_vaq_top100mrs_75 <- dset[vaq_top100mrs_75,]

vaq_top100mrs_75fit  = lmFit(exp_vaq_top100mrs_75, design)
vaq_top100mrs_75fit2 = contrasts.fit(vaq_top100mrs_75fit, cont.matrix)
vaq_top100mrs_75fit2 = eBayes(vaq_top100mrs_75fit2)

vaq_top100mrs_75topTable <- topTable(vaq_top100mrs_75fit2, coef=1, adjust = "fdr", n = length(vaq_top100mrs_75))

write.table(vaq_top100mrs_75topTable,file="./lists/p75/vaq_top100mrs_75_LfCh.txt",quote=FALSE,sep="\t")
###############################################################################################

write.table(ShA_top100mrs_75,file="./lists/p75/p75_ShA_top100mrs.txt",quote=FALSE,sep="\t",row.names=F,col.names=F)
write.table(vaq_top100mrs_75,file="./lists/p75/p75_vaq_top100mrs.txt",quote=FALSE,sep="\t",row.names=F,col.names=F)
write.table(hugo_top100mrs_75,file="./lists/p75/p75_hugo_top100mrs.txt",quote=FALSE,sep="\t",row.names=F,col.names=F)
write.table(ani_top100mrs_75,file="./lists/p75/p75_ani_top100mrs.txt",quote=FALSE,sep="\t",row.names=F,col.names=F)
write.table(all_top100mrs_75,file="./lists/p75/p75_all_top100.txt",quote=FALSE,sep="\t",row.names=F,col.names=F)









summary(ShA_mrshadow_75,mrs=100)$MARINA.results[,1:5]
summary(hugo_mrshadow_75,mrs=100)$MARINA.results[,1:5]
summary(vaq_mrshadow_75,mrs=100)$MARINA.results[,1:5]
summary(ani_mrshadow_75,mrs=100)$MARINA.results[,1:5]