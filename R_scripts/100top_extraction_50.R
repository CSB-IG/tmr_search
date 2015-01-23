### Estract Results Of MARINA Shadow and Synergy

library(ssmarina)
library(limma)

load("p50_results_MARINA.RData")
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

ShA_top100mrs_50 <- as.vector(summary(ShA_mrshadow_50,mrs=100)$MARINA.results$Regulon)

hugo_top100mrs_50 <- as.vector(summary(hugo_mrshadow_50,mrs=100)$MARINA.results$Regulon)

vaq_top100mrs_50 <- as.vector(summary(vaq_mrshadow_50,mrs=100)$MARINA.results$Regulon)

ani_top100mrs_50 <- as.vector(summary(ani_mrshadow_50,mrs=100)$MARINA.results$Regulon)

all_top100mrs_50 <- unique(c(ani_top100mrs_50,vaq_top100mrs_50,hugo_top100mrs_50,ShA_top100mrs_50))


exp_all_top100mrs_50 <- dset[all_top100mrs_50,]

all_top100mrs_50fit  = lmFit(exp_all_top100mrs_50, design)
all_top100mrs_50fit2 = contrasts.fit(all_top100mrs_50fit, cont.matrix)
all_top100mrs_50fit2 = eBayes(all_top100mrs_50fit2)

all_top100mrs_50topTable <- topTable(all_top100mrs_50fit2, coef=1, adjust = "fdr", n = length(all_top100mrs_50))

write.table(all_top100mrs_50topTable,file="./lists/p50/all_top100mrs_50_LfCh.txt",quote=FALSE,sep="\t")
###############################################################################################

exp_ShA_top100mrs_50 <- dset[ShA_top100mrs_50,]

ShA_top100mrs_50fit  = lmFit(exp_ShA_top100mrs_50, design)
ShA_top100mrs_50fit2 = contrasts.fit(ShA_top100mrs_50fit, cont.matrix)
ShA_top100mrs_50fit2 = eBayes(ShA_top100mrs_50fit2)

ShA_top100mrs_50topTable <- topTable(ShA_top100mrs_50fit2, coef=1, adjust = "fdr", n = length(ShA_top100mrs_50))

write.table(ShA_top100mrs_50topTable,file="./lists/p50/ShA_top100mrs_50_LfCh.txt",quote=FALSE,sep="\t")
###############################################################################################
exp_hugo_top100mrs_50 <- dset[hugo_top100mrs_50,]

hugo_top100mrs_50fit  = lmFit(exp_hugo_top100mrs_50, design)
hugo_top100mrs_50fit2 = contrasts.fit(hugo_top100mrs_50fit, cont.matrix)
hugo_top100mrs_50fit2 = eBayes(hugo_top100mrs_50fit2)

hugo_top100mrs_50topTable <- topTable(hugo_top100mrs_50fit2, coef=1, adjust = "fdr", n = length(hugo_top100mrs_50))

write.table(hugo_top100mrs_50topTable,file="./lists/p50/hugo_top100mrs_50_LfCh.txt",quote=FALSE,sep="\t")
###############################################################################################
exp_ani_top100mrs_50 <- dset[ani_top100mrs_50,]

ani_top100mrs_50fit  = lmFit(exp_ani_top100mrs_50, design)
ani_top100mrs_50fit2 = contrasts.fit(ani_top100mrs_50fit, cont.matrix)
ani_top100mrs_50fit2 = eBayes(ani_top100mrs_50fit2)

ani_top100mrs_50topTable <- topTable(ani_top100mrs_50fit2, coef=1, adjust = "fdr", n = length(ani_top100mrs_50))

write.table(ani_top100mrs_50topTable,file="./lists/p50/ani_top100mrs_50_LfCh.txt",quote=FALSE,sep="\t")
###############################################################################################
exp_vaq_top100mrs_50 <- dset[vaq_top100mrs_50,]

vaq_top100mrs_50fit  = lmFit(exp_vaq_top100mrs_50, design)
vaq_top100mrs_50fit2 = contrasts.fit(vaq_top100mrs_50fit, cont.matrix)
vaq_top100mrs_50fit2 = eBayes(vaq_top100mrs_50fit2)

vaq_top100mrs_50topTable <- topTable(vaq_top100mrs_50fit2, coef=1, adjust = "fdr", n = length(vaq_top100mrs_50))

write.table(vaq_top100mrs_50topTable,file="./lists/p50/vaq_top100mrs_50_LfCh.txt",quote=FALSE,sep="\t")
###############################################################################################

write.table(ShA_top100mrs_50,file="./lists/p50/p50_ShA_top100mrs.txt",quote=FALSE,sep="\t",row.names=F,col.names=F)
write.table(vaq_top100mrs_50,file="./lists/p50/p50_vaq_top100mrs.txt",quote=FALSE,sep="\t",row.names=F,col.names=F)
write.table(hugo_top100mrs_50,file="./lists/p50/p50_hugo_top100mrs.txt",quote=FALSE,sep="\t",row.names=F,col.names=F)
write.table(ani_top100mrs_50,file="./lists/p50/p50_ani_top100mrs.txt",quote=FALSE,sep="\t",row.names=F,col.names=F)
write.table(all_top100mrs_50,file="./lists/p50/p50_all_top100.txt",quote=FALSE,sep="\t",row.names=F,col.names=F)









summary(ShA_mrshadow_50,mrs=100)$MARINA.results[,1:5]
summary(hugo_mrshadow_50,mrs=100)$MARINA.results[,1:5]
summary(vaq_mrshadow_50,mrs=100)$MARINA.results[,1:5]
summary(ani_mrshadow_50,mrs=100)$MARINA.results[,1:5]














## Shadow top 20 pairs

ShA_top20_shadowpairs <- ShA_mrshadow_50$shadow[which(ShA_mrshadow_50$shadow[,1] %in% ShA_top20mrs),]

# Leading-edge of top 20 mrs

ShA_top20_ledges <- ShA_mrshadow_50$ledge[which(names(ShA_mrshadow_50$ledge) %in% ShA_top20mrs)]

## Extract leading-edges from every one gene

# for (i in 1:length(ShA_top20_ledges)){

#   n <-  ShA_top20mrs[i] 
#   v <- matrix(unlist(ShA_top20_ledges[n]))
#   colnames(v) <- n
#   write.table(v, file=paste(n, "ShA_top20_ledges.txt",sep="_"), row.names=F, quote=F)
# }

## All togueter leading edges

ShA_All_ledges <- unique(unlist(ShA_top20_ledges))



###############################################################################################
###############################################################################################

###############################################################################################
###############################################################################################
## Hugo Top 20 Masters regulators




## Shadow top 20 pairs

hugo_top20_shadowpairs <- hugo_mrshadow_50$shadow[which(hugo_mrshadow_50$shadow[,1] %in% hugo_top20mrs),]

# Leading-edge of top 20 mrs

hugo_top20_ledges <- hugo_mrshadow_50$ledge[which(names(hugo_mrshadow_50$ledge) %in% hugo_top20mrs)]

## Extract leading-edges from every one gene

# for (i in 1:length(hugo_top20_ledges)){

#   n <-  hugo_top20mrs[i] 
#   v <- matrix(unlist(hugo_top20_ledges[n]))
#   colnames(v) <- n
#   write.table(v, file=paste(n, "hugo_top20_ledges.txt",sep="_"), row.names=F, quote=F)
# }

## All togueter leading edges

hugo_All_ledges <- unique(unlist(hugo_top20_ledges))


hugo_ledges_exp <- dset[hugo_All_ledges,]

hugo_lg_fit  = lmFit(hugo_ledges_exp, design)
hugo_lg_fit2 = contrasts.fit(hugo_lg_fit, cont.matrix)
hugo_lg_fit2 = eBayes(hugo_lg_fit2)

hugo_lg_topTable <- topTable(hugo_lg_fit2, coef=1, adjust = "fdr", n = 100)

write.table(hugo_lg_topTable,file="/lists/leading_edges_top20/hugo_lg_l-f-Ch_top20.txt",quote=FALSE,sep="\t")

###############################################################################################
###############################################################################################

###############################################################################################
###############################################################################################
## Vaqueriza Top 20 Masters regulators




## Shadow top 20 pairs

vaq_top20_shadowpairs <- vaq_mrshadow_50$shadow[which(vaq_mrshadow_50$shadow[,1] %in% vaq_top20mrs),]

# Leading-edge of top 20 mrs

vaq_top20_ledges <- vaq_mrshadow_50$ledge[which(names(vaq_mrshadow_50$ledge) %in% vaq_top20mrs)]

## Extract leading-edges from every one gene

# for (i in 1:length(vaq_top20_ledges)){

#   n <-  vaq_top20mrs[i] 
#   v <- matrix(unlist(vaq_top20_ledges[n]))
#   colnames(v) <- n
#   write.table(v, file=paste(n, "vaq_top20_ledges.txt",sep="_"), row.names=F, quote=F)
# }

## All togueter leading edges

vaq_All_ledges <- unique(unlist(vaq_top20_ledges))


vaq_ledges_exp <- dset[vaq_All_ledges,]

vaq_lg_fit  = lmFit(vaq_ledges_exp, design)
vaq_lg_fit2 = contrasts.fit(vaq_lg_fit, cont.matrix)
vaq_lg_fit2 = eBayes(vaq_lg_fit2)

vaq_lg_topTable <- topTable(vaq_lg_fit2, coef=1, adjust = "fdr", n = 100)

write.table(vaq_lg_topTable,file="/lists/leading_edges_top20/vaq_lg_l-f-Ch_top20.txt",quote=FALSE,sep="\t")

###############################################################################################
###############################################################################################

###############################################################################################
###############################################################################################
## Animal TFDB Top 20 Masters regulators



summary(ani_mrshadow_50,mrs=20)$MARINA.results[,1:5]

## Shadow top 20 pairs

ani_top20_shadowpairs <- ani_mrshadow_50$shadow[which(ani_mrshadow_50$shadow[,1] %in% ani_top20mrs),]

# Leading-edge of top 20 mrs

ani_top20_ledges <- ani_mrshadow_50$ledge[which(names(ani_mrshadow_50$ledge) %in% ani_top20mrs)]

## Extract leading-edges from every one gene

# for (i in 1:length(ani_top20_ledges)){

#   n <-  ani_top20mrs[i] 
#   v <- matrix(unlist(ani_top20_ledges[n]))
#   colnames(v) <- n
#   write.table(v, file=paste(n, "ani_top20_ledges.txt",sep="_"), row.names=F, quote=F)
# }

## All togueter leading edges

ani_All_ledges <- unique(unlist(ani_top20_ledges))


ani_ledges_exp <- dset[ani_All_ledges,]

ani_lg_fit  = lmFit(ani_ledges_exp, design)
ani_lg_fit2 = contrasts.fit(ani_lg_fit, cont.matrix)
ani_lg_fit2 = eBayes(ani_lg_fit2)

ani_lg_topTable <- topTable(ani_lg_fit2, coef=1, adjust = "fdr", n = 100)

write.table(ani_lg_topTable,file="/lists/leading_edges_top20/ani_lg_l-f-Ch_top20.txt",quote=FALSE,sep="\t")

###############################################################################################
###############################################################################################

## Saving
save.image("results_extractio.RData")





