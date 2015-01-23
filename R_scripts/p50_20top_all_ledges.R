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
## Shimoni & AlvarezTop 20 Masters regulators

ShA_top20mrs <- summary(ShA_mrshadow_50,mrs=20)$MARINA.results$Regulon
ShA_top20mrs <- as.vector(ShA_top20mrs)

summary(ShA_mrshadow_50,mrs=20)$MARINA.results[,1:5]

## Shadow top 20 pairs

ShA_top20_shadowpairs <- ShA_mrshadow_50$shadow[which(ShA_mrshadow_50$shadow[,1] %in% ShA_top20mrs),]

# Leading-edge of top 20 mrs

ShA_top20_ledges <- ShA_mrshadow_50$ledge[which(names(ShA_mrshadow_50$ledge) %in% ShA_top20mrs)]

## All togueter leading edges

ShA_All_ledges <- unique(unlist(ShA_top20_ledges))


ShA_ledges_exp <- dset[ShA_All_ledges,]

ShA_lg_fit  = lmFit(ShA_ledges_exp, design)
ShA_lg_fit2 = contrasts.fit(ShA_lg_fit, cont.matrix)
ShA_lg_fit2 = eBayes(ShA_lg_fit2)

ShA_lg_topTable <- topTable(ShA_lg_fit2, coef=1, adjust = "fdr", n = 100)

write.table(ShA_lg_topTable,file="/lists/p50_/p50_leading_edges_top20/p50_sha_lg_l-f-Ch_top20.txt",quote=FALSE,sep="\t")

###############################################################################################
###############################################################################################

###############################################################################################
###############################################################################################
## Hugo Top 20 Masters regulators

hugo_top20mrs <- summary(hugo_mrshadow_50,mrs=20)$MARINA.results$Regulon
hugo_top20mrs <- as.vector(hugo_top20mrs)

summary(hugo_mrshadow_50,mrs=20)$MARINA.results[,1:5]

## Shadow top 20 pairs

hugo_top20_shadowpairs <- hugo_mrshadow_50$shadow[which(hugo_mrshadow_50$shadow[,1] %in% hugo_top20mrs),]

# Leading-edge of top 20 mrs

hugo_top20_ledges <- hugo_mrshadow_50$ledge[which(names(hugo_mrshadow_50$ledge) %in% hugo_top20mrs)]

## All togueter leading edges

hugo_All_ledges <- unique(unlist(hugo_top20_ledges))


hugo_ledges_exp <- dset[hugo_All_ledges,]

hugo_lg_fit  = lmFit(hugo_ledges_exp, design)
hugo_lg_fit2 = contrasts.fit(hugo_lg_fit, cont.matrix)
hugo_lg_fit2 = eBayes(hugo_lg_fit2)

hugo_lg_topTable <- topTable(hugo_lg_fit2, coef=1, adjust = "fdr", n = 100)

write.table(hugo_lg_topTable,file="/lists/lp50_/p50_leading_edges_top20/p50_hugo_lg_l-f-Ch_top20.txt",quote=FALSE,sep="\t")

###############################################################################################
###############################################################################################

###############################################################################################
###############################################################################################
## Vaqueriza Top 20 Masters regulators

vaq_top20mrs <- summary(vaq_mrshadow_50,mrs=20)$MARINA.results$Regulon
vaq_top20mrs <- as.vector(vaq_top20mrs)

summary(vaq_mrshadow_50,mrs=20)$MARINA.results[,1:5]

## Shadow top 20 pairs

vaq_top20_shadowpairs <- vaq_mrshadow_50$shadow[which(vaq_mrshadow_50$shadow[,1] %in% vaq_top20mrs),]

# Leading-edge of top 20 mrs

vaq_top20_ledges <- vaq_mrshadow_50$ledge[which(names(vaq_mrshadow_50$ledge) %in% vaq_top20mrs)]

## All togueter leading edges

vaq_All_ledges <- unique(unlist(vaq_top20_ledges))


vaq_ledges_exp <- dset[vaq_All_ledges,]

vaq_lg_fit  = lmFit(vaq_ledges_exp, design)
vaq_lg_fit2 = contrasts.fit(vaq_lg_fit, cont.matrix)
vaq_lg_fit2 = eBayes(vaq_lg_fit2)

vaq_lg_topTable <- topTable(vaq_lg_fit2, coef=1, adjust = "fdr", n = 100)

write.table(vaq_lg_topTable,file="/lists/p50_/p50_leading_edges_top20/p50_vaq_lg_l-f-Ch_top20.txt",quote=FALSE,sep="\t")

###############################################################################################
###############################################################################################

###############################################################################################
###############################################################################################
## Animal TFDB Top 20 Masters regulators

ani_top20mrs <- summary(ani_mrshadow_50,mrs=20)$MARINA.results$Regulon
ani_top20mrs <- as.vector(ani_top20mrs)

summary(ani_mrshadow_50,mrs=20)$MARINA.results[,1:5]

## Shadow top 20 pairs

ani_top20_shadowpairs <- ani_mrshadow_50$shadow[which(ani_mrshadow_50$shadow[,1] %in% ani_top20mrs),]

# Leading-edge of top 20 mrs

ani_top20_ledges <- ani_mrshadow_50$ledge[which(names(ani_mrshadow_50$ledge) %in% ani_top20mrs)]

## All togueter leading edges

ani_All_ledges <- unique(unlist(ani_top20_ledges))


ani_ledges_exp <- dset[ani_All_ledges,]

ani_lg_fit  = lmFit(ani_ledges_exp, design)
ani_lg_fit2 = contrasts.fit(ani_lg_fit, cont.matrix)
ani_lg_fit2 = eBayes(ani_lg_fit2)

ani_lg_topTable <- topTable(ani_lg_fit2, coef=1, adjust = "fdr", n = 100)

write.table(ani_lg_topTable,file="/lists/p50_/p50_leading_edges_top20/p50_ani_lg_l-f-Ch_top20.txt",quote=FALSE,sep="\t")

###############################################################################################
###############################################################################################

## Saving
save.image("results_extractio.RData")





