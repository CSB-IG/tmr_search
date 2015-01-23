### Estract Results Of MARINA Shadow and Synergy

library(ssmarina)
library(limma)

load("p10_results_MARINA.RData")
load("Rdatas/dset.RData")

design = matrix(rep(0,1760), nrow=880)
colnames(design) = c("enf","san")

design[enfMuest,1]=1
design[sanMuest,2]=1

cont.matrix = makeContrasts("enf - san", levels=design)

toclean <- ls()

###############################################################################################
###############################################################################################

# Top 100 masters regulators

ani_top100mrs <- summary(ani_mrshadow_10,mrs=100)$Regulon
ani_top100mrs <- as.vector(ani_top100mrs)

hugo_top100mrs <- summary(hugo_mrshadow_10,mrs=100)$Regulon
hugo_top100mrs <- as.vector(hugo_top100mrs)

ShA_top100mrs <- summary(ShA_mrshadow_10,mrs=100)$Regulon
ShA_top100mrs <- as.vector(ShA_top100mrs)

vaq_top100mrs <- summary(vaq_mrshadow_10,mrs=100)$Regulon
vaq_top100mrs <- as.vector(vaq_top100mrs)

all_top100mrs <- unique(c(ani_top100mrs,hugo_top100mrs,ShA_top100mrs,vaq_top100mrs))

###############################################################################################
###############################################################################################

# Leading-edge of top 100 mrs

# ani_top100mrs_ledges <- ani_mrshadow_10$ledge[which(names(ani_mrshadow_10$ledge) %in% ani_top100mrs)]
# hugo_top100mrs_ledges <- hugo_mrshadow_10$ledge[which(names(hugo_mrshadow_10$ledge) %in% hugo_top100mrs)]
# ShA_top100mrs_ledges <- ShA_mrshadow_10$ledge[which(names(ShA_mrshadow_10$ledge) %in% ShA_top100mrs)]
# vaq_top100mrs_ledges <- vaq_mrshadow_10$ledge[which(names(vaq_mrshadow_10$ledge) %in% vaq_top100mrs)]

ani_All_ledges <- unique(unlist(ani_mrshadow_10$ledge[ani_top100mrs]))
hugo_All_ledges <- unique(unlist(hugo_mrshadow_10$ledge[hugo_top100mrs]))
ShA_All_ledges <- unique(unlist(ShA_mrshadow_10$ledge[ShA_top100mrs]))
vaq_All_ledges <- unique(unlist(vaq_mrshadow_10$ledge[vaq_top100mrs]))

all_All_ledges <- unique(c(ani_All_ledges, hugo_All_ledges, ShA_All_ledges, vaq_All_ledges))

###############################################################################################

# Log fold-change

all_ledges_exp <- dset[all_All_ledges,]

all_lg_fit  = lmFit(all_ledges_exp, design)
all_lg_fit2 = contrasts.fit(all_lg_fit, cont.matrix)
all_lg_fit2 = eBayes(all_lg_fit2)

all_lg_topTable <- topTable(all_lg_fit2, coef=1, adjust = "fdr", n = length(all_All_ledges))

write.table(all_lg_topTable,file="./lists/p10/p10_top100mrs_ledges/p10_all_top100mrs_ledges_LFCh.txt",quote=FALSE,sep="\t")

###############################################################################################

# Log fold-change

ani_ledges_exp <- dset[ani_All_ledges,]

ani_lg_fit  = lmFit(ani_ledges_exp, design)
ani_lg_fit2 = contrasts.fit(ani_lg_fit, cont.matrix)
ani_lg_fit2 = eBayes(ani_lg_fit2)

ani_lg_topTable <- topTable(ani_lg_fit2, coef=1, adjust = "fdr", n = length(ani_All_ledges))

write.table(ani_lg_topTable,file="./lists/p10/p10_top100mrs_ledges/p10_ani_top100mrs_ledges_LFCh.txt",quote=FALSE,sep="\t")

###############################################################################################

# Log fold-change

hugo_ledges_exp <- dset[hugo_All_ledges,]

hugo_lg_fit  = lmFit(hugo_ledges_exp, design)
hugo_lg_fit2 = contrasts.fit(hugo_lg_fit, cont.matrix)
hugo_lg_fit2 = eBayes(hugo_lg_fit2)

hugo_lg_topTable <- topTable(hugo_lg_fit2, coef=1, adjust = "fdr", n = 100)

write.table(hugo_lg_topTable,file="./lists/p10/p10_top100mrs_ledges/p10_hugo_top100mrs_ledges_LFCh.txt",quote=FALSE,sep="\t")

###############################################################################################

# Log fold-change

ShA_ledges_exp <- dset[ShA_All_ledges,]

ShA_lg_fit  = lmFit(ShA_ledges_exp, design)
ShA_lg_fit2 = contrasts.fit(ShA_lg_fit, cont.matrix)
ShA_lg_fit2 = eBayes(ShA_lg_fit2)

ShA_lg_topTable <- topTable(ShA_lg_fit2, coef=1, adjust = "fdr", n = 100)

write.table(ShA_lg_topTable,file="./lists/p10/p10_top100mrs_ledges/p10_sha_top100mrs_ledges_LFCh.txt",quote=FALSE,sep="\t")

###############################################################################################

# Log fold-change

vaq_ledges_exp <- dset[vaq_All_ledges,]

vaq_lg_fit  = lmFit(vaq_ledges_exp, design)
vaq_lg_fit2 = contrasts.fit(vaq_lg_fit, cont.matrix)
vaq_lg_fit2 = eBayes(vaq_lg_fit2)

vaq_lg_topTable <- topTable(vaq_lg_fit2, coef=1, adjust = "fdr", n = 100)

write.table(vaq_lg_topTable,file="./lists/p10/p10_top100mrs_ledges/p10_vaq_top100mrs_ledges_LFCh.txt",quote=FALSE,sep="\t")

###############################################################################################
###############################################################################################

## Shadow top 100 pairs

ani_top100mrs_shadowpairs <- ani_mrshadow_10$shadow[which(ani_mrshadow_10$shadow[,1] %in% ani_top100mrs),]
ani_top100mrs_shadowpairs <- unique(as.vector(ani_top100mrs_shadowpairs))

hugo_top100mrs_shadowpairs <- hugo_mrshadow_10$shadow[which(hugo_mrshadow_10$shadow[,1] %in% hugo_top100mrs),]
hugo_top100mrs_shadowpairs <- unique(as.vector(hugo_top100mrs_shadowpairs))

ShA_top100mrs_shadowpairs <- ShA_mrshadow_10$shadow[which(ShA_mrshadow_10$shadow[,1] %in% ShA_top100mrs),]
ShA_top100mrs_shadowpairs <- unique(as.vector(ShA_top100mrs_shadowpairs))

vaq_top100mrs_shadowpairs <- vaq_mrshadow_10$shadow[which(vaq_mrshadow_10$shadow[,1] %in% vaq_top100mrs),]
vaq_top100mrs_shadowpairs <- unique(as.vector(vaq_top100mrs_shadowpairs))

###############################################################################################
###############################################################################################
