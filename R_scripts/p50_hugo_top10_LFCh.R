
setwd("~/MARINa")
library(ssmarina)
library(limma)

load("p50_results_MARINA.RData")
load("Rdatas/dset.RData")

design = matrix(rep(0,1760), nrow=880)
colnames(design) = c("enf","san")

design[enfMuest,1]=1
design[sanMuest,2]=1

cont.matrix = makeContrasts("enf - san", levels=design)

hugo_top10mrs <- summary(hugo_mrshadow_50,mrs=10)$MARINA.results$Regulon
hugo_top10mrs <- as.vector(hugo_top10mrs)

hugo_top10_exp <- dset[hugo_top10mrs,]

hugo_lg_fit  = lmFit(hugo_top10_exp, design)
hugo_lg_fit2 = contrasts.fit(hugo_lg_fit, cont.matrix)
hugo_lg_fit2 = eBayes(hugo_lg_fit2)

hugo_lg_top10_Table <- topTable(hugo_lg_fit2, coef=1, adjust = "fdr", n = length(hugo_top10mrs))

write.table(hugo_lg_top10_Table,file="./top10_logfull-Ch.txt",quote=FALSE,sep="\t")
