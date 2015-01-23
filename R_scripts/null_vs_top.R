


library(ssmarina)
library(limma)

load("p75_results_MARINA.RData")
load('Rdatas/dset.RData')

hugo_top100mrs_75 <- as.vector(summary(hugo_mrshadow_75,mrs=100)$MARINA.results$Regulon)

design = matrix(rep(0,1760), nrow=880)
colnames(design) = c('enf','san')

design[enfMuest,1]=1
design[sanMuest,2]=1

cont.matrix = makeContrasts('enf - san', levels=design)

toclean <- ls()


hugo_neighbors <- read.table(file="lists/p25/p25_hugo_neighbors.txt", stringsAsFactors = FALSE)
hugo_neighbors <- hugo_neighbors$V1

hugo_p25_null <- sample(hugo_neighbors,100)

write.table(hugo_p25_null,file='./lists/p25/p25_hugo_null_list.txt',quote=FALSE,sep='\t',row.names=F,col.names=F)





TF_p75_hugo <- read.table(file="./lists/p75/p75_hugo_TFs.txt", stringsAsFactors = FALSE)
TF_p75_hugo <- TF_p75_hugo$V1

hugo_null100_p75 <- sample(TF_p75_hugo,100)

write.table(hugo_null100_p75,file='./lists/p75/p75_hugo_null100_list.txt',quote=FALSE,sep='\t',row.names=F,col.names=F)

compare()





hugo_top20mrs_75 <- as.vector(summary(hugo_mrshadow_75,mrs=20)$MARINA.results$Regulon)

write.table(hugo_top20mrs_75,file='./lists/p75/p75_hugo_top20.txt',quote=FALSE,sep='\t',row.names=F,col.names=F)

hugo_null20_p75 <- sample(TF_p75_hugo,20)

write.table(hugo_null20_p75,file='./lists/p75/p75_hugo_null20_list.txt',quote=FALSE,sep='\t',row.names=F,col.names=F)


hugo_top10mrs_75 <- as.vector(summary(hugo_mrshadow_75,mrs=10)$MARINA.results$Regulon)

write.table(hugo_top10mrs_75,file='./lists/p75/p75_hugo_top10.txt',quote=FALSE,sep='\t',row.names=F,col.names=F)

hugo_null10_p75 <- sample(TF_p75_hugo,10)

write.table(hugo_null10_p75,file='./lists/p75/p75_hugo_null10_list.txt',quote=FALSE,sep='\t',row.names=F,col.names=F)