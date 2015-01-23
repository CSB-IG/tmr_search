

library(limma)
load('Rdatas/dset.RData')

design = matrix(rep(0,1760), nrow=880)
colnames(design) = c('enf','san')

design[enfMuest,1]=1
design[sanMuest,2]=1

cont.matrix = makeContrasts('enf - san', levels=design)

toclean <- ls()

###############################################################################################
### GROUPS:

ani_top100mrs_neighbors <- read.table(file="lists/p25/p25_ani_top100mrs&neighbors.txt", stringsAsFactors = FALSE)
ani_top100mrs_neighbors <- ani_top100mrs_neighbors$V1

hugo_top100mrs_neighbors <- read.table(file="lists/p25/p25_hugo_top100mrs&neighbors.txt", stringsAsFactors = FALSE)
hugo_top100mrs_neighbors <- hugo_top100mrs_neighbors$V1

sha_top100mrs_neighbors <- read.table(file="lists/p25/p25_sha_top100mrs&neighbors.txt", stringsAsFactors = FALSE)
sha_top100mrs_neighbors <- sha_top100mrs_neighbors$V1


vaq_top100mrs_neighbors <- read.table(file="lists/p25/p25_vaq_top100mrs&neighbors.txt", stringsAsFactors = FALSE)
vaq_top100mrs_neighbors <- vaq_top100mrs_neighbors$V1


all_top100mrs_neighbors <- unique(c(ani_top100mrs_neighbors,vaq_top100mrs_neighbors,hugo_top100mrs_neighbors,sha_top100mrs_neighbors))


# Cuenta falsos:
table(sha_top100mrs_neighbors %in% row.names(dset))["FALSE"]
table(ani_top100mrs_neighbors %in% row.names(dset))["FALSE"]
table(hugo_top100mrs_neighbors %in% row.names(dset))["FALSE"]
table(vaq_top100mrs_neighbors %in% row.names(dset))["FALSE"]

table(all_top100mrs_neighbors %in% row.names(dset))["FALSE"]

# Quita falsos "ZNF816"

ani_top100mrs_neighbors<-ani_top100mrs_neighbors[which(ani_top100mrs_neighbors %in% row.names(dset), arr.ind = FALSE, useNames = TRUE)]
hugo_top100mrs_neighbors<-hugo_top100mrs_neighbors[which(hugo_top100mrs_neighbors %in% row.names(dset), arr.ind = FALSE, useNames = TRUE)]
sha_top100mrs_neighbors<-sha_top100mrs_neighbors[which(sha_top100mrs_neighbors %in% row.names(dset), arr.ind = FALSE, useNames = TRUE)]
vaq_top100mrs_neighbors<-vaq_top100mrs_neighbors[which(vaq_top100mrs_neighbors %in% row.names(dset), arr.ind = FALSE, useNames = TRUE)]

###############################################################################################

exp_ani_top100mrs_neighbors <- dset[ani_top100mrs_neighbors,]

ani_top100mrs_neighborsfit  = lmFit(exp_ani_top100mrs_neighbors, design)
ani_top100mrs_neighborsfit2 = contrasts.fit(ani_top100mrs_neighborsfit, cont.matrix)
ani_top100mrs_neighborsfit2 = eBayes(ani_top100mrs_neighborsfit2)

ani_top100mrs_neighborstopTable <- topTable(ani_top100mrs_neighborsfit2, coef=1, adjust = 'fdr', n = length(ani_top100mrs_neighbors))

write.table(ani_top100mrs_neighborstopTable,file='./lists/p25/p25_ani_top100mrs&neighbors_LfCh.txt',quote=FALSE,sep='\t')
###############################################################################################

exp_hugo_top100mrs_neighbors <- dset[hugo_top100mrs_neighbors,]

hugo_top100mrs_neighborsfit  = lmFit(exp_hugo_top100mrs_neighbors, design)
hugo_top100mrs_neighborsfit2 = contrasts.fit(hugo_top100mrs_neighborsfit, cont.matrix)
hugo_top100mrs_neighborsfit2 = eBayes(hugo_top100mrs_neighborsfit2)

hugo_top100mrs_neighborstopTable <- topTable(hugo_top100mrs_neighborsfit2, coef=1, adjust = 'fdr', n = length(hugo_top100mrs_neighbors))

write.table(hugo_top100mrs_neighborstopTable,file='./lists/p25/p25_hugo_top100mrs&neighbors_LfCh.txt',quote=FALSE,sep='\t')
###############################################################################################


exp_sha_top100mrs_neighbors <- dset[sha_top100mrs_neighbors,]

sha_top100mrs_neighborsfit  = lmFit(exp_sha_top100mrs_neighbors, design)
sha_top100mrs_neighborsfit2 = contrasts.fit(sha_top100mrs_neighborsfit, cont.matrix)
sha_top100mrs_neighborsfit2 = eBayes(sha_top100mrs_neighborsfit2)

sha_top100mrs_neighborstopTable <- topTable(sha_top100mrs_neighborsfit2, coef=1, adjust = 'fdr', n = length(sha_top100mrs_neighbors))

write.table(sha_top100mrs_neighborstopTable,file='./lists/p25/p25_sha_top100mrs&neighbors_LfCh.txt',quote=FALSE,sep='\t')
###############################################################################################

exp_vaq_top100mrs_neighbors <- dset[vaq_top100mrs_neighbors,]

vaq_top100mrs_neighborsfit  = lmFit(exp_vaq_top100mrs_neighbors, design)
vaq_top100mrs_neighborsfit2 = contrasts.fit(vaq_top100mrs_neighborsfit, cont.matrix)
vaq_top100mrs_neighborsfit2 = eBayes(vaq_top100mrs_neighborsfit2)

vaq_top100mrs_neighborstopTable <- topTable(vaq_top100mrs_neighborsfit2, coef=1, adjust = 'fdr', n = length(vaq_top100mrs_neighbors))

write.table(vaq_top100mrs_neighborstopTable,file='./lists/p25/p25_vaq_top100mrs&neighbors_LfCh.txt',quote=FALSE,sep='\t')
###############################################################################################

exp_all_top100mrs_neighbors <- dset[all_top100mrs_neighbors,]

all_top100mrs_neighborsfit  = lmFit(exp_all_top100mrs_neighbors, design)
all_top100mrs_neighborsfit2 = contrasts.fit(all_top100mrs_neighborsfit, cont.matrix)
all_top100mrs_neighborsfit2 = eBayes(all_top100mrs_neighborsfit2)

all_top100mrs_neighborstopTable <- topTable(all_top100mrs_neighborsfit2, coef=1, adjust = 'fdr', n = length(all_top100mrs_neighbors))

write.table(all_top100mrs_neighborstopTable,file='./lists/p25/p25_all_top100mrs&neighbors_LfCh.txt',quote=FALSE,sep='\t')
###############################################################################################