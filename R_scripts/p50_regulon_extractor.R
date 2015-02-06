### Extracción del Regulon Vaquerizas

setwd("~/MARINa")
library(ssmarina)
load("p50_results_MARINA.RData")

vaquerizas_regulon_50 <- vaq_mrs_50$regulon

vaq_regulon_DF <- do.call(rbind, lapply(seq_along(vaquerizas_regulon_50), function(i){data.frame(TFs=i, vaquerizas_regulon_50[[i]])}))

r <- unlist(vaq_regulon_DF[1])

vaq_regulon_DF[1] <- names(vaquerizas_regulon_50[r])

head(vaq_regulon_DF, 100)

write.table(vaq_regulon_DF,file="./regulones/vaq_p50_regulon.txt", quote=F)



### Extracción del Regulon Hugo.

setwd("~/MARINa")
library(ssmarina)
load("p50_results_MARINA.RData")

hugo_regulon_50 <- hugo_mrs_50$regulon

hugo_regulon_50_DF <- do.call(rbind, lapply(seq_along(hugo_regulon_50), function(i){data.frame(TFs=i, hugo_regulon_50[[i]], check.names = FALSE)}))

r <- unlist(hugo_regulon_50_DF[1])

hugo_regulon_50_DF[1] <- names(hugo_regulon_50[r])

length(class(row.names(hugo_regulon_50_DF))

write.table(hugo_regulon_50_DF,file="./regulones/hugo_p50_regulon.txt", quote=F)