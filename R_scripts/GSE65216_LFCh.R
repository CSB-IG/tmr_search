###Log Full Change

setwd("~/MARINa")
library(limma)

# Expression data
dset <- read.table('/mnt/c/hachepunto/MARINA_backup/validation/GSE65216/affy_GSE65216_hgu133a_norm.txt',header=TRUE,row.names=1)
dset <- as.matrix(dset)
head(dset)
length(dset)

# Groups list

enfermos <- c("GSM1588970", "GSM1588972", "GSM1588974", "GSM1588975", "GSM1588976", "GSM1588977", "GSM1588979", "GSM1588980", "GSM1588981", "GSM1588982", "GSM1588984", "GSM1588985", "GSM1588986", "GSM1588988", "GSM1588989", "GSM1588991", "GSM1588993", "GSM1588995", "GSM1588997", "GSM1588998", "GSM1589000", "GSM1589002", "GSM1589003", "GSM1589004", "GSM1589005", "GSM1589006", "GSM1589008", "GSM1589009", "GSM1589011", "GSM1589012", "GSM1589013", "GSM1589014", "GSM1589016", "GSM1589017", "GSM1589018", "GSM1589020", "GSM1589021", "GSM1589022", "GSM1589024", "GSM1589025", "GSM1589026", "GSM1589027", "GSM1589028", "GSM1589030", "GSM1589031", "GSM1589033", "GSM1589034", "GSM1589035", "GSM1589037", "GSM1589038", "GSM1589039", "GSM1589041", "GSM1589042", "GSM1589043", "GSM1589044", "GSM1589045", "GSM1589046", "GSM1589047", "GSM1589049", "GSM1589050", "GSM1589051", "GSM1589052", "GSM1589053", "GSM1589054", "GSM1589056", "GSM1589057", "GSM1589058", "GSM1589060", "GSM1589061", "GSM1589062", "GSM1589063", "GSM1589064", "GSM1589065", "GSM1589066", "GSM1589067", "GSM1589068", "GSM1589069", "GSM1589070", "GSM1589071", "GSM1589072", "GSM1589073", "GSM1589074", "GSM1589075", "GSM1589076", "GSM1589077", "GSM1589078", "GSM1589079", "GSM1589080", "GSM1589081", "GSM1589082", "GSM1589083", "GSM1589084", "GSM1589085", "GSM1589086", "GSM1589087", "GSM1589088", "GSM1589089", "GSM1589090", "GSM1589091", "GSM1589092", "GSM1589093", "GSM1589094", "GSM1589095", "GSM1589097", "GSM1589098", "GSM1589099", "GSM1589101", "GSM1589103", "GSM1589105", "GSM1589107", "GSM1589109", "GSM1589110", "GSM1589111", "GSM1589112", "GSM1589113", "GSM1589114", "GSM1589115", "GSM1589116", "GSM1589117", "GSM1589118", "GSM1589119", "GSM1589120", "GSM1589121", "GSM1589122", "GSM1589123", "GSM1589124", "GSM1589125", "GSM1589126", "GSM1589127", "GSM1589128")
sanos <- c("GSM1589130", "GSM1589132", "GSM1589135", "GSM1589136", "GSM1589139", "GSM1589142", "GSM1589144", "GSM1589145", "GSM1589148", "GSM1589150", "GSM1589151")
grupos <- list(sanos=sanos,enfermos=enfermos)

# Expresion matrix separation
sanMuest <- which(colnames(dset) %in% grupos[["sanos"]])
enfMuest <- which(colnames(dset) %in% grupos[["enfermos"]])

design = matrix(rep(0,282), nrow=141)
colnames(design) = c("enf","san")

design[enfMuest,1]=1
design[sanMuest,2]=1

cont.matrix = makeContrasts("enf - san", levels=design)


# Log fold-change

fit  = lmFit(dset, design)
fit2 = contrasts.fit(fit, cont.matrix)
fit2 = eBayes(fit2)

topTable <- topTable(fit2, coef=1, adjust = "fdr",sort.by="none", n = length(rownames(dset)))

head(topTable)
str(topTable)

head(rownames(dset),10)
head(rownames(topTable),10)

write.table(topTable,file="./validation_analisys/LFCh_GSE65216_hgu133a.txt",quote=FALSE,sep="\t")