

###### Log fold-change of the top 100 MTF and null model ######

# Vector con los nombres de los top 100 genes.

library(limma)
ls()
load("general_MARINa.RData")

top100 <- c("IRF9", "IRF7", "SREBF1", "NCOR1", "ZBTB17", "E4F1", "ZKSCAN4", "CREB3L2", "JUNB", "FOXE1", "HNF1A", "FOXL1", "POU4F1", "ZNF157", "MEOX2", "HESX1", "ZNF80", "SOX3", "ZFX", "POU2F3", "GMEB2", "ZNF668", "ZNF142", "DEAF1", "ATF4", "FOXO3", "TFCP2", "FOXK2", "PITX1", "HOXA1", "DMRT1", "SOX30", "ST18", "RREB1", "ZKSCAN5", "RFX4", "ESR2", "SATB2", "MSX2", "PBRM1", "SIM1", "MYT1L", "ZNF506", "ZNF586", "MYB", "ZNF419", "ZNF780B", "NPAS3", "ZNF667", "ZBTB43", "ZBED4", "DR1", "ZNF234", "ZNF674", "PHB2", "CIZ1", "VPS72", "PDLIM1", "SUPT6H", "REXO4", "ZNF428", "HCFC1", "CSDE1", "NOTCH1", "ZBED1", "SUB1", "ZNF43", "NCOA1", "CEBPB", "CREB3", "ZNF175", "ZNF135", "ZNF33B", "DLX2", "HLF", "ZNF639", "IRX4", "ZNF334", "ZNF254", "BNC1", "ZNF682", "ZNF225", "ZNF613", "CAMTA1", "MYST2", "NEUROG2", "SCAND2", "TOX", "CHD7", "DAXX", "SKIL", "ZNF345", "BLZF1", "TAF9B", "PIAS2", "PTH", "TIGD6", "ALS2CR8", "TLE4", "C2orf3", "SLC11A1", "INHBA", "ZDHHC7", "PARP1", "HLA-A", "PLXND1", "CPSF4", "TRRAP", "SIGIRR", "TICAM1", "ZMYM6", "JMJD6", "ASXL1", "RING1", "ACP5", "DAP", "C14orf169", "HIF1AN", "PNN", "HLA-C", "SMARCD2", "PPRC1", "RPS19", "TRIB3", "MKRN1", "MYH9", "PYCARD", "RBM9", "SUPT16H", "COBRA1", "PDCD11", "TOE1", "FLII", "HIST1H1C", "HMGN2", "C17orf81", "CTNNBIP1", "DNMT1", "ENG", "HIST1H3G", "KDM1", "LASS2", "MST1R", "TERF2IP", "HIST1H2BD", "MATR3", "FLNA", "NOLC1", "ZC3H11A", "TRIM28", "HELZ", "TRIB1", "PHF21A", "SPG7", "TAF10", "ZNF593", "PFDN5", "BOLA1", "RPL27", "UBR4", "DUSP12", "TPT1", "CSDA", "HNRNPK", "NPAS2", "HDAC5", "IRF3", "SP110", "USP7", "MBD1", "THAP7", "ELF1", "HMGB1", "TGIF1", "ZNF212", "RFX5", "ZNF410", "LHX2", "TCF21", "ZNF460", "HOXB6", "HOXB9", "IKZF5", "ZNF41", "AFF2", "ZNF493", "CSRNP3", "ZNF669", "ZNF549", "ZNF471", "KDM4C", "LOC100130354", "WDHD1", "BRPF1", "MKRN3", "ZFR", "SNAPC5", "YOD1", "ZNF614", "RAG1", "ZMAT3", "ZRSR1", "IL2", "ELSPBP1", "KCNH6", "IKBKAP", "RQCD1", "RPA4", "IL17A", "SLC4A10", "ZFAND5", "H2AFJ", "HIST1H2AL", "USP22", "EREG", "FUT1", "IL11", "MSTN", "HTATIP2", "ZNF185", "C11orf30", "SUV39H2", "BRWD1", "RNF2", "RCHY1", "TMF1")

exp_top100 <- dset[top100,]

design = matrix(rep(0,1760), nrow=880)
colnames(design) = c("enf","san")

design[enfMuest,1]=1
design[sanMuest,2]=1

cont.matrix = makeContrasts("enf - san", levels=design)

fit  = lmFit(exp_top100, design)
fit2 = contrasts.fit(fit, cont.matrix)
fit2 = eBayes(fit2)

contrastes = ("enf vs san")

topTable_top100 <- topTable(fit2, coef=1, adjust = "fdr", n = 226)



############################### Null 226 ###############################

genesym_dset <- row.names(dset[-names_at,])

nullM <- sample(genesym_dset,226)
nullmodel <- dset[nullM,] 


null_design = matrix(rep(0,1760), nrow=880)
colnames(null_design) = c("enf_null","san_null")

null_design[enfMuest,1]=1
null_design[sanMuest,2]=1

null_cont.matrix = makeContrasts("enf_null - san_null", levels=null_design)

null_fit  = lmFit(nullmodel, null_design)
null_fit2 = contrasts.fit(null_fit, null_cont.matrix)
null_fit2 = eBayes(null_fit2)

null_contrastes = ("enf_null vs san_null")

null_topTable <- topTable(null_fit2, coef=1, adjust = "fdr", n = 226)


head(null_topTable)
head(topTable_top100)


##### Volcano plot ######
#volcanoplot(fit2, coef=1, main = contrastes[1],  highlight=100, names=top100)
#volcanoplot(null_fit2, coef=1, main = null_contrastes[1],  highlight=100, names=nullM)



write.table(topTable_top100,file="log_f-Ch_top100.txt",quote=FALSE,sep="\t")
write.table(null_topTable,file="log_f-Ch_nullModel_226.txt",quote=FALSE,sep="\t")


############################### Null 100 ###############################

nullM_100 <- sample(genesym_dset,100)
nullmodel_100 <- dset[nullM,] 

null_100_fit  = lmFit(nullmodel_100, design)
null_100_fit2 = contrasts.fit(null_100_fit, cont.matrix)
null_100_fit2 = eBayes(null_100_fit2)

null_100_topTable <- topTable(null_100_fit2, coef=1, adjust = "fdr", n = 100)

write.table(null_100_topTable,file="log_f-Ch_nullModel_100.txt",quote=FALSE,sep="\t")


###########################################################################
############################### Top 100 by list ###############################
###########################################################################


############################### VAQUERIZAS ###############################

vaq_names_t100 <- c("GMEB2", "ZNF668", "ZNF142", "PHB2", "CIZ1", "PLXND1", "CPSF4", "VPS72", "IRF9", "IRF7", "DEAF1", "SREBF1", "ZBTB17", "E4F1", "MKRN1", "ZNF428", "NCOR1", "TOE1", "ZKSCAN4", "CREB3L2", "HIST1H1C", "CSDE1", "TERF2IP", "ZBED1", "MST1R", "MATR3", "ZC3H11A", "ATF4", "FOXO3", "HELZ", "TFCP2", "FOXK2", "ZNF593", "JUNB", "ZNF43", "PHF21A", "UBR4", "BOLA1", "NCOA1", "DUSP12", "PITX1", "HOXA1", "ZNF334", "ZNF254", "DMRT1", "WDHD1", "SOX30", "ZNF682", "MKRN3", "BRPF1", "ZNF225", "ST18", "BNC1", "ZNF613", "ZFR", "ZKSCAN5", "RREB1", "CAMTA1", "RFX4", "SATB2", "ESR2", "YOD1", "PBRM1", "RAG1", "ZMAT3", "MSX2", "ZRSR1", "SIM1", "FOXE1", "NEUROG2", "MYT1L", "ZNF506", "ZNF586", "MYB", "ZNF419", "HNF1A", "TOX", "CHD7", "ZNF780B", "NPAS3", "RPA4", "FOXL1", "ZNF667", "SKIL", "POU4F1", "ZNF157", "HESX1", "MEOX2", "SLC4A10", "ZBTB43", "ZNF345", "ZBED4", "ZNF80", "SOX3", "TIGD6", "DR1", "ZNF234", "ZNF674", "ZFX", "POU2F3")
bt_vaq_names_t100 <- c("GMEB2", "ZNF668", "ZNF142", "PHB2", "CIZ1", "PLXND1", "CPSF4", "VPS72", "IRF7", "IRF9", "DEAF1", "E4F1", "SREBF1", "ZBTB17", "NCOR1", "ZNF428", "MKRN1", "HIST1H1C", "TOE1", "ZBED1", "ZKSCAN4", "TERF2IP", "CREB3L2", "MST1R", "CSDE1", "FOXK2", "ATF4", "MATR3", "ZC3H11A", "JUNB", "FOXO3", "HELZ", "ZNF43", "TFCP2", "ZNF593", "PHF21A", "UBR4", "DUSP12", "BOLA1", "ADAR", "NCOA1", "CREB3", "BRPF1", "ZNF639", "BNC1", "WDHD1", "ZNF254", "ST18", "SOX30", "ZNF682", "ZNF225", "DMRT1", "ZFR", "ZNF613", "MKRN3", "ZKSCAN5", "RREB1", "RFX4", "ZRSR1", "NEUROG2", "YOD1", "CAMTA1", "ZNF506", "PBRM1", "ESR2", "MYT1L", "FOXE1", "SATB2", "SIM1", "RAG1", "ZMAT3", "MSX2", "ZNF419", "ZNF586", "HNF1A", "CHD7", "MYB", "ZNF780B", "FOXL1", "NPAS3", "TOX", "RPA4", "SKIL", "ZNF667", "ZNF157", "ZBTB43", "POU4F1", "SLC4A10", "MEOX2", "HESX1", "ZNF345", "ZBED4", "ZNF80", "SOX3", "DR1", "TIGD6", "ZNF674", "ZNF234", "ZFX", "POU2F3")

vaq_top100 <- dset[vaq_names_t100,]
bt_vaq_top100 <- dset[bt_vaq_names_t100,]

vaq_fit  = lmFit(vaq_top100, design)
vaq_fit2 = contrasts.fit(vaq_fit, cont.matrix)
vaq_fit2 = eBayes(vaq_fit2)

bt_vaq_fit  = lmFit(bt_vaq_top100, design)
bt_vaq_fit2 = contrasts.fit(bt_vaq_fit, cont.matrix)
bt_vaq_fit2 = eBayes(bt_vaq_fit2)


vaq_topTable <- topTable(vaq_fit2, coef=1, adjust = "fdr", n = 226)
bt_vaq_topTable <- topTable(bt_vaq_fit2, coef=1, adjust = "fdr", n = 226)

# volcanoplot(vaq_fit2, coef=1, main = contrastes[1],  highlight=100, names=vaq_names_t100)
# volcanoplot(bt_vaq_fit2, coef=1, main = contrastes[1],  highlight=100, names=bt_vaq_names_t100)

write.table(vaq_topTable,file="vaq_l-f-Ch_top100.txt",quote=FALSE,sep="\t")
write.table(bt_vaq_topTable,file="bt_vaq_l-f-Ch_top100.txt",quote=FALSE,sep="\t")


############################### ANIMAL TFDB ###############################

TFDB_names_top100 <- c("GMEB2", "ZNF668", "ZNF142", "IRF9", "IRF7", "DEAF1", "E4F1", "SREBF1", "ZBTB17", "ZKSCAN4", "NCOR1", "CSDE1", "CREB3L2", "ZBED1", "ATF4", "FOXO3", "SUB1", "TFCP2", "FOXK2", "ZNF43", "PITX1", "CREB3", "JUNB", "NPAS2", "CEBPB", "SP110", "IRF3", "MBD1", "THAP7", "TGIF1", "ZNF212", "ELF1", "HMGB1", "ZNF410", "RFX5", "TCF21", "LHX2", "ZNF460", "ZNF41", "HOXB9", "HOXB6", "ZNF175", "IKZF5", "AFF2", "ZNF135", "ZNF493", "ZNF33B", "ZNF669", "ZNF549", "HLF", "DLX2", "ZNF471", "ZNF639", "IRX4", "HOXA1", "DMRT1", "ZNF334", "ZNF254", "ST18", "SOX30", "ZNF225", "ZNF682", "ZNF613", "CAMTA1", "RREB1", "ZKSCAN5", "RFX4", "SATB2", "ESR2", "ZNF614", "PBRM1", "SIM1", "MSX2", "FOXE1", "MYT1L", "ZNF586", "ZNF506", "ZNF419", "MYB", "NEUROG2", "HNF1A", "TOX", "ZNF780B", "NPAS3", "FOXL1", "ZNF667", "POU4F1", "ZBTB43", "HESX1", "MEOX2", "ZNF157", "ZNF345", "PIAS2", "ZBED4", "ZNF80", "SOX3", "ZNF674", "ZNF234", "ZFX", "POU2F3")

TFDB_top100 <- dset[TFDB_names_top100,]

TFDB_fit  = lmFit(TFDB_top100, design)
TFDB_fit2 = contrasts.fit(TFDB_fit, cont.matrix)
TFDB_fit2 = eBayes(TFDB_fit2)

TFDB_topTable <- topTable(TFDB_fit2, coef=1, adjust = "fdr", n = 100)

write.table(TFDB_topTable,file="TFDB_l-f-Ch_top100.txt",quote=FALSE,sep="\t")

############################### Shimony & Álvarez ###############################

ShA_names_top100 <- c("GMEB2", "ZNF668", "SLC11A1", "INHBA", "ZDHHC7", "ZNF142", "PHB2", "CIZ1", "TRRAP", "VPS72", "SIGIRR", "ZMYM6", "TICAM1", "JMJD6", "NCOR1", "ASXL1", "SREBF1", "RING1", "DEAF1", "HIF1AN", "C14orf169", "PNN", "SMARCD2", "IRF9", "IRF7", "SUPT6H", "PPRC1", "ZBTB17", "TRIB3", "PDLIM1", "REXO4", "E4F1", "PYCARD", "SUPT16H", "CREB3L2", "COBRA1", "ZNF428", "HMGN2", "FLII", "ZKSCAN4", "HCFC1", "JUNB", "C17orf81", "ENG", "CTNNBIP1", "NOTCH1", "IL2", "PBRM1", "FOXE1", "ZNF586", "ZNF506", "KCNH6", "MYST2", "ZNF419", "RQCD1", "IKBKAP", "CHD7", "ZNF780B", "SCAND2", "DAXX", "HNF1A", "NPAS3", "POU4F1", "FOXL1", "ZNF667", "SKIL", "ZNF157", "MEOX2", "HESX1", "ZBTB43", "IL17A", "ZFAND5", "TAF9B", "USP22", "BLZF1", "PTH", "PIAS2", "EREG", "IL11", "ZBED4", "ZNF80", "SOX3", "TIGD6", "DR1", "MSTN", "ALS2CR8", "HTATIP2", "ZNF185", "C11orf30", "TLE4", "SUV39H2", "ZNF234", "ZNF674", "RNF2", "BRWD1", "ZFX", "POU2F3", "RCHY1", "C2orf3", "TMF1")

ShA_top100 <- dset[ShA_names_top100,]

ShA_fit  = lmFit(ShA_top100, design)
ShA_fit2 = contrasts.fit(ShA_fit, cont.matrix)
ShA_fit2 = eBayes(ShA_fit2)

ShA_topTable <- topTable(ShA_fit2, coef=1, adjust = "fdr", n = 100)

write.table(ShA_topTable,file="ShA_l-f-Ch_top100.txt",quote=FALSE,sep="\t")

############################### Hugo ###############################

Hugo_names_top100 <- c("PARP1", "HLA-A", "SREBF1", "DAP", "IRF9", "ACP5", "NCOR1", "IRF7", "HLA-C", "PDLIM1", "RPS19", "ZBTB17", "E4F1", "RBM9", "ZKSCAN4", "SUPT6H", "MYH9", "PDCD11", "HCFC1", "REXO4", "JUNB", "CREB3L2", "DNMT1", "HIST1H3G", "KDM1", "NOTCH1", "LASS2", "HIST1H2BD", "NOLC1", "ATF4", "FLNA", "FOXO3", "TRIM28", "TRIB1", "CEBPB", "SUB1", "TFCP2", "FOXK2", "TAF10", "SPG7", "PFDN5", "NCOA1", "RPL27", "CREB3", "TPT1", "PITX1", "CSDA", "HNRNPK", "HDAC5", "USP7", "ZNF33B", "DLX2", "ZNF135", "ZNF175", "CSRNP3", "HLF", "ZNF639", "HOXA1", "IRX4", "LOC100130354", "KDM4C", "SOX30", "BNC1", "ST18", "DMRT1", "RREB1", "RFX4", "SATB2", "ESR2", "ZKSCAN5", "SNAPC5", "MSX2", "SIM1", "MYT1L", "FOXE1", "ELSPBP1", "MYST2", "MYB", "SCAND2", "HNF1A", "DAXX", "FOXL1", "POU4F1", "ZNF157", "MEOX2", "HESX1", "H2AFJ", "HIST1H2AL", "BLZF1", "TAF9B", "PTH", "FUT1", "ZNF80", "SOX3", "ALS2CR8", "DR1", "TLE4", "ZFX", "POU2F3", "C2orf3")

Hugo_top100 <- dset[Hugo_names_top100,]

Hugo_fit  = lmFit(Hugo_top100, design)
Hugo_fit2 = contrasts.fit(Hugo_fit, cont.matrix)
Hugo_fit2 = eBayes(Hugo_fit2)

Hugo_topTable <- topTable(Hugo_fit2, coef=1, adjust = "fdr", n = 100)

write.table(Hugo_topTable,file="Hugo_l-f-Ch_top100.txt",quote=FALSE,sep="\t")

##############################################################
######################### S A V  I N G ############################
##############################################################

save.image("/Users/hachepunto/GitHub/tmr_search/results/top100_bioserch/log_fold_change.RData")