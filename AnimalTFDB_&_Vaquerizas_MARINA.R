#########################################################################################
################################## M  A  R  I  N  a #####################################
#########################################################################################

######## LOADS ########

library(ssmarina)
library(mixtools)
load("~/MARINa/general_MARINA.RData")


########  R E G U L O N ######## 
### One proces by TF list

# adj from ARACNe imput data
TFDB_adjfile <- file.path("/home/hachepunto/transfac_network/AnimalTFDB_p1e-10.adj")
Vaquerizas_adjfile <- file.path("/home/hachepunto/transfac_network/Vaquerizas_p1e-10.adj")

# Regulon generation

TFDB_regulon <- aracne2regulon(TFDB_adjfile, dset)
Vaquerizas_regulon <- aracne2regulon(Vaquerizas_adjfile, dset)


########  M A R I N A ######## 

TFDB_mrs <- marina(signature, TFDB_regulon, nullmodel)
Vaquerizas_mrs <- marina(signature, Vaquerizas_regulon, nullmodel)

# Leading-edge analysis
TFDB_mrs <- ledge(TFDB_mrs)
Vaquerizas_mrs <- ledge(Vaquerizas_mrs)


########  B O O T S T R A P E D   M A R I N A ######## 
TFDB_mrs_BT <- marina(BTsignature, TFDB_regulon, nullmodel)
Vaquerizas_mrs_BT <- marina(BTsignature, Vaquerizas_regulon, nullmodel)

TFDB_mrs_BT <- bootstrapMarina(mrs, "mode")
Vaquerizas_mrs_BT <- bootstrapMarina(mrs, "mode")


#######  S H A D O W and S Y N E R G Y A N A L Y S I S #########

# Shadow analysis
TFDB_mrshadow <- shadow(TFDB_mrs, pval=25)
Vaquerizas_mrshadow <- shadow(Vaquerizas_mrs, pval=25)

# Synergy analysis
TFDB_mrs <- marinaCombinatorial(TFDB_mrs, regulators=25) 
Vaquerizas_mrs <- marinaCombinatorial(Vaquerizas_mrs, regulators=25)

# Bootstraped Shadow analysis
TFDB_mrshadow_BT <- shadow(TFDB_mrs_BT, pval=25)
Vaquerizas_mrshadow_BT <- shadow(Vaquerizas_mrs_BT, pval=25)

# Bootstraped Synergy analysis
TFDB_mrs_BT <- marinaCombinatorial(TFDB_mrs_BT, regulators=25) 
Vaquerizas_mrs_BT <- marinaCombinatorial(Vaquerizas_mrs_BT, regulators=25)

###################################################################################


##################################################


save.image("~/MARINa/TFDB&Vaq_MARINa.RData")
savehistory("~/MARINa/TFDB&Vaq_MARINa.Rhistory")















dset_sin_CLIC1_FTL <- dset[-c(6468,8309),]


regulon_sin_CLIC1_FTL <- aracne2regulon(adjfile_ShA, dset_sin_CLIC1_FTL)
### ***Did't work for the no convergent warning*** ###


 #### Signature without problem genes:

sanMuest_sin_CLIC1_FTL <- which(colnames(dset_sin_CLIC1_FTL) %in% grupos[["sanos"]])
enfMuest_sin_CLIC1_FTL <- which(colnames(dset_sin_CLIC1_FTL) %in% grupos[["enfermos"]])
signature_sin_CLIC1_FTL <- rowTtest(dset_sin_CLIC1_FTL[, enfMuest_sin_CLIC1_FTL], dset_sin_CLIC1_FTL[, sanMuest_sin_CLIC1_FTL])
signature_sin_CLIC1_FTL <- (qnorm(signature_sin_CLIC1_FTL$p.value/2, lower.tail=F) * sign(signature_sin_CLIC1_FTL$statistic))[, 1]

#### Generar el modelo nulo sin genes problema:
nullmodel_sin_CLIC1_FTL <- ttestNull(dset_sin_CLIC1_FTL[, enfMuest_sin_CLIC1_FTL], dset_sin_CLIC1_FTL[, sanMuest_sin_CLIC1_FTL], per=1000, repos=T)

#### Generar el objeto marina sin genes problema:
mrs_sin_CLIC1_FTL <- marina(signature_sin_CLIC1_FTL, regulon_sin_CLIC1_FTL, nullmodel_sin_CLIC1_FTL)

#### Respalda el objeto marina sin genes problema:
mrs_sin_CLIC1_FTL_wo_ledge <- mrs_sin_CLIC1_FTL

#### Generar el objeto marina con ledge sin genes problema:
mrs_sin_CLIC1_FTL <- ledge(mrs_sin_CLIC1_FTL)
#### *** FUNCIONÓ *** #####
# ya permite hacer el análisis de Ledge

