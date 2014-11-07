#########################################################################################
################################## M  A  R  I  N  a #####################################
#########################################################################################

######## LOADS ########

library(ssmarina)
library(mixtools)
load("/home/hachepunto/MARINa/general_MARINA.RData")


########  R E G U L O N ######## 
### One proces by TF list

# adj from ARACNe imput data
TFDB_adjfile <- file.path("/home/hachepunto/MARINa/AnimalTFDB_p1e-10.adj")

# Regulon generation

TFDB_regulon <- aracne2regulon(TFDB_adjfile, dset)


########  M A R I N A ######## 

TFDB_mrs <- marina(signature, TFDB_regulon, nullmodel)


# Leading-edge analysis
TFDB_mrs <- ledge(TFDB_mrs)



########  B O O T S T R A P E D   M A R I N A ######## 
TFDB_mrs_BT <- marina(BTsignature, TFDB_regulon, nullmodel)


TFDB_mrs_BT <- bootstrapMarina(mrs, "mode")



#######  S H A D O W and S Y N E R G Y A N A L Y S I S #########

# Shadow analysis
TFDB_mrshadow <- shadow(TFDB_mrs, pval=25)


# Synergy analysis
TFDB_mrs <- marinaCombinatorial(TFDB_mrs, regulators=25) 


# Bootstraped Shadow analysis
TFDB_mrshadow_BT <- shadow(TFDB_mrs_BT, pval=25)


# Bootstraped Synergy analysis
TFDB_mrs_BT <- marinaCombinatorial(TFDB_mrs_BT, regulators=25) 


###################################################################################


##################################################


save.image("/home/hachepunto/MARINa/TFDB_MARINa.RData")
savehistory("/home/hachepunto/MARINa/TFDB_MARINa.Rhistory")


