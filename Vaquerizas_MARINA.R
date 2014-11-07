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
Vaquerizas_adjfile <- file.path("/home/hachepunto/MARINa/Vaquerizas_p1e-10.adj")

# Regulon generation

Vaquerizas_regulon <- aracne2regulon(Vaquerizas_adjfile, dset)


########  M A R I N A ######## 

Vaquerizas_mrs <- marina(signature, Vaquerizas_regulon, nullmodel)

# Leading-edge analysis
Vaquerizas_mrs <- ledge(Vaquerizas_mrs)


########  B O O T S T R A P E D   M A R I N A ######## 

Vaquerizas_mrs_BT <- marina(BTsignature, Vaquerizas_regulon, nullmodel)

Vaquerizas_mrs_BT <- bootstrapMarina(mrs, "mode")


#######  S H A D O W and S Y N E R G Y A N A L Y S I S #########

# Shadow analysis
Vaquerizas_mrshadow <- shadow(Vaquerizas_mrs, pval=25)

# Synergy analysis
Vaquerizas_mrs <- marinaCombinatorial(Vaquerizas_mrs, regulators=25)

# Bootstraped Shadow analysis
Vaquerizas_mrshadow_BT <- shadow(Vaquerizas_mrs_BT, pval=25)

# Bootstraped Synergy analysis
Vaquerizas_mrs_BT <- marinaCombinatorial(Vaquerizas_mrs_BT, regulators=25)

###################################################################################


#################### S A V E  ######################
##################################################


save.image("home/hachepunto/MARINa/Vaquerizas_MARINa.RData")
savehistory("home/hachepunto/MARINa/Vaquerizas_MARINa.Rhistory")