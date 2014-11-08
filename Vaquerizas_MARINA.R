#########################################################################################
################################## M  A  R  I  N  a #####################################
#########################################################################################

######## LOADS ########

setwd(~/MARINa/)

library(ssmarina)
library(mixtools)
load("general_MARINA.RData")
to_clean <- ls()

print("Source loaded")

########  R E G U L O N ######## 
### One proces by TF list

# adj from ARACNe imput data
Vaquerizas_adjfile <- file.path("/home/hachepunto/MARINa/Vaquerizas_p1e-10.adj")

print("adj loaded")


print("proced to calculet regulon")

Vaquerizas_regulon <- aracne2regulon(Vaquerizas_adjfile, dset)

print("regulon calculated")


print("proced to calculet  master regulators")

########  M A R I N A ######## 

Vaquerizas_mrs <- marina(signature, Vaquerizas_regulon, nullmodel)

# Leading-edge analysis
Vaquerizas_mrs <- ledge(Vaquerizas_mrs)

print("master regulators calculated")
summary(Vaquerizas_mrs)


Vaquerizas_top_100 <-summary(Vaquerizas_mrs,mrs=100)
write.table(Vaquerizas_top_100,file = "Vaquerizas_top_100mrs.txt",quote = FALSE, sep = "\t")

print("Top 100 master regulators saved")


print("proced to calculet  bootstraped master regulators")

########  B O O T S T R A P E D   M A R I N A ######## 

Vaquerizas_BT_mrs <- marina(BTsignature, Vaquerizas_regulon, nullmodel)

Vaquerizas_BT_mrs <- bootstrapMarina(Vaquerizas_BT_mrs , "mode")

print("master regulators bootstraped calculated")
summary(Vaquerizas_BT_mrs)

TFDB_BT_top_100 <-summary(Vaquerizas_BT_mrs,mrs=100)
write.table(TFDB_BT_top_100,file = "TFDB_BT_top_100mrs.txt",quote = FALSE, sep = "\t")

print("Top 100 bootstraped saved")

print("proced to Shadow and Synergy analysis")
print("proced to Shadow analysis")
#######  S H A D O W and S Y N E R G Y A N A L Y S I S #########

# Shadow analysis
Vaquerizas_mrshadow <- shadow(Vaquerizas_mrs, pval=25)

# Bootstraped Shadow analysis
Vaquerizas_mrshadow_BT <- shadow(Vaquerizas_BT_mrs, pval=25)

print("Shadow analysis done")


print("proced to Synergy analysis")

# Synergy analysis
Vaquerizas_mrs <- marinaCombinatorial(Vaquerizas_mrs, regulators=25)


# Bootstraped Synergy analysis
Vaquerizas_BT_mrs <- marinaCombinatorial(Vaquerizas_BT_mrs, regulators=25)

print("Synergy analysis done")


print("proced to Save")
###################################################################################


#################### S A V E  ######################
##################################################

rm(list=to_clean,to_clean) #cleaning

save.image("Vaquerizas_MARINa.RData")
savehistory("Vaquerizas_MARINa.Rhistory")

print("All done")