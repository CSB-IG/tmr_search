#########################################################################################
################################## M  A  R  I  N  a #####################################
#########################################################################################

######## LOADS ########

setwd("~/MARINa")
library(ssmarina)
library(mixtools)
load("general_MARINA.RData")
to_clean <- ls()
print("Source loaded")


########  R E G U L O N ######## 
### One proces by TF list

# adj from ARACNe imput data
TFDB_adjfile <- file.path("/home/hachepunto/MARINa/AnimalTFDB_p1e-10.adj")
print("adj loaded")


print("proced to calculet regulon")
# Regulon generation

TFDB_regulon <- aracne2regulon(TFDB_adjfile, dset)


print("regulon calculated")
########  M A R I N A ######## 
print("proced to calculet  master regulators")
TFDB_mrs <- marina(signature, TFDB_regulon, nullmodel)


# Leading-edge analysis
TFDB_mrs <- ledge(TFDB_mrs)

print("master regulators calculated")

summary(TFDB_mrs)

TFDB_top_100 <-summary(TFDB_mrs,mrs=100)
write.table(TFDB_top_100,file = "TFDB_top_100mrs.txt",quote = FALSE, sep = "\t")

print("Top 100 saved")
########  B O O T S T R A P E D   M A R I N A ######## 
print("proced to calculet  bootstraped master regulators")
TFDB_BT_mrs <- marina(BTsignature, TFDB_regulon, nullmodel)


TFDB_BT_mrs <- bootstrapMarina(TFDB_BT_mrs, "mode")
print("master regulators bootstraped calculated")
summary(TFDB_BT_mrs)

TFDB_BT_top_100 <-summary(TFDB_BT_mrs,mrs=100)
write.table(TFDB_BT_top_100,file = "TFDB_BT_top_100mrs.txt",quote = FALSE, sep = "\t")

print("Top 100 bootstraped saved")

print("proced to Shadow and Synergy analysis")
#######  S H A D O W and S Y N E R G Y A N A L Y S I S #########

# Shadow analysis
TFDB_mrshadow <- shadow(TFDB_mrs, pval=25)
# Bootstraped Shadow analysis
TFDB_mrshadow_BT <- shadow(TFDB_BT_mrs, pval=25)

print("Shadow analysis done")


print("proced to Synergy analysis")

# Synergy analysis
TFDB_mrs <- marinaCombinatorial(TFDB_mrs, regulators=25) 

# Bootstraped Synergy analysis
TFDB_BT_mrs <- marinaCombinatorial(TFDB_BT_mrs, regulators=25) 
print("Synergy analysis done")


print("proced to Save")

###################################################################################


##################################################
rm(list=to_clean,to_clean) #cleaning
save.image("TFDB_MARINa.RData")
savehistory("TFDB_MARINa.Rhistory")

print("All done")
