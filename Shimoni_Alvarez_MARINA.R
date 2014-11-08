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
ShA_adjfile <- file.path("/home/hachepunto/MARINa/Shimoni_Alvarez_p1e-10.adj")
print("adj loaded")


print("proced to calculet regulon")
# Regulon generation

ShA_regulon <- aracne2regulon(ShA_adjfile, dset)


print("regulon calculated")
########  M A R I N A ######## 
print("proced to calculet  master regulators")
ShA_mrs <- marina(signature, ShA_regulon, nullmodel)


# Leading-edge analysis
ShA_mrs <- ledge(ShA_mrs)

print("master regulators calculated")

summary(ShA_mrs)

ShA_top_100 <-summary(ShA_mrs,mrs=100)
write.table(ShA_top_100,file = "ShA_top_100mrs.txt",quote = FALSE, sep = "\t")

print("Top 100 saved")
########  B O O T S T R A P E D   M A R I N A ######## 
print("proced to calculet  bootstraped master regulators")
ShA_BT_mrs <- marina(BTsignature, ShA_regulon, nullmodel)


ShA_BT_mrs <- bootstrapMarina(ShA_BT_mrs, "mode")
print("master regulators bootstraped calculated")
summary(ShA_BT_mrs)

ShA_BT_top_100 <-summary(ShA_BT_mrs,mrs=100)
write.table(ShA_BT_top_100,file = "ShA_BT_top_100mrs.txt",quote = FALSE, sep = "\t")

print("Top 100 bootstraped saved")

print("proced to Shadow and Synergy analysis")
#######  S H A D O W and S Y N E R G Y A N A L Y S I S #########

# Shadow analysis
ShA_mrshadow <- shadow(ShA_mrs, pval=25)
# Bootstraped Shadow analysis
ShA_mrshadow_BT <- shadow(ShA_BT_mrs, pval=25)

print("Shadow analysis done")


print("proced to Synergy analysis")

# Synergy analysis
ShA_mrs <- marinaCombinatorial(ShA_mrs, regulators=25) 

# Bootstraped Synergy analysis
ShA_BT_mrs <- marinaCombinatorial(ShA_BT_mrs, regulators=25) 
print("Synergy analysis done")


print("proced to Save")

###################################################################################


##################################################
rm(list=to_clean,to_clean) #cleaning
save.image("ShA_MARINa.RData")
savehistory("ShA_MARINa.Rhistory")

print("All done")
