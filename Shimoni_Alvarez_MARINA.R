#########################################################################################
################################## M  A  R  I  N  a #####################################
#########################################################################################

######## LOADS ########

setwd("~/MARINa")
library(ssmarina)
library(mixtools)
load("general_MARINa.RData")
to_clean <- ls()
print("Source and data set loaded")


########  R E G U L O N ######## 
### One proces by TF list

# adj from ARACNe imput data
ShA_adjfile <- file.path("/home/hachepunto/MARINa/Shimoni_Alvarez_p1e-10.adj")
print("adj loaded")


# Regulon generation
print("proced to compute the regulon")
ShA_regulon <- aracne2regulon(ShA_adjfile, dset)


print("regulon computed")

########  M A R I N A ######## 
print("proced to compute the master regulators")
ShA_mrs <- marina(signature, ShA_regulon, nullmodel)


# Leading-edge analysis
ShA_mrs <- ledge(ShA_mrs)

print("master regulators computed")

print(summary(ShA_mrs))

ShA_top_100 <-summary(ShA_mrs,mrs=100)
write.table(ShA_top_100,file = "ShA_top_100mrs.txt",quote = FALSE, sep = "\t")

print("Top 100 Master Regulators saved in 'ShA_top_100mrs.txt' file")

########  B O O T S T R A P E D   M A R I N A ######## 

print("proced to compute the bootstraped master regulators")
ShA_BT_mrs <- marina(BTsignature, ShA_regulon, nullmodel)


ShA_BT_mrs <- bootstrapMarina(ShA_BT_mrs, "mode")

print("master regulators bootstraped calculated")

print(summary(ShA_BT_mrs))

ShA_BT_top_100 <-summary(ShA_BT_mrs,mrs=100)
write.table(ShA_BT_top_100,file = "ShA_BT_top_100mrs.txt",quote = FALSE, sep = "\t")

print("Top 100 bootstraped Master Regulators saved in 'ShA_BT_top_100mrs.txt' file")

#######  S H A D O W and S Y N E R G Y A N A L Y S I S #########
print("proced to run Shadow and Synergy analysis")

print("proced to compute the Shadow analysis")
# Shadow analysis
ShA_mrshadow <- shadow(ShA_mrs, pval=25)
# Bootstraped Shadow analysis
ShA_mrshadow_BT <- shadow(ShA_BT_mrs, pval=25)

print("Shadow analysis done")


print("proced to compute Synergy analysis")

# Synergy analysis
ShA_mrs <- marinaCombinatorial(ShA_mrs, regulators=25) 

# Bootstraped Synergy analysis
ShA_BT_mrs <- marinaCombinatorial(ShA_BT_mrs, regulators=25) 
print("Synergy analysis done")



##################################################
#################### S A V E  ######################
##################################################
print("proced to Save")
rm(list=to_clean,to_clean) #cleaning
save.image("ShA_MARINa.RData")
savehistory("ShA_MARINa.Rhistory")

print("All done!")
