############ ShA 50 ###########

######## LOADS ########

setwd("~/MARINa")
library(ssmarina)
library(mixtools)
load("./Rdatas/general_MARINa.RData")
load("./Rdatas/ShA_regulon_25.RData")
to_clean <- ls()
print("Source and data set loaded")

########  M A R I N A ######## 
print("proced to compute the master regulators")
ShA_mrs_25 <- marina(signature, ShA_regulon_25, nullmodel)


# Leading-edge analysis
ShA_mrs_25 <- ledge(ShA_mrs_25)

print("master regulators computed")

print(summary(ShA_mrs_25))

ShA_top100_25 <-summary(ShA_mrs_25,mrs=100)
write.table(ShA_top100_25,file = "./lists/ShA_top100mrs_25.txt",quote = FALSE, sep = "\t")

print("Top 100 Master Regulators saved in 'ShA_top100mrs_25.txt' file")

print("proced to Save (mrs)")
save.image("ShA_MARINa_25.RData")

#######  S H A D O W and S Y N E R G Y A N A L Y S I S #########
print("proced to run Shadow and Synergy analysis")

print("proced to compute the Shadow analysis")
# Shadow analysis
ShA_mrshadow_25 <- shadow(ShA_mrs_25, pval=25)

print("Shadow analysis done")

print("proced to Save (shadow)")
save.image("ShA_MARINa_25.RData")


print("proced to compute Synergy analysis")

# Synergy analysis
ShA_mrs_25 <- marinaCombinatorial(ShA_mrs_25, regulators=25)

######### F I N A L  S A V E  ##########

print("proced to Save")
rm(list=to_clean,to_clean) #cleaning
save.image("ShA_MARINa_25.RData")

print("All done!")