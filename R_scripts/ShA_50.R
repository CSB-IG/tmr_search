############ ShA 50 ###########

######## LOADS ########

setwd("~/MARINa")
library(ssmarina)
library(mixtools)
load("./Rdatas/general_MARINa.RData")
load("./Rdatas/ShA_regulon_50.RData")
to_clean <- ls()
print("Source and data set loaded")

########  M A R I N A ######## 
print("proced to compute the master regulators")
ShA_mrs_50 <- marina(signature, ShA_regulon_50, nullmodel)


# Leading-edge analysis
ShA_mrs_50 <- ledge(ShA_mrs_50)

print("master regulators computed")

print(summary(ShA_mrs_50))

ShA_top100_50 <-summary(ShA_mrs_50,mrs=100)
write.table(ShA_top100_50,file = "./lists/ShA_top100mrs_50.txt",quote = FALSE, sep = "\t")

print("Top 100 Master Regulators saved in 'ShA_top100mrs_50.txt' file")

print("proced to Save (mrs)")
save.image("ShA_MARINa_50.RData")

#######  S H A D O W and S Y N E R G Y A N A L Y S I S #########
print("proced to run Shadow and Synergy analysis")

print("proced to compute the Shadow analysis")
# Shadow analysis
ShA_mrshadow_50 <- shadow(ShA_mrs_50, pval=25)

print("Shadow analysis done")

print("proced to Save (shadow)")
save.image("ShA_MARINa_50.RData")


print("proced to compute Synergy analysis")

# Synergy analysis
ShA_mrs_50 <- marinaCombinatorial(ShA_mrs_50, regulators=25)

######### F I N A L  S A V E  ##########

print("proced to Save")
rm(list=to_clean,to_clean) #cleaning
save.image("ShA_MARINa_50.RData")

print("All done!")