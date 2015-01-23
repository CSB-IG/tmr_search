############ ShA2 ###########

######## LOADS ########

setwd("~/MARINa")
library(ssmarina)
library(mixtools)
load("./Rdatas/general_MARINa.RData")
load("./Rdatas/ShA_regulon_100.RData")
to_clean <- ls()
print("Source and data set loaded")

########  M A R I N A ######## 
print("proced to compute the master regulators")
ShA_mrs_100 <- marina(signature, ShA_regulon_100, nullmodel)


# Leading-edge analysis
ShA_mrs_100 <- ledge(ShA_mrs_100)

print("master regulators computed")

print(summary(ShA_mrs_100))

ShA_top_100 <-summary(ShA_mrs_100,mrs=100)
write.table(ShA_top_100,file = "ShA_top_100mrs_100.txt",quote = FALSE, sep = "\t")

print("Top 100 Master Regulators saved in 'ShA_top_100mrs_100.txt' file")

print("proced to Save (mrs)")
save.image("ShA_MARINa_100.RData")

#######  S H A D O W and S Y N E R G Y A N A L Y S I S #########
print("proced to run Shadow and Synergy analysis")

print("proced to compute the Shadow analysis")
# Shadow analysis
ShA_mrshadow_100 <- shadow(ShA_mrs_100, pval=25)

print("Shadow analysis done")

print("proced to Save (shadow)")
save.image("ShA_MARINa_100.RData")


print("proced to compute Synergy analysis")

# Synergy analysis
ShA_mrs_100 <- marinaCombinatorial(ShA_mrs_100, regulators=25)

######### F I N A L  S A V E  ##########

print("proced to Save")
rm(list=to_clean,to_clean) #cleaning
save.image("ShA_MARINa_100.RData")

print("All done!")