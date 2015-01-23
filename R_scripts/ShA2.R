############ ShA2 ###########

######## LOADS ########

setwd("~/MARINa")
library(ssmarina)
library(mixtools)
load("general_MARINa.RData")
load("ShA_regulon.RData")
to_clean <- ls()
print("Source and data set loaded")

########  M A R I N A ######## 
print("proced to compute the master regulators")
ShA_mrs <- marina(signature, ShA_regulon, nullmodel)


# Leading-edge analysis
ShA_mrs <- ledge(ShA_mrs)

print("master regulators computed")

print(summary(ShA_mrs))

ShA_top_100 <-summary(ShA_mrs,mrs=100)
write.table(ShA_top_100,file = "ShA_top_100mrs2.txt",quote = FALSE, sep = "\t")

print("Top 100 Master Regulators saved in 'ShA_top_100mrs2s.txt' file")

print("proced to Save (mrs)")
save.image("ShA_MARINa2.RData")

#######  S H A D O W and S Y N E R G Y A N A L Y S I S #########
print("proced to run Shadow and Synergy analysis")

print("proced to compute the Shadow analysis")
# Shadow analysis
ShA_mrshadow <- shadow(ShA_mrs, pval=10)

print("Shadow analysis done")

print("proced to Save (shadow)")
save.image("ShA_MARINa2.RData")


print("proced to compute Synergy analysis")

# Synergy analysis
ShA_mrs <- marinaCombinatorial(ShA_mrs, regulators=5)

######### F I N A L  S A V E  ##########

print("proced to Save")
rm(list=to_clean,to_clean) #cleaning
save.image("ShA_MARINa2.RData")

print("All done!")