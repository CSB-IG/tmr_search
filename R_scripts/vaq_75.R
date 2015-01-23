############ vaq p50###########

######## LOADS ########

setwd("~/MARINa")
library(ssmarina)
library(mixtools)
load("./Rdatas/general_MARINa.RData")
load("./Rdatas/vaquerizas_regulon_75.RData")

to_clean <- ls()

print("Source and data set loaded")

########  M A R I N A ######## 
print("proced to compute the master regulators")
vaq_mrs_75 <- marina(signature, vaquerizas_regulon_75, nullmodel)


# Leading-edge analysis
vaq_mrs_75 <- ledge(vaq_mrs_75)

print("master regulators computed")

print(summary(vaq_mrs_75))

vaq_top100_75 <-summary(vaq_mrs_75,mrs=100)
write.table(vaq_top100_75,file = "./lists/vaq_top_100mrs_75.txt",quote = FALSE, sep = "\t")

print("Top 100 Master Regulators saved in 'vaq_top100mrs_75.txt' file")

print("proced to Save (mrs)")
save.image("vaq_MARINa_75.RData")

#######  S H A D O W and S Y N E R G Y A N A L Y S I S #########
print("proced to run Shadow and Synergy analysis")

print("proced to compute the Shadow analysis")
# Shadow analysis
vaq_mrshadow_75 <- shadow(vaq_mrs_75, pval=25)

print("Shadow analysis done")

print("proced to Save (shadow)")
save.image("vaq_MARINa_75.RData")


print("proced to compute Synergy analysis")

# Synergy analysis
vaq_mrs_75 <- marinaCombinatorial(vaq_mrs_75, regulators=15)

######### F I N A L  S A V E  ##########

print("proced to Save")
rm(list=to_clean,to_clean) #cleaning
save.image("vaq_MARINa_75.RData")

print("All done!")