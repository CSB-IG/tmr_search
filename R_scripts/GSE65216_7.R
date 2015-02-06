############ GSE65216 p7 ###########

######## LOADS ########

setwd("~/MARINa")
library(ssmarina)
library(mixtools)
load("./validation_analisys/general_GSE65216.Rdata")
load("./validation_analisys/GSE65216_regulon_p7.Rdata")

to_clean <- ls()

print("Source and data set loaded")

########  M A R I N A ######## 
print("proced to compute the master regulators")
GSE65216_mrs_7 <- marina(signature, GSE65216_regulon_7, nullmodel)

GSE65216_mrs_noledges_7 <- GSE65216_mrs_7

# Leading-edge analysis
GSE65216_mrs_7 <- ledge(GSE65216_mrs_7)

print("master regulators computed")

print(summary(GSE65216_mrs_7))

GSE65216_top100_7 <-summary(GSE65216_mrs_7,mrs=100)
write.table(GSE65216_top100_7,file = "./lists/GSE65216_top100mrs_7.txt",quote = FALSE, sep = "\t")

print("Top 100 Master Regulators saved in 'GSE65216_top100mrs_7.txt' file")

print("proced to Save (mrs)")
save.image("GSE65216_MARINa_7.RData")

#######  S H A D O W and S Y N E R G Y A N A L Y S I S #########
print("proced to run Shadow and Synergy analysis")

print("proced to compute the Shadow analysis")
# Shadow analysis
GSE65216_mrshadow_7 <- shadow(GSE65216_mrs_7, pval=25)

print("Shadow analysis done")

print("proced to Save (shadow)")
save.image("GSE65216_MARINa_7.RData")


print("proced to compute Synergy analysis")

# Synergy analysis
GSE65216_mrs_7 <- marinaCombinatorial(GSE65216_mrs_7, regulators=25)

######### F I N A L  S A V E  ##########

print("cleaning")
rm(list=to_clean,to_clean) #cleaning
print("proced to Save")
save.image("GSE65216_MARINa_7.RData")

print("All done!")