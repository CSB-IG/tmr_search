# Synergy analysis
setwd("~/MARINa")
library(ssmarina)
load("vaq_MARINa_75.RData")
print("Source and data set loaded")
print("proced to compute Synergy analysis")

vaq_mrs_75 <- marinaCombinatorial(vaq_mrs_75, regulators=25)

######### F I N A L  S A V E  ##########

print("proced to Save")
rm(list=to_clean,to_clean) #cleaning
save.image("vaq_MARINa_25_25reg.RData")

print("All done!")
