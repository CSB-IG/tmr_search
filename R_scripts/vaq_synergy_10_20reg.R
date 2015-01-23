# Synergy analysis
setwd("~/MARINa")
library(ssmarina)
load("vaq_MARINa_10.RData")
print("Source and data set loaded")
print("proced to compute Synergy analysis")

vaq_mrs_10 <- marinaCombinatorial(vaq_mrs_10, regulators=20)

######### F I N A L  S A V E  ##########

print("proced to Save")
rm(list=to_clean,to_clean) #cleaning
save.image("vaq_MARINa_25_20reg.RData")

print("All done!")
