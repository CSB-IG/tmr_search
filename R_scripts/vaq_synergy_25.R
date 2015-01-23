# Synergy analysis
setwd("~/MARINa")
library(ssmarina)
load("vaq_MARINa_25.RData")
print("Source and data set loaded")
print("proced to compute Synergy analysis")

vaq_mrs_25 <- marinaCombinatorial(vaq_mrs_25, regulators=15)

######### F I N A L  S A V E  ##########

print("proced to Save")
rm(list=to_clean,to_clean) #cleaning
save.image("vaq_MARINa_25.RData")

print("All done!")