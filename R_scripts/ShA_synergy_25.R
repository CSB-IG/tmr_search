# Synergy analysis
setwd("~/MARINa")
library(ssmarina)
load("ShA_MARINa_25.RData")
print("Source and data set loaded")
print("proced to compute Synergy analysis")

ShA_mrs_25 <- marinaCombinatorial(ShA_mrs_25, regulators=15)

######### F I N A L  S A V E  ##########

print("proced to Save")
rm(list=to_clean,to_clean) #cleaning
save.image("ShA_MARINa_25.RData")

print("All done!")