# Synergy analysis
setwd("~/MARINa")
library(ssmarina)
load("hugo_MARINa_25.RData")
print("Source and data set loaded")
print("proced to compute Synergy analysis")

hugo_mrs_25 <- marinaCombinatorial(hugo_mrs_25, regulators=15)

######### F I N A L  S A V E  ##########

print("proced to Save")
rm(list=to_clean,to_clean) #cleaning
save.image("hugo_MARINa_25.RData")

print("All done!")