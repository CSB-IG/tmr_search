# Synergy analysis
setwd("~/MARINa")
library(ssmarina)
load("ani_MARINa_75.RData")
print("Source and data set loaded")
print("proced to compute Synergy analysis")

ani_mrs_75 <- marinaCombinatorial(ani_mrs_75, regulators=20)

######### F I N A L  S A V E  ##########

print("proced to Save")
rm(list=to_clean,to_clean) #cleaning
save.image("ani_MARINa_75_20reg.RData")

print("All done!")
