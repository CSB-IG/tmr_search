############ GSE65216 p7 ###########

######## LOADS ########

setwd("~/MARINa")
library(ssmarina)
library(mixtools)
load("./validation_analisys/GSE4611/general_GSE4611.Rdata")
load("./validation_analisys/GSE4611/GSE4611_regulon.Rdata")

to_clean <- ls()

print("Source and data set loaded")

########  M A R I N A ######## 
print("proced to compute the master regulators")
GSE4611_mrs <- marina(signature, GSE4611_regulon, nullmodel)

GSE4611_mrs_noledges <- GSE4611_mrs

# Leading-edge analysis
GSE4611_mrs <- ledge(GSE4611_mrs)

print("master regulators computed")

print(summary(GSE4611_mrs))

GSE4611_top100 <-summary(GSE4611_mrs,mrs=100)
write.table(GSE4611_top100,file = "./lists/GSE4611_top100mrs.txt",quote = FALSE, sep = "\t")

print("Top 100 Master Regulators saved in 'GSE4611_top100mrs.txt' file")

print("proced to Save (mrs)")
save.image("GSE4611_MARINa.RData")

#######  S H A D O W and S Y N E R G Y A N A L Y S I S #########
print("proced to run Shadow and Synergy analysis")

print("proced to compute the Shadow analysis")
# Shadow analysis
GSE4611_mrshadow <- shadow(GSE4611_mrs, pval=25)

print("Shadow analysis done")

print("proced to Save (shadow)")
save.image("GSE4611_MARINa.RData")


print("proced to compute Synergy analysis")

# Synergy analysis
GSE4611_mrs <- marinaCombinatorial(GSE4611_mrs, regulators=20)

######### F I N A L  S A V E  ##########

print("cleaning")
rm(list=to_clean,to_clean) #cleaning
print("proced to Save")
save.image("GSE4611_MARINa.RData")

print("All done!")