############ vaq2 ###########

######## LOADS ########

setwd("~/MARINa")
library(ssmarina)
library(mixtools)
load("./Rdatas/general_MARINa.RData")
load("./Rdatas/vaquerizas_regulon.RData")

to_clean <- ls()

print("Source and data set loaded")

########  M A R I N A ######## 
print("proced to compute the master regulators")
vaq_mrs <- marina(signature, vaquerizas_regulon, nullmodel)


# Leading-edge analysis
vaq_mrs <- ledge(vaq_mrs)

print("master regulators computed")

print(summary(vaq_mrs))

vaq_top_100 <-summary(vaq_mrs,mrs=100)
write.table(vaq_top_100,file = "vaq_top_100mrs2.txt",quote = FALSE, sep = "\t")

print("Top 100 Master Regulators saved in 'vaq_top_100mrs2s.txt' file")

print("proced to Save (mrs)")
save.image("vaq_MARINa2.RData")

#######  S H A D O W and S Y N E R G Y A N A L Y S I S #########
print("proced to run Shadow and Synergy analysis")

print("proced to compute the Shadow analysis")
# Shadow analysis
vaq_mrshadow <- shadow(vaq_mrs, pval=25)

print("Shadow analysis done")

print("proced to Save (shadow)")
save.image("vaq_MARINa2.RData")


print("proced to compute Synergy analysis")

# Synergy analysis
vaq_mrs <- marinaCombinatorial(vaq_mrs, regulators=15)

######### F I N A L  S A V E  ##########

print("proced to Save")
rm(list=to_clean,to_clean) #cleaning
save.image("vaq_MARINa2.RData")

print("All done!")