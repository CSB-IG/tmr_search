############ ani 50 ###########

######## LOADS ########

setwd("~/MARINa")
library(ssmarina)
library(mixtools)
load("./Rdatas/general_MARINa.RData")
load("./Rdatas/ani_regulon_50.RData")

to_clean <- ls()

print("Source and data set loaded")

########  M A R I N A ######## 
print("proced to compute the master regulators")
ani_mrs_50 <- marina(signature, ani_regulon_50, nullmodel)


# Leading-edge analysis
ani_mrs_50 <- ledge(ani_mrs_50)

print("master regulators computed")

print(summary(ani_mrs_50))

ani_top_100_50 <-summary(ani_mrs_50,mrs=100)
write.table(ani_top_100_50,file = "lists/ani_top_100mrs_50.txt",quote = FALSE, sep = "\t")

print("Top 100 Master Regulators saved in 'ani_top_100mrs_50.txt' file")

print("proced to Save (mrs)")
save.image("ani_MARINa_50.RData")

#######  S H A D O W and S Y N E R G Y A N A L Y S I S #########
print("proced to run Shadow and Synergy analysis")

print("proced to compute the Shadow analysis")
# Shadow analysis
ani_mrshadow_50 <- shadow(ani_mrs_50, pval=25)

print("Shadow analysis done")

print("proced to Save (shadow)")
save.image("ani_MARINa_50.RData")


print("proced to compute Synergy analysis")

# Synergy analysis
ani_mrs_50 <- marinaCombinatorial(ani_mrs_50, regulators=25)

######### F I N A L  S A V E  ##########

print("proced to Save")
rm(list=to_clean,to_clean) #cleaning
save.image("ani_MARINa_50.RData")

print("All done!")