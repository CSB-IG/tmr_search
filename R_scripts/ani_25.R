############ ani 50 ###########

######## LOADS ########

setwd("~/MARINa")
library(ssmarina)
library(mixtools)
load("./Rdatas/general_MARINa.RData")
load("./Rdatas/ani_regulon_25.RData")

to_clean <- ls()

print("Source and data set loaded")

########  M A R I N A ######## 
print("proced to compute the master regulators")
ani_mrs_25 <- marina(signature, ani_regulon_25, nullmodel)


# Leading-edge analysis
ani_mrs_25 <- ledge(ani_mrs_25)

print("master regulators computed")

print(summary(ani_mrs_25))

ani_top_100_25 <-summary(ani_mrs_25,mrs=100)
write.table(ani_top_100_25,file = "lists/ani_top_100mrs_25.txt",quote = FALSE, sep = "\t")

print("Top 100 Master Regulators saved in 'ani_top_100mrs_25.txt' file")

print("proced to Save (mrs)")
save.image("ani_MARINa_25.RData")

#######  S H A D O W and S Y N E R G Y A N A L Y S I S #########
print("proced to run Shadow and Synergy analysis")

print("proced to compute the Shadow analysis")
# Shadow analysis
ani_mrshadow_25 <- shadow(ani_mrs_25, pval=25)

print("Shadow analysis done")

print("proced to Save (shadow)")
save.image("ani_MARINa_25.RData")


print("proced to compute Synergy analysis")

# Synergy analysis
ani_mrs_25 <- marinaCombinatorial(ani_mrs_25, regulators=15)

######### F I N A L  S A V E  ##########

rm(list=to_clean,to_clean) #cleaning
print("celaned")
print("Saving...")
save.image("ani_MARINa_25.RData")

print("All done!")