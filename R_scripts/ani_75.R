############ ani 50 ###########

######## LOADS ########

setwd("~/MARINa")
library(ssmarina)
library(mixtools)
load("./Rdatas/general_MARINa.RData")
load("./Rdatas/ani_regulon_75.RData")

to_clean <- ls()

print("Source and data set loaded")

########  M A R I N A ######## 
print("proced to compute the master regulators")
ani_mrs_75 <- marina(signature, ani_regulon_75, nullmodel)


# Leading-edge analysis
ani_mrs_75 <- ledge(ani_mrs_75)

print("master regulators computed")

print(summary(ani_mrs_75))

ani_top_100_75 <-summary(ani_mrs_75,mrs=100)
write.table(ani_top_100_75,file = "lists/ani_top_100mrs_75.txt",quote = FALSE, sep = "\t")

print("Top 100 Master Regulators saved in 'ani_top_100mrs_75.txt' file")

print("proced to Save (mrs)")
save.image("ani_MARINa_75.RData")

#######  S H A D O W and S Y N E R G Y A N A L Y S I S #########
print("proced to run Shadow and Synergy analysis")

print("proced to compute the Shadow analysis")
# Shadow analysis
ani_mrshadow_75 <- shadow(ani_mrs_75, pval=25)

print("Shadow analysis done")

print("proced to Save (shadow)")
save.image("ani_MARINa_75.RData")


print("proced to compute Synergy analysis")

# Synergy analysis
ani_mrs_75 <- marinaCombinatorial(ani_mrs_75, regulators=25)

######### F I N A L  S A V E  ##########

rm(list=to_clean,to_clean) #cleaning
print("celaned")
print("Saving...")
save.image("ani_MARINa_75.RData")

print("All done!")