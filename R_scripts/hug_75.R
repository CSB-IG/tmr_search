############ hug 50 ###########

######## LOADS ########

setwd("~/MARINa")
library(ssmarina)
library(mixtools)
load("./Rdatas/general_MARINa.RData")
load("./Rdatas/hugo_regulon_75.RData")

to_clean <- ls()

print("Source and data set loaded")

########  M A R I N A ######## 
print("proced to compute the master regulators")
hugo_mrs_75 <- marina(signature, hugo_regulon_75, nullmodel)


# Leading-edge analysis
hugo_mrs_75 <- ledge(hugo_mrs_75)

print("master regulators computed")

print(summary(hugo_mrs_75))

hugo_top_100_75 <-summary(hugo_mrs_75,mrs=100)
write.table(hugo_top_100_75,file = "./lists/hugo_top_100mrs_75.txt",quote = FALSE, sep = "\t")

print("Top 100 Master Regulators saved in 'hugo_top_100mrs_75.txt' file")

print("proced to Save (mrs)")
save.image("hugo_MARINa_75.RData")

#######  S H A D O W and S Y N E R G Y A N A L Y S I S #########
print("proced to run Shadow and Synergy analysis")

print("proced to compute the Shadow analysis")
# Shadow analysis
hugo_mrshadow_75 <- shadow(hugo_mrs_75, pval=25)

print("Shadow analysis done")

print("proced to Save (shadow)")
save.image("hugo_MARINa_75.RData")


print("proced to compute Synergy analysis")

# Synergy analysis
hugo_mrs_75 <- marinaCombinatorial(hugo_mrs_75, regulators=25)

######### F I N A L  S A V E  ##########

print("proced to Save")
rm(list=to_clean,to_clean) #cleaning
save.image("hugo_MARINa_75.RData")

print("All done!")