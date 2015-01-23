############ hug 50 ###########

######## LOADS ########

setwd("~/MARINa")
library(ssmarina)
library(mixtools)
load("./Rdatas/general_MARINa.RData")
load("./Rdatas/hugo_regulon_50.RData")

to_clean <- ls()

print("Source and data set loaded")

########  M A R I N A ######## 
print("proced to compute the master regulators")
hugo_mrs_50 <- marina(signature, hugo_regulon_50, nullmodel)


# Leading-edge analysis
hugo_mrs_50 <- ledge(hugo_mrs_50)

print("master regulators computed")

print(summary(hugo_mrs_50))

hugo_top_100_50 <-summary(hugo_mrs_50,mrs=100)
write.table(hugo_top_100_50,file = "./lists/hugo_top_100mrs_50.txt",quote = FALSE, sep = "\t")

print("Top 100 Master Regulators saved in 'hugo_top_100mrs_50.txt' file")

print("proced to Save (mrs)")
save.image("hugo_MARINa_50.RData")

#######  S H A D O W and S Y N E R G Y A N A L Y S I S #########
print("proced to run Shadow and Synergy analysis")

print("proced to compute the Shadow analysis")
# Shadow analysis
hugo_mrshadow_50 <- shadow(hugo_mrs_50, pval=25)

print("Shadow analysis done")

print("proced to Save (shadow)")
save.image("hugo_MARINa_50.RData")


print("proced to compute Synergy analysis")

# Synergy analysis
hugo_mrs_50 <- marinaCombinatorial(hugo_mrs_50, regulators=25)

######### F I N A L  S A V E  ##########

print("proced to Save")
rm(list=to_clean,to_clean) #cleaning
save.image("hugo_MARINa_50.RData")

print("All done!")