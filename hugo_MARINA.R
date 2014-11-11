#########################################################################################
################################## M  A  R  I  N  a #####################################
#########################################################################################

######## LOADS ########

setwd("~/MARINa")
library(ssmarina)
library(mixtools)
load("general_MARINa.RData")
to_clean <- ls()
print("Source and data set loaded")


########  R E G U L O N ######## 
### One proces by TF list

# adj from ARACNe imput data
hugo_adjfile <- file.path("/home/hachepunto/MARINa/hugo_p1e-10.adj")
print("adj loaded")


# Regulon generation
print("proced to compute the regulon")
hugo_regulon <- aracne2regulon(hugo_adjfile, dset)

print("regulon computed")


########  M A R I N A ######## 
print("proced to estimate  the Master Regulators")
hugo_mrs <- marina(signature, hugo_regulon, nullmodel)


# Leading-edge analysis
hugo_mrs <- ledge(hugo_mrs)

print("Master Regulators estimated")

print(summary(hugo_mrs))

hugo_top_100 <-summary(hugo_mrs,mrs=100)
write.table(hugo_top_100,file = "hugo_top_100mrs.txt",quote = FALSE, sep = "\t")

print("Top 100 Master Regulators saved in 'hugo_top_100mrs.txt' file")

########  B O O T S T R A P E D   M A R I N A ######## 

print("proced to calculet  bootstraped master regulators")
hugo_BT_mrs <- marina(BTsignature, hugo_regulon, nullmodel)

hugo_BT_mrs <- bootstrapMarina(hugo_BT_mrs, "mode")

print("master regulators bootstraped calculated")
print(summary(hugo_BT_mrs))

hugo_BT_top_100 <-summary(hugo_BT_mrs,mrs=100)
write.table(hugo_BT_top_100,file = "hugo_BT_top_100mrs.txt",quote = FALSE, sep = "\t")

print("Top 100 bootstraped Master Regulators saved in 'hugo_BT_top_100mrs.txt' file")


#######  S H A D O W and S Y N E R G Y A N A L Y S I S #########
print("proced to compute the Shadow and Synergy analysis")

print("proced to compute the Shadow analysis")
# Shadow analysis
hugo_mrshadow <- shadow(hugo_mrs, pval=25)

# Bootstraped Shadow analysis
hugo_mrshadow_BT <- shadow(hugo_BT_mrs, pval=25)

print("Shadow analysis done")


print("proced to compute the Synergy analysis")

# Synergy analysis
hugo_mrs <- marinaCombinatorial(hugo_mrs, regulators=25) 

# Bootstraped Synergy analysis
hugo_BT_mrs <- marinaCombinatorial(hugo_BT_mrs, regulators=25) 

print("Synergy analysis done")



##################################################
#################### S A V E  ######################
##################################################
print("proced to Save")
rm(list=to_clean,to_clean) #cleaning
save.image("hugo_MARINa.RData")
savehistory("hugo_MARINa.Rhistory")

print("All done!")
