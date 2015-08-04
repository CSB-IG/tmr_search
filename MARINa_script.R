#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(mixtools))
library(ssmarina)

parser <- ArgumentParser(description="MARINa Master regulators inference from .adj file and Expresion Matrix")
parser$add_argument("--exp_set",	help="RData with the expression matrix (eset) and case (case) and control (control) groups objects")
parser$add_argument("--adj", default="none",	choices=adjfile, help="ARACNE .adj file")
parser$add_argument("--regulon", default="none", help="RData with a regulon object")
parser$add_argument("--sig_null", default="none", help="RData with a signature and nullmodel object")
parser$add_argument("--output", default="getwd()", help="output folder")
args   <- parser$parse_args()

setwd(args$output)
load(args$exp_mat)

# Regulon generation

regulon <- aracne2regulon(args$adj, eset)

save(regulon, file = "./regulon.RData")

#  Signature generation

# Row by row T test comparation of the two matrix
signature <- rowTtest(eset[, case], eset[, control])
# z-score values for the GES
signature <- (qnorm(signature$p.value/2, lower.tail=F) * sign(signature$statistic))[, 1]

# Null model by sample permutations
nullmodel <- ttestNull(eset[, case], eset[, control], per=1000, repos=T)

save(signature,nullmodel,file="./sign_null.RData")

to_clean <- ls()

########  M A R I N A ######## 

mrs <- marina(signature, regulon, nullmodel)

mrs_noledges <- mrs

# Leading-edge analysis
mrs <- ledge(mrs)

# Ploting masters regulators
pdf("Top10mrs.pdf",width=6,height=7)
plot(mrs_noledges, cex=.7)
dev.off()

print(summary(mrs))

top100 <-summary(mrs,mrs=100)
write.table(top100,file = "./top100mrs.txt",quote = FALSE, sep = "\t")


#save(mrs,top100,mrs_noledges,file="MARINa.RData")

#######  S H A D O W and S Y N E R G Y A N A L Y S I S #########

# Shadow analysis
mrshadow <- shadow(mrs, pval=25)

write.table(summary(mrshadow)$Shadow.pairs, file="shadow_pairs.sif", quote=FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

# Synergy analysis
mrs <- marinaCombinatorial(mrs, regulators=15)

mrs <- marinaSynergy(mrs)

# ploting synergy regulators
pdf("Synergy.pdf",width=11,height=7)
plot(mrs,mrs=10, cex=.7)
dev.off()

######### F I N A L  S A V E  ##########
rm(list=to_clean,to_clean) #cleaning
save(mrs, mrs_noledges, mrshadow,file="MARINa.RData")

print("All done!")