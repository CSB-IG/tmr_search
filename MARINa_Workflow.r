#########################################################################################
################################## M  A  R  I  N  a #####################################
#########################################################################################

######## Loads ########

loadhistory(file="MARINa.Rhistory")
load("MARINa.RData")

library(ssmarina)
library(mixtools)

#########################################################################################
#############################   Vignette example   ######################################

### Basic ###
library(bcellExample)
data(bcellExample)
adjfile <- file.path(find.package("bcellExample"), "aracne", "bcellaracne.adj")
regul <- aracne2regulon_debug(adjfile, dset)
naiveBcell <- which(colnames(dset) %in% normalSamples[["N"]])
GCBcell <- which(colnames(dset) %in% c(normalSamples[["CB"]], normalSamples[["CC"]]))
signature <- rowTtest(dset[, GCBcell], dset[, naiveBcell])
signature <- (qnorm(signature$p.value/2, lower.tail=F) * sign(signature$statistic))[, 1]
nullmodel <- ttestNull(dset[, GCBcell], dset[, naiveBcell], per=1000, repos=T)
mrs <- marina(signature, regulon, nullmodel)
summary(mrs)
plot(mrs, cex=.7)
mrs <- ledge(mrs)
summary(mrs)

### Beyond ###
# Bootstraped
signature <- bootstrapTtest(dset[, GCBcell], dset[, naiveBcell], per=100)
mrs <- marina(signature, regulon, nullmodel)
mrs <- bootstrapMarina(mrs, "mode")
plot(mrs, cex=.7)
# Shadow
mrshadow <- shadow(mrs, pval=25)
summary(mrshadow)
# Synergy
mrs <- marinaCombinatorial(mrs, regulators=25)
mrs <- marinaSynergy(mrs)
summary(mrs)

##########################################################################################



#########################################################################################
############################ M A R I N a  W o r k f l o w ###############################
#########################################################################################



# generate the regulon object from the ARACNe output data
adjfile <- file.path("/home/hachepunto/MARINa/transfac_10.adj")

# and the expression data
dset <- read.table("/home/hachepunto/MARINa/Expression_Matrix_genesym_todos_sin_duplicados_colapsed.txt",header=TRUE,row.names=1)
dset = as.matrix(dset)

# Regulon is generate

# W A R N I N G ->
#regulon <- aracne2regulon(adjfile, dset)
#regulon_ShA <- aracne2regulon(adjfile_ShA, dset) #Para la lista de Califano
# <- W A R N I N G

# This "Warning" was solved by increasing the number of iterations defoult (1000 to 5000 while the necessary 
# number round 2500 iterations). Modifications and details on file /Users/hachepunto/MARINa/aracne_debug.r

regulon_iter <- aracne2regulon_debug(adjfile_ShA, dset)


# List of cases
sanos <- c("GSM158644.CEL", "GSM158645.CEL", "GSM158646.CEL","GSM241999.CEL","GSM242000.CEL","GSM242001.CEL","GSM242002.CEL","GSM242003.CEL","GSM242004.CEL","GSM242005.CEL","GSM242006.CEL","GSM242007.CEL","GSM242008.CEL","GSM242009.CEL","GSM242010.CEL","GSM242011.CEL","GSM242012.CEL","GSM242013.CEL","GSM398074.CEL","GSM398076.CEL","GSM398078.CEL","GSM398080.CEL","GSM398082.CEL","GSM398084.CEL","GSM398086.CEL","GSM398088.CEL","GSM398090.CEL","GSM398092.CEL","GSM398094.CEL","GSM398096.CEL","GSM398098.CEL","GSM398100.CEL","GSM398102.CEL","GSM398104.CEL","GSM398106.CEL","GSM398108.CEL","GSM398110.CEL","GSM398112.CEL","GSM398114.CEL","GSM398116.CEL","GSM398118.CEL","GSM398120.CEL","GSM398122.CEL","GSM398124.CEL","GSM398126.CEL","GSM398128.CEL","GSM398130.CEL","GSM398132.CEL","GSM398134.CEL","GSM398136.CEL","GSM398138.CEL","GSM398140.CEL","GSM398142.CEL","GSM398144.CEL","GSM398146.CEL","GSM398148.CEL","GSM398150.CEL","GSM398152.CEL","GSM398154.CEL","GSM398156.CEL","GSM398158.CEL")
enfermos <- c("GSM107072.CEL","GSM107073.CEL","GSM107074.CEL","GSM107075.CEL","GSM107076.CEL","GSM107077.CEL","GSM107078.CEL","GSM107079.CEL","GSM107080.CEL","GSM107081.CEL","GSM107082.CEL","GSM107083.CEL","GSM107084.CEL","GSM107085.CEL","GSM107086.CEL","GSM107087.CEL","GSM107088.CEL","GSM107089.CEL","GSM107090.CEL","GSM107091.CEL","GSM107092.CEL","GSM107093.CEL","GSM107094.CEL","GSM107095.CEL","GSM107096.CEL","GSM107097.CEL","GSM107098.CEL","GSM107099.CEL","GSM107100.CEL","GSM107101.CEL","GSM107102.CEL","GSM107103.CEL","GSM107104.CEL","GSM107105.CEL","GSM107106.CEL","GSM107107.CEL","GSM107108.CEL","GSM107109.CEL","GSM107110.CEL","GSM107111.CEL","GSM107112.CEL","GSM107113.CEL","GSM107114.CEL","GSM107115.CEL","GSM107117.CEL","GSM107118.CEL","GSM107119.CEL","GSM107120.CEL","GSM107121.CEL","GSM107122.CEL","GSM107123.CEL","GSM107124.CEL","GSM107125.CEL","GSM107126.CEL","GSM107127.CEL","GSM107128.CEL","GSM107129.CEL","GSM107130.CEL","GSM107131.CEL","GSM107132.CEL","GSM107133.CEL","GSM107134.CEL","GSM107135.CEL","GSM107136.CEL","GSM107137.CEL","GSM107138.CEL","GSM107139.CEL","GSM107140.CEL","GSM107141.CEL","GSM107142.CEL","GSM107143.CEL","GSM107144.CEL","GSM107145.CEL","GSM107146.CEL","GSM107147.CEL","GSM107148.CEL","GSM107149.CEL","GSM107150.CEL","GSM107151.CEL","GSM107152.CEL","GSM107153.CEL","GSM107154.CEL","GSM107155.CEL","GSM107156.CEL","GSM107157.CEL","GSM107158.CEL","GSM107159.CEL","GSM107160.CEL","GSM107161.CEL","GSM107162.CEL","GSM107163.CEL","GSM107164.CEL","GSM107165.CEL","GSM107166.CEL","GSM107167.CEL","GSM107168.CEL","GSM107169.CEL","GSM107170.CEL","GSM107171.CEL","GSM107172.CEL","GSM107173.CEL","GSM107174.CEL","GSM107175.CEL","GSM107176.CEL","GSM107177.CEL","GSM107178.CEL","GSM107179.CEL","GSM107180.CEL","GSM107181.CEL","GSM107182.CEL","GSM107183.CEL","GSM107184.CEL","GSM107185.CEL","GSM107186.CEL","GSM107187.CEL","GSM107188.CEL","GSM107189.CEL","GSM107190.CEL","GSM107191.CEL","GSM107192.CEL","GSM107193.CEL","GSM107194.CEL","GSM107195.CEL","GSM107196.CEL","GSM107197.CEL","GSM107198.CEL","GSM107199.CEL","GSM107200.CEL","GSM107201.CEL","GSM107202.CEL","GSM107203.CEL","GSM107204.CEL","GSM107205.CEL","GSM107206.CEL","GSM107207.CEL","GSM107208.CEL","GSM107209.CEL","GSM107210.CEL","GSM107211.CEL","GSM107212.CEL","GSM107213.CEL","GSM107214.CEL","GSM107215.CEL","GSM107216.CEL","GSM107217.CEL","GSM107218.CEL","GSM107219.CEL","GSM107220.CEL","GSM107221.CEL","GSM107222.CEL","GSM107223.CEL","GSM107224.CEL","GSM107225.CEL","GSM107226.CEL","GSM107227.CEL","GSM107228.CEL","GSM107229.CEL","GSM107230.CEL","GSM107231.CEL","GSM110625.CEL","GSM110626.CEL","GSM110627.CEL","GSM110628.CEL","GSM110629.CEL","GSM110630.CEL","GSM110631.CEL","GSM110632.CEL","GSM110633.CEL","GSM110634.CEL","GSM110635.CEL","GSM110636.CEL","GSM110637.CEL","GSM110638.CEL","GSM110639.CEL","GSM110640.CEL","GSM110641.CEL","GSM110642.CEL","GSM110643.CEL","GSM110644.CEL","GSM110645.CEL","GSM110646.CEL","GSM110647.CEL","GSM110648.CEL","GSM110649.CEL","GSM110650.CEL","GSM110651.CEL","GSM110652.CEL","GSM110653.CEL","GSM110654.CEL","GSM110655.CEL","GSM110656.CEL","GSM110657.CEL","GSM110658.CEL","GSM110659.CEL","GSM110660.CEL","GSM110661.CEL","GSM110662.CEL","GSM110663.CEL","GSM110664.CEL","GSM110665.CEL","GSM110666.CEL","GSM110667.CEL","GSM110668.CEL","GSM110669.CEL","GSM110670.CEL","GSM110671.CEL","GSM110672.CEL","GSM110673.CEL","GSM110674.CEL","GSM110675.CEL","GSM110676.CEL","GSM110677.CEL","GSM110678.CEL","GSM110679.CEL","GSM110680.CEL","GSM110681.CEL","GSM110682.CEL","GSM110683.CEL","GSM110684.CEL","GSM110685.CEL","GSM110686.CEL","GSM110687.CEL","GSM110688.CEL","GSM110689.CEL","GSM110690.CEL","GSM110691.CEL","GSM110692.CEL","GSM110693.CEL","GSM110694.CEL","GSM110695.CEL","GSM110696.CEL","GSM110697.CEL","GSM110698.CEL","GSM110699.CEL","GSM110700.CEL","GSM110701.CEL","GSM110702.CEL","GSM110703.CEL","GSM110704.CEL","GSM110705.CEL","GSM110706.CEL","GSM110707.CEL","GSM110708.CEL","GSM110709.CEL","GSM110710.CEL","GSM110711.CEL","GSM110712.CEL","GSM110713.CEL","GSM110714.CEL","GSM110715.CEL","GSM110716.CEL","GSM110717.CEL","GSM110718.CEL","GSM110719.CEL","GSM110720.CEL","GSM110721.CEL","GSM110722.CEL","GSM110723.CEL","GSM110724.CEL","GSM110725.CEL","GSM110726.CEL","GSM110727.CEL","GSM110728.CEL","GSM110729.CEL","GSM110730.CEL","GSM110731.CEL","GSM110732.CEL","GSM110733.CEL","GSM110734.CEL","GSM110735.CEL","GSM110736.CEL","GSM110737.CEL","GSM110738.CEL","GSM110739.CEL","GSM110740.CEL","GSM110741.CEL","GSM110742.CEL","GSM110743.CEL","GSM110744.CEL","GSM110745.CEL","GSM110746.CEL","GSM110747.CEL","GSM110748.CEL","GSM110749.CEL","GSM110750.CEL","GSM110751.CEL","GSM110752.CEL","GSM110753.CEL","GSM110754.CEL","GSM110755.CEL","GSM110756.CEL","GSM110757.CEL","GSM110758.CEL","GSM110759.CEL","GSM110760.CEL","GSM110761.CEL","GSM110762.CEL","GSM110763.CEL","GSM110764.CEL","GSM110765.CEL","GSM110766.CEL","GSM110767.CEL","GSM110768.CEL","GSM110769.CEL","GSM110770.CEL","GSM110771.CEL","GSM110772.CEL","GSM110773.CEL","GSM110774.CEL","GSM110775.CEL","GSM110776.CEL","GSM110777.CEL","GSM110778.CEL","GSM110779.CEL","GSM110780.CEL","GSM110781.CEL","GSM110782.CEL","GSM110783.CEL","GSM110784.CEL","GSM110785.CEL","GSM110786.CEL","GSM110787.CEL","GSM110788.CEL","GSM110789.CEL","GSM110790.CEL","GSM110791.CEL","GSM110792.CEL","GSM110793.CEL","GSM110794.CEL","GSM110795.CEL","GSM110796.CEL","GSM110797.CEL","GSM110798.CEL","GSM110799.CEL","GSM110800.CEL","GSM110801.CEL","GSM110802.CEL","GSM110803.CEL","GSM110804.CEL","GSM110805.CEL","GSM110806.CEL","GSM110807.CEL","GSM110808.CEL","GSM110809.CEL","GSM110810.CEL","GSM110811.CEL","GSM110812.CEL","GSM110813.CEL","GSM110814.CEL","GSM110815.CEL","GSM110816.CEL","GSM110817.CEL","GSM110818.CEL","GSM110819.CEL","GSM110820.CEL","GSM110821.CEL","GSM110822.CEL","GSM110823.CEL","GSM110824.CEL","GSM110825.CEL","GSM110826.CEL","GSM110827.CEL","GSM110828.CEL","GSM110829.CEL","GSM110830.CEL","GSM110831.CEL","GSM110832.CEL","GSM110833.CEL","GSM110834.CEL","GSM110835.CEL","GSM110836.CEL","GSM110837.CEL","GSM110838.CEL","GSM110839.CEL","GSM110840.CEL","GSM110841.CEL","GSM110842.CEL","GSM110843.CEL","GSM110844.CEL","GSM110845.CEL","GSM110846.CEL","GSM110847.CEL","GSM110848.CEL","GSM110849.CEL","GSM110850.CEL","GSM110851.CEL","GSM110852.CEL","GSM110853.CEL","GSM110854.CEL","GSM110855.CEL","GSM110856.CEL","GSM110857.CEL","GSM110858.CEL","GSM110859.CEL","GSM110860.CEL","GSM110861.CEL","GSM110862.CEL","GSM110863.CEL","GSM110864.CEL","GSM110865.CEL","GSM110866.CEL","GSM110867.CEL","GSM110868.CEL","GSM110869.CEL","GSM110870.CEL","GSM110871.CEL","GSM110872.CEL","GSM110873.CEL","GSM177885.cel","GSM177886.cel","GSM177887.cel","GSM177888.cel","GSM177889.cel","GSM177890.cel","GSM177891.cel","GSM177892.cel","GSM177893.cel","GSM177894.cel","GSM177895.cel","GSM177896.cel","GSM177897.cel","GSM177898.cel","GSM177899.cel","GSM177900.cel","GSM177901.cel","GSM177902.cel","GSM177903.cel","GSM177904.cel","GSM177905.cel","GSM177906.cel","GSM177907.cel","GSM177908.cel","GSM177909.cel","GSM177910.cel","GSM177911.cel","GSM177912.cel","GSM177913.cel","GSM177914.cel","GSM177915.cel","GSM177916.cel","GSM177917.cel","GSM177918.cel","GSM177919.cel","GSM177920.cel","GSM177921.cel","GSM177922.cel","GSM177923.cel","GSM177924.cel","GSM177925.cel","GSM177926.cel","GSM177927.cel","GSM177928.cel","GSM177929.cel","GSM177930.cel","GSM177931.cel","GSM177932.cel","GSM177933.cel","GSM177934.cel","GSM177935.cel","GSM177936.cel","GSM177937.cel","GSM177938.cel","GSM177939.cel","GSM177940.cel","GSM177941.cel","GSM177942.cel","GSM177943.cel","GSM177944.cel","GSM177945.cel","GSM177946.cel","GSM177947.cel","GSM177948.cel","GSM177949.cel","GSM177950.cel","GSM177951.cel","GSM177952.cel","GSM177953.cel","GSM177954.cel","GSM177955.cel","GSM177956.cel","GSM177957.cel","GSM177958.cel","GSM177959.cel","GSM177960.cel","GSM177961.cel","GSM177962.cel","GSM177963.cel","GSM177964.cel","GSM177965.cel","GSM177966.cel","GSM177967.cel","GSM177968.cel","GSM177969.cel","GSM177970.cel","GSM177971.cel","GSM177972.cel","GSM177973.cel","GSM177974.cel","GSM177975.cel","GSM177976.cel","GSM177977.cel","GSM177978.cel","GSM177979.cel","GSM177980.cel","GSM177981.cel","GSM177982.cel","GSM177983.cel","GSM177984.cel","GSM177985.cel","GSM177986.cel","GSM177987.cel","GSM177988.cel","GSM177989.cel","GSM177990.cel","GSM177991.cel","GSM177992.cel","GSM177993.cel","GSM177994.cel","GSM177995.cel","GSM177996.cel","GSM177997.cel","GSM177998.cel","GSM177999.cel","GSM178000.cel","GSM178001.cel","GSM178002.cel","GSM178003.cel","GSM178004.cel","GSM178005.cel","GSM178006.cel","GSM178007.cel","GSM178008.cel","GSM178009.cel","GSM178010.cel","GSM178011.cel","GSM178012.cel","GSM178013.cel","GSM178014.cel","GSM178015.cel","GSM178016.cel","GSM178017.cel","GSM178018.cel","GSM178019.cel","GSM178020.cel","GSM178021.cel","GSM178022.cel","GSM178023.cel","GSM178024.cel","GSM178025.cel","GSM178026.cel","GSM178027.cel","GSM178028.cel","GSM178029.cel","GSM178030.cel","GSM178031.cel","GSM178032.cel","GSM178033.cel","GSM178034.cel","GSM178035.cel","GSM178036.cel","GSM178037.cel","GSM178038.cel","GSM178039.cel","GSM178040.cel","GSM178041.cel","GSM178042.cel","GSM178043.cel","GSM178044.cel","GSM178045.cel","GSM178046.cel","GSM178047.cel","GSM178048.cel","GSM178049.cel","GSM178050.cel","GSM178051.cel","GSM178052.cel","GSM178053.cel","GSM178054.cel","GSM178055.cel","GSM178056.cel","GSM178057.cel","GSM178058.cel","GSM178059.cel","GSM178060.cel","GSM178061.cel","GSM178062.cel","GSM178063.cel","GSM178064.cel","GSM178065.cel","GSM178066.cel","GSM178067.cel","GSM178068.cel","GSM178069.cel","GSM178070.cel","GSM178071.cel","GSM178072.cel","GSM178073.cel","GSM178074.cel","GSM178075.cel","GSM178076.cel","GSM178077.cel","GSM178078.cel","GSM178079.cel","GSM178080.cel","GSM178081.cel","GSM178082.cel","gsm26804.cel","gsm26867.cel","gsm26868.cel","gsm26869.cel","gsm26870.cel","gsm26871.cel","gsm26872.cel","gsm26873.cel","gsm26874.cel","gsm26875.cel","gsm26876.cel","gsm26877.cel","gsm26878.cel","gsm26879.cel","gsm26880.cel","gsm26881.cel","gsm26882.cel","gsm26883.cel","gsm26884.cel","gsm26885.cel","gsm26886.cel","gsm26887.cel","gsm26888.cel","gsm26889.cel","gsm26890.cel","gsm26891.cel","gsm26892.cel","gsm26893.cel","gsm26894.cel","gsm26895.cel","gsm26896.cel","gsm26897.cel","gsm26898.cel","gsm26899.cel","gsm26900.cel","gsm26901.cel","gsm26902.cel","gsm26903.cel","gsm26904.cel","gsm26905.cel","gsm26906.cel","gsm26907.cel","gsm26908.cel","gsm26909.cel","gsm26910.cel","gsm26911.cel","gsm26912.cel","gsm26913.cel","gsm26914.cel","gsm50034.cel","gsm50035.cel","gsm50036.cel","gsm50037.cel","gsm50038.cel","gsm50039.cel","gsm50040.cel","gsm50041.cel","gsm50042.cel","gsm50043.cel","gsm50044.cel","gsm50045.cel","gsm50046.cel","gsm50047.cel","gsm50048.cel","gsm50049.cel","gsm50050.cel","gsm50051.cel","gsm50052.cel","gsm50053.cel","gsm50054.cel","gsm50055.cel","gsm50056.cel","gsm50057.cel","gsm50058.cel","gsm50059.cel","gsm50060.cel","gsm50061.cel","gsm50062.cel","gsm50063.cel","gsm50064.cel","gsm50065.cel","gsm50066.cel","gsm50067.cel","gsm50068.cel","gsm50069.cel","gsm50070.cel","gsm50071.cel","gsm50072.cel","gsm50073.cel","gsm50074.cel","gsm50075.cel","gsm50076.cel","gsm50077.cel","gsm50078.cel","gsm50079.cel","gsm50080.cel","gsm50081.cel","gsm50082.cel","gsm50083.cel","gsm50084.cel","gsm50085.cel","gsm50086.cel","gsm50087.cel","gsm50088.cel","gsm50089.cel","gsm50090.cel","gsm50091.cel","gsm50092.cel","gsm50093.cel","gsm50094.cel","gsm50095.cel","gsm50096.cel","gsm50097.cel","gsm50098.cel","gsm50099.cel","gsm50100.cel","gsm50101.cel","gsm50102.cel","gsm50103.cel","gsm50104.cel","gsm50105.cel","gsm50106.cel","gsm50107.cel","gsm50108.cel","gsm50109.cel","gsm50110.cel","gsm50111.cel","gsm50112.cel","gsm50113.cel","gsm50114.cel","gsm50115.cel","gsm50116.cel","gsm50117.cel","gsm50118.cel","gsm50119.cel","gsm50120.cel","gsm50121.cel","gsm50122.cel","gsm50123.cel","gsm50124.cel","gsm50125.cel","gsm50126.cel","gsm50127.cel","gsm50128.cel","gsm50129.cel","gsm50130.cel","gsm50131.cel","gsm50132.cel","gsm65820.cel","gsm65821.cel","gsm65822.cel","gsm65823.cel","gsm65824.cel","gsm65825.cel","gsm65826.cel","gsm65827.cel","gsm65828.cel","gsm65829.cel","gsm65830.cel","gsm65831.cel","gsm65832.cel","gsm65833.cel","gsm65834.cel","gsm65835.cel","gsm65836.cel","gsm65837.cel","gsm65838.cel","gsm65839.cel","gsm65840.cel","gsm65841.cel","gsm65842.cel","gsm65843.cel","gsm65844.cel","gsm65845.cel","gsm65846.cel","gsm65847.cel","gsm65848.cel","gsm65849.cel","gsm65850.cel","gsm65851.cel","gsm65852.cel","gsm65853.cel","gsm65854.cel","gsm65855.cel","gsm65856.cel","gsm65857.cel","gsm65858.cel","gsm65859.cel","gsm65860.cel","gsm65861.cel","gsm65862.cel","gsm65863.cel","gsm65864.cel","gsm65865.cel","gsm65866.cel","gsm65867.cel","gsm65868.cel","gsm65869.cel","gsm65870.cel","gsm65871.cel","gsm65872.cel","gsm65873.cel","gsm65874.cel","gsm65875.cel","gsm65876.cel","gsm65877.cel","gsm65878.cel","gsm65879.cel","gsm65880.cel","gsm79149.cel","gsm79240.cel","gsm79284.cel","gsm79307.cel")
grupos <- list(sanos=sanos,enfermos=enfermos)

# Expresion matrix separation
sanMuest <- which(colnames(dset) %in% grupos[["sanos"]])
enfMuest <- which(colnames(dset) %in% grupos[["enfermos"]])

# Row by row T test comparation of the two matrix
signature <- rowTtest(dset[, enfMuest], dset[, sanMuest])

# z-score values for the GES
signature <- (qnorm(signature$p.value/2, lower.tail=F) * sign(signature$statistic))[, 1]

# NULL model by sample permutations
nullmodel <- ttestNull(dset[, enfMuest], dset[, sanMuest], per=1000, repos=T)




#### *MARINa* ### Without Bootstrap 
mrs <- marina(signature, regulon, nullmodel)

# Leading-edge analysis
mrs <- ledge(mrs)

#######  ** mrs2 = mrs_sin_CLIC1_FTL ** #########

# Shadow analysis
mrshadow <- shadow(mrs2, pval=25)

# Synergy analysis
mrs2 <- marinaCombinatorial(mrs2, regulators=25) # Los 25 reguladores más altos

#### *MARINa* Bootstrap

signature <- bootstrapTtest(dset[, enfMuest], dset[, sanMuest])
mrs_ShA_BT <- marina(signature, regul_ShA, nullmodel)
mrs_ShA_BT <- bootstrapMarina(mrs, "mode")
plot(mrs, cex=.7)


###################################################################################
# Names code:
# *ShA* = A list of putative transcription factors as defined by GO annotation.
# Go annotations that were chosen include transcriptional regulation, DNA binding,
# and co-TF. by Yishai Shimoni and Mariano Alvarez
# http://figshare.com/articles/TF_list_/871524
#
# *wo* = Without
# *sin* = Without
# *w* = With
#
# mrs2 = mrs_sin_CLIC1_FTL



####################################################################################
################################### D E B U G I N G ################################
####################################################################################

# W A R N I N G :
#
# Loading the dataset...
# Generating the regulon objects...WARNING! NOT CONVERGENT!
# number of iterations= 1000
# WARNING! NOT CONVERGENT!
# number of iterations= 1000
# number of iterations= 423
# mu: -0.226710225185765, 0, 0.236399216367234. sigma: 0.116390379759305, 0.266600530589323, 0.160242547263094.

# Cols with 0 on signature: CLIC1:6468, FTL:8309
# Rows erased with the command:

dset_sin_CLIC1_FTL <- dset[-c(6468,8309),]


regulon_sin_CLIC1_FTL <- aracne2regulon(adjfile_ShA, dset_sin_CLIC1_FTL)
### ***Did't work for the no convergent warning*** ###


 #### Signature without problem genes:

sanMuest_sin_CLIC1_FTL <- which(colnames(dset_sin_CLIC1_FTL) %in% grupos[["sanos"]])
enfMuest_sin_CLIC1_FTL <- which(colnames(dset_sin_CLIC1_FTL) %in% grupos[["enfermos"]])
signature_sin_CLIC1_FTL <- rowTtest(dset_sin_CLIC1_FTL[, enfMuest_sin_CLIC1_FTL], dset_sin_CLIC1_FTL[, sanMuest_sin_CLIC1_FTL])
signature_sin_CLIC1_FTL <- (qnorm(signature_sin_CLIC1_FTL$p.value/2, lower.tail=F) * sign(signature_sin_CLIC1_FTL$statistic))[, 1]

#### Generar el modelo nulo sin genes problema:
nullmodel_sin_CLIC1_FTL <- ttestNull(dset_sin_CLIC1_FTL[, enfMuest_sin_CLIC1_FTL], dset_sin_CLIC1_FTL[, sanMuest_sin_CLIC1_FTL], per=1000, repos=T)

#### Generar el objeto marina sin genes problema:
mrs_sin_CLIC1_FTL <- marina(signature_sin_CLIC1_FTL, regulon_sin_CLIC1_FTL, nullmodel_sin_CLIC1_FTL)

#### Respalda el objeto marina sin genes problema:
mrs_sin_CLIC1_FTL_wo_ledge <- mrs_sin_CLIC1_FTL

#### Generar el objeto marina con ledge sin genes problema:
mrs_sin_CLIC1_FTL <- ledge(mrs_sin_CLIC1_FTL)
#### *** FUNCIONÓ *** #####
# ya permite hacer el análisis de Ledge


##################################################


save.image("~/MARINa/MARINa.RData")
savehistory("~/MARINa/MARINa.Rhistory")