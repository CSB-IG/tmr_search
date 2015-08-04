#Programmer:  Pierre R. Bushel
#Location: NIEHS
#email  Bushel@niehs.nih.gov
#Code: R
#Program name: pvca.R
#Date: May 26, 2009

# Modified for two treatments without time series 
# and lme4 v 1.1-8
# Hugo Tovar 
# https://github.com/hachepunto 

########## load libraries ##########
library(lme4) 
####### Edit these variables according to user defined parameters and the path to your data and data files names ##########

# myPath <- "../../Data/"
# theGene_expression_file <- "ge_data_transformed_mvi_tab_delimited.TXT"
# theExperiment_data_file <- "expinfo_tab_delimited2.TXT"

load("exp_set.RData")
setwd( "~/working/directory/")

pct_threshold = 0.8 # Amount of variability desired to be explained by the principal components.  Set to match the results in book chapter and SAS code.  User can adjust this to a higher (>= 0.8) number but < 1.0

### In addition, be sure to modify the mixed linear model by adding the appropriate random effects terms in the model

n_sample = 880

################################################

# theGEDFilePath = paste(myPath,theGene_expression_file, sep="")
# theExpDataFilePath = paste(myPath,theExperiment_data_file, sep="")

########## Load data ##########

theDataMatrix <- eset
dataRowN <- nrow(theDataMatrix)
dataColN <- ncol(theDataMatrix)

########## Center the data (center rows) ##########
theDataMatrixCentered <- matrix(data = 0, nrow = dataRowN, ncol = dataColN)
theDataMatrixCentered_transposed = apply(theDataMatrix, 1, scale, center = TRUE, scale = FALSE)
theDataMatrixCentered = t(theDataMatrixCentered_transposed)


exp_design = as.data.frame(matrix(rep(0,(n_sample*4),nrow=n_sample))
exp_design[,1]<-c(1:n_sample)
rownames(exp_design)<-exp_design[,1]
colnames(exp_design) <-c("sample","Treatment","Batch","columnname")
exp_design[,2]<-tratamiento
exp_design[,3]<-batch
exp_design[,4]<-GSMnames

expDesignRowN <- nrow(exp_design)
expDesignColN <- ncol(exp_design)
myColNames <- names(exp_design)


########## Compute correlation matrix ##########

theDataCor <- cor(theDataMatrixCentered)

########## Obtain eigenvalues ##########

eigenData <- eigen(theDataCor)
eigenValues = eigenData$values
ev_n <- length(eigenValues)
eigenVectorsMatrix = eigenData$vectors
eigenValuesSum = sum(eigenValues)
percents_PCs = eigenValues /eigenValuesSum

########## Merge experimental file and eigenvectors for n components ##########

my_counter_2 = 0
my_sum_2 = 1
for (i in ev_n:1){
my_sum_2  = my_sum_2 - percents_PCs[i]
	if ((my_sum_2) <= pct_threshold ){
		my_counter_2 = my_counter_2 + 1
	}

}
if (my_counter_2 < 3){
	pc_n  = 3

}else {
	pc_n = my_counter_2 
}

# pc_n is the number of principal components to model

pc_data_matrix <- matrix(data = 0, nrow = (expDesignRowN*pc_n), ncol = 1)
mycounter = 0
for (i in 1:pc_n){
	for (j in 1:expDesignRowN){
	mycounter <- mycounter + 1
		pc_data_matrix[mycounter,1] = eigenVectorsMatrix[j,i]

	}
}

AAA <- exp_design[rep(1:expDesignRowN,pc_n),]

Data <- cbind(AAA,pc_data_matrix)

####### Edit these variables according to your factors #######

Data$Treatment <- as.factor(Data$Treatment)
Data$Batch <- as.factor(Data$Batch)

########## Mixed linear model ##########
op <- options(warn = (-1)) 
effects_n = (expDesignColN - 2) + ((expDesignColN - 2)*(((expDesignColN - 2)-1)))/2 + 1
randomEffectsMatrix <- matrix(data = 0, nrow = pc_n, ncol = effects_n)


for (i in 1:pc_n){
	y = (((i-1)*expDesignRowN)+1)
	#randomEffects <- (summary(Rm1ML <- lmer(pc_data_matrix ~ (1|Time) + (1|Treatment) + (1|Batch) + (1|Time:Treatment) + (1|Time:Batch) + (1|Treatment:Batch), Data[y:(((i-1)*expDesignRowN)+expDesignRowN),], REML = TRUE, control=lmerControl(maxIter = 1000000, msMaxIter=1000000, singular.ok=TRUE,tolerance=1e-4, returnObject=TRUE),verbose = FALSE, na.action = na.omit))@REmat)
	randomEffects <- as.data.frame(summary(Rm1ML <- lmer(pc_data_matrix ~ (1|Treatment) + (1|Batch) + (1|Treatment:Batch), Data[y:(((i-1)*expDesignRowN)+expDesignRowN),], REML = TRUE, verbose = FALSE, na.action = na.omit))$varcor)
	randomEffectsMatrix[i,] = randomEffects[,4]
	
}

effectsNames <- randomEffects[,1]

########## Standardize Variance ##########

randomEffectsMatrixStdze <- matrix(data = 0, nrow = pc_n, ncol = effects_n)
for (i in 1:pc_n){
	mySum = sum(randomEffectsMatrix[i,])
	for (j in 1:effects_n){
		randomEffectsMatrixStdze[i,j] = randomEffectsMatrix[i,j]/mySum	
	}
}

########## Compute Weighted Proportions ##########

randomEffectsMatrixWtProp <- matrix(data = 0, nrow = pc_n, ncol = effects_n)
for (i in 1:pc_n){
	weight = eigenValues[i]/eigenValuesSum
	for (j in 1:effects_n){
		randomEffectsMatrixWtProp[i,j] = randomEffectsMatrixStdze[i,j]*weight
	}
}

########## Compute Weighted Ave Proportions ##########

randomEffectsSums <- matrix(data = 0, nrow = 1, ncol = effects_n)
randomEffectsSums <-colSums(randomEffectsMatrixWtProp)
totalSum = sum(randomEffectsSums)
randomEffectsMatrixWtAveProp <- matrix(data = 0, nrow = 1, ncol = effects_n)

for (j in 1:effects_n){
	randomEffectsMatrixWtAveProp[j] = randomEffectsSums[j]/totalSum 	
	
}

pdf("pvca.pdf",width=5,height=7)
bp <- barplot(randomEffectsMatrixWtAveProp,  xlab = "Effects", ylab = "Weighted average proportion variance", ylim= c(0,1.1),col = c("blue"), las=2)

axis(1, at = bp, labels = effectsNames, xlab = "Effects", cex.axis = 0.5, las=2)
values = randomEffectsMatrixWtAveProp
new_values = round(values , 3)
text(bp,randomEffectsMatrixWtAveProp,labels = new_values, pos=3, cex = 0.8) # place numbers on top of bars 
dev.off()
