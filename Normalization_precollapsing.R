
# open librarys

library(affy) 
library(frma)
library(sva)
library(annotate)
library(hgu133a.db)
library(limma)


# Read CEL files from the Working Directory
Data = ReadAffy()

# Parameters for graphics

N=length(Data@phenoData@data$sample)
pm.mm=0
for (i in 1:N) {pm.mm[i] = mean(mm(Data[,i])>pm(Data[,i]))}
mycolors = rep(c("blue","red","green", "magenta"), each = 2)

# Control plot for raw data

pdf("DataGraphs.pdf",width=7,height=5)
hist(Data, col=mycolors, main="Raw data distribution")
boxplot(Data,col=mycolors, main="Raw data distribution")
plot(100*pm.mm, type='h', main='Percent of MMs > PMs', ylab="%",xlab="Microarrays", ylim=c(0,50), col="red", lwd=5 )
grid(nx = NULL, ny = 6, col = "blue", lty = "dotted",lwd = par("lwd"), equilogs = TRUE)
dev.off()

# Frozen Robust Multi-Array summarization/normalization
frmaData <- frma(Data, summarize="robust_weighted_average") #tardado

# Extracting the expresion matrix
edata<-exprs(frmaData)

# Control plot for summarized data

pdf("frmaNormalized.pdf",width=7,height=5)
mycolors = rep(c("blue","red","green", "magenta"), each = 2)
plotDensity(edata, col=mycolors, main="frma normalization")
boxplot(edata,col=mycolors, main="Normalized data distribution")
dev.off()

################# BATCH FILE #########################


GSMnames <- colnames(edata)

GSMnames <- data.frame(lapply(GSMnames, function(v) {
  if (is.character(v)) return(toupper(v))
  else return(v)
}))

GSMnames <- as.character(GSMnames)

colnames(edata) <- GSMnames

batch <- c((rep(0,N)))

# Enfermos

GSE1456<-read.table("/lists_of_names/all_GSE1456.txt", colClasses = "character") # This .txt files contain the list of GSM ID for any GSE
GSE1561<-read.table("/lists_of_names/all_GSE1561.txt", colClasses = "character") # This .txt files contain the list of GSM ID for any GSE
GSE2603<-read.table("/lists_of_names/all_GSE2603.txt", colClasses = "character") # This .txt files contain the list of GSM ID for any GSE
GSE2990<-read.table("/lists_of_names/all_GSE2990.txt", colClasses = "character") # This .txt files contain the list of GSM ID for any GSE
GSE3494<-read.table("/lists_of_names/all_GSE3494.txt", colClasses = "character") # This .txt files contain the list of GSM ID for any GSE
GSE4922<-read.table("/lists_of_names/all_GSE4922.txt", colClasses = "character") # This .txt files contain the list of GSM ID for any GSE
GSE7390<-read.table("/lists_of_names/all_GSE7390.txt", colClasses = "character") # This .txt files contain the list of GSM ID for any GSE

GSE1456 <- GSE1456[[1]]
GSE1561 <- GSE1561[[1]]
GSE2603 <- GSE2603[[1]]
GSE2990 <- GSE2990[[1]]
GSE3494 <- GSE3494[[1]]
GSE4922 <- GSE4922[[1]]
GSE7390 <- GSE7390[[1]]

n1456 <- which(GSMnames %in% GSE1456)
n1561 <- which(GSMnames %in% GSE1561)
n2603 <- which(GSMnames %in% GSE2603)
n2990 <- which(GSMnames %in% GSE2990)
n3494 <- which(GSMnames %in% GSE3494)
n4922 <- which(GSMnames %in% GSE4922)
n7390 <- which(GSMnames %in% GSE7390)

batch[n1456] = 1
batch[n1561] = 2
batch[n2603] = 3
batch[n2990] = 4
batch[n3494] = 5
batch[n4922] = 6
batch[n7390] = 7

case <- c(n1456, n1561, n2603, n2990, n3494, n4922, n7390)

# Sanos

GSE15852 <-read.table("/lists_of_names/all_GSE15852.txt", colClasses = "character") # This .txt files contain the list of GSM ID for any GSE
GSE6883 <-read.table("/lists_of_names/all_GSE6883.txt", colClasses = "character") # This .txt files contain the list of GSM ID for any GSE
GSE9574 <-read.table("/lists_of_names/all_GSE9574.txt", colClasses = "character") # This .txt files contain the list of GSM ID for any GSE

GSE15852 <- GSE15852[[1]]
GSE6883 <- GSE6883 [[1]]
GSE9574 <- GSE9574 [[1]]

n15852 <- which(GSMnames %in% GSE15852)
n6883  <- which(GSMnames %in% GSE6883)
n9574  <- which(GSMnames %in% GSE9574)

batch[n15852] = 8
batch[n6883] = 9
batch[n9574] = 10

control <- c(n15852, n6883, n9574)

##  ComBat for separated treatments ##

# split of the expression matrix
caseExp <- edata[,-control]
healtExp <- edata[,-case]

# split of the batch file
caseBatch <- c(batch[-control])
healtBatch <- c(batch[-case])

# ComBat for bouth groups
case_combat = ComBat(dat=caseExp, batch=caseBatch, mod=NULL)
healt_combat = ComBat(dat=healtExp, batch=healtBatch, mod=NULL)

# Join of the tu matrix
combat_2ways <- matrix(rep(0,(22283*N)), ncol=N)
colnames(combat_2ways) <- colnames(edata)
rownames(combat_2ways) <- rownames(edata)
combat_2ways[,control] <- healt_combat
combat_2ways[,case] <- case_combat

# Control plot for jioned matrix

pdf("pre2waysCombat_robust_weighted_average.pdf",width=7,height=5)
mycolors = rep(c("blue","red","green", "magenta"), each = 2)
plotDensity(combat_2ways, col=mycolors, main="Pair combat normalization", sub="background=none summarize=random_effect")
boxplot(combat_2ways,col=mycolors, main="Normalized data distribution")
dev.off()

# Normalization of both Combat normalized subsets
n_combat_2ways <- normalizeBetweenArrays(combat_2ways, method="cyclicloess")



condition <- c(1:880)
condition[case]=1
condition[control]=2
nsample <- N
sample <- 1:nsample
pdata <- data.frame(sample, batch, condition)
modmatrix = model.matrix(~as.factor(condition), data=pdata)

batchqc_heatmap(n_combat_2ways, batch=batch, mod=modmatrix)

batchQC(edata, batch=batch, mod=modmatrix, report_file="batchqc_report_raw.html", report_dir=".")


# Control plot for jioned cyclic loess normalized matrix

pdf("2waysCombat_robust_weighted_average_cyclicloess.pdf",width=7,height=5)
mycolors = rep(c("blue","red","green", "magenta"), each = 2)
plotDensity(n_combat_2ways, col=mycolors, main="Pair combat normalization", sub="background=none summarize=random_effect")
boxplot(n_combat_2ways,col=mycolors, main="Normalized data distribution")
dev.off()

### Pre-collapse preparation  ###

# Design and contingence matrix
design = matrix(rep(0,(N*2)), nrow=N)
colnames(design) = c('case','healt')
rownames(design) = colnames(combat_edata)
design[n1456,1]=1
design[n1561,1]=1
design[n2603,1]=1
design[n2990,1]=1
design[n3494,1]=1
design[n4922,1]=1
design[n7390,1]=1
design[n15852,2]=1
design[n6883,2]=1
design[n9574,2]=1
cont.matrix = makeContrasts('case - healt', levels=design)

# limma for Differential Expression's analysis

fit = lmFit(n_combat_2ways, design)
fit2 = contrasts.fit(fit, cont.matrix)
fit2 = eBayes(fit2)

n_combat_2ways_statistics <- topTable(fit2n, coef=1, adjust='fdr', n = length(row.names(n_combat_2ways)), sort = "none")

# Control plot for log Fold Change

pdf(file="2ways_combat_cyclicloess_logFC_statisitics.pdf")
boxplot(n_combat_2ways_statistics[,1])
dev.off()

# Control plot for B statistic

pdf(file="2ways_combat_cyclicloess_B_statisitics.pdf")
boxplot(n_combat_2ways_statistics[,6])
dev.off()

# Read my oun annotation of the chip 
my_annotation<-read.table(file="my_annotation_hgu133A.txt", header = TRUE, row.names=1,colClasses = "character")

# Expresion matrix whith B statistic for collapse
precolaps_2ways_combat <- cbind(my_annotation, n_combat_2ways_statistics[,6],n_combat_2ways)
colnames(precolaps_2ways_combat)[2] <- c("b")
write.table(precolaps_2ways_combat, file="exp_matrix_normalized_precolaps.txt", quote = FALSE, sep = "\t", row.names = FALSE)


                 ########################################################################################
                 #															#
                 #                                                                    External processing							#
                 #															#
                 #	Python script collapser: "microarray_colapser.py"								#
                 #	input:  exp_matrix_normalized_precolaps.txt	                         						#
                 #	output: exp_matrix_normalized_precolaps_collapsed.txt							#
                 #	$ cut --complement -f2 exp_matrix_normalized_precolaps_colapsed.txt > exp_collapsed.txt 		#		                           #
                 #															#
                 ########################################################################################



eset <- read.table("exp_collapsed.txt",header=TRUE, sep="\t", row.names=1)


save(eset, case, control, file="exp_set.RData")
