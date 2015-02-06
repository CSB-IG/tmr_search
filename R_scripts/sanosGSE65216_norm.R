library(affy)
library(limma,gplots)
library(GO.db)
library(annotate)
library(hgu133plus2.db)

annDB = "hgu133plus2.db"
datos = ReadAffy()

pdf("Graphs.pdf",width=7,height=5)
mycolors = rep(c("blue","red","green", "magenta"), each = 2)
hist(datos, col=mycolors, main="Datos Crudos")
boxplot(datos,col=mycolors, main="DatosCrudos")
dev.off()

eset.quantiles <- expresso(datos, bgcorrect.method="none",
                           normalize.method="quantiles" , pmcorrect.method="pmonly",
                           summary.method="medianpolish")
eset<-exprs(eset.quantiles)

pdf("Grafica_normalizada.pdf",width=7,height=5)
mycolors = rep(c("blue","red","green", "magenta"), each = 2)
plotDensity(exprs(eset.quantiles), col=mycolors, main="Despues de normalizar")
boxplot(exprs(eset.quantiles),col=mycolors, main="DistribuciondeDatosCrudos")
dev.off()

#
affys<-rownames(eset)
genesymbols<-getSYMBOL(as.character(affys), "hgu133plus2.db")
x<-cbind(genesymbols, eset)
x<-ifelse(is.na(x), as.vector(rownames(x)), as.vector(x))
#
write.table(x, file="sanosGSE65216_expresion.txt",quote = FALSE, sep = "\t")

str(x)
