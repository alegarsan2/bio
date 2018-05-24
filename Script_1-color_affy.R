#Master in Bioinformatics 2018
#Normalization and diferential expression for Affymetrix arrays


#setwd("mypath") ### Establecer directorio de trabajo

##1. Load libraries

source("https://bioconductor.org/biocLite.R")

biocLite(c("affy","genefilter","limma"))

library("affy")
library("limma")
library("genefilter")
f

#2. Import targets.txt file
targets <- readTargets("list.txt", row.names="FileName")


#3. Import .CEL files
data <- ReadAffy(filenames=targets$FileName) 		#Import intensities from Affymetrix arrays (.CEL)
													# data is an object of AffyBatch class

#4. Normalize with RMA 
#generates object eset (class ExprSet), 
#expresso function provides intensities in log scale
eset <- expresso(data,
 				bg.correct = TRUE, 
                bgcorrect.method="rma",
                normalize = TRUE, 
                normalize.method="quantiles", 
                pmcorrect.method="pmonly", 
                summary.method="medianpolish",
                verbose = TRUE,
				 ) 




#5. Generate BOXPLOTS before and after normalization

#boxplot for raw data
boxplot(data,
		 main="Boxplot Before Normalization",
		 col = "lightgrey")
		

#boxplot for normalized data
exprseset <- as.data.frame(exprs(eset))		
boxplot(data.frame(exprseset),
		 main="Boxplot After Normalization (log scale)",
		 col = "white")


#6.Data filtering using IQR.
esetIQR <- varFilter(eset, var.func=IQR, var.cutoff=0.5, filterByQuantile=TRUE)




#######Differential expression analysis.#######
#7. Design matrix.
design<-cbind(Control=c(1,1,1,1,1,0,0,0,0,0), Nanog_RNAi=c(0,0,0,0,0,1,1,1,1,1))
rownames(design)<-targets$FileName

#8. Contrasts matrix.
cont.matrix<-makeContrasts(Nanog_RNAivsControl=Nanog_RNAi-Control,levels=design) 


#9. Obtaining differentially expressed genes (DEGs)
#Linear model and eBayes 
fit<-lmFit(esetIQR,design)  ##getting DEGs from IQR 
fit2<-contrasts.fit(fit, cont.matrix)
fit2<-eBayes(fit2)

#Table with DEGs results
toptableIQR<-topTable(fit2, number=dim(exprs(esetIQR))[1], adjust.method="BH", sort.by="p")


##10. Save results
save(toptableIQR,file="MyResults.RData")


biocLite("hgu95av2.db")
library("hgu95av2.db")

#Master in Bioinformatics 2018
#Annotating gene lists using Bioconductor




setwd("mypath") ### Working directory

###Hands-On 4
##1. Load annotation library
library("hgu95av2.db")
hgu95av2()			##info annotations available in library
hgu95av2GENENAME	##class Bimap
hgu95av2GO
hgu95av2CHRLOC

hgu95av2_dbInfo()	##more info


##2. Retrieve Gene symbol, chromosomic coordinates and KEGG IDs for probe 1003_s_at.

get("1003_s_at",hgu95av2SYMBOL)		##Gene symbol
head(toTable(hgu95av2SYMBOL))		##from Genes in Bimap to data.frame
get("1003_s_at",hgu95av2CHRLOC)		##Chromosomic coordinates
head(toTable(hgu95av2CHRLOC))		##from Chromosome to data.frame
get("1003_s_at",hgu95av2PATH)		##IDs KEGG



###Hands-On 5
#Get Gene Symbols from DEGs (FDR<0.05) obtained in Hands-On 3.
load("MyResults.RData")  
biocLite("mouse4302.db")
library("mouse4302.db")				#Array Affymetrix Mouse 430 v2
ID.fdr.005.table<-subset(toptableIQR, toptableIQR$adj.P.Val<=0.05)


probenames.fdr.005<-as.character(rownames(ID.fdr.005.table))
list.GeneSymbol.fdr.005<-mget(probenames.fdr.005, mouse4302SYMBOL)
char.GeneSymbol.fdr.005<- as.character(list.GeneSymbol.fdr.005)
toptable.annotated<-cbind(ID.fdr.005.table,char.GeneSymbol.fdr.005)


#Alternatively you may use getSYMBOL from annotate library
library("annotate")
GeneSymbol2.fdr.005<-getSYMBOL(probenames.fdr.005,"mouse4302.db")		

