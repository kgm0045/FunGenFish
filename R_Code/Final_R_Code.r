####  DESeq2 Script for Differential Gene Expression Analysis ############
#install.packages("Rtools")

########You will need to set your working directory to the location you have your data##################
# You can do this by using  the Session menu to set working directory To Source File Directory
#Example: setwd("C:/Users/snook/Documents/Functional_Genomics/Final_Project")

## Make sure you have an empty folder for the working directory 
###containing only the gene/transcript count matrices and the PHENO_DATA files

#### Install the DESeq2 package if you have not already
# install.packages("BiocManager")
#BiocManager::install("DESeq2")

## Load the DESeq2 library 
library(DESeq2)

##########   1.3 Input data   ##############

### Input the count data, the gene(/transcript) count matrix and labels
countdata <- as.matrix(read.csv("gene_count_matrix.csv", row.names="gene_id"))
#Check to make sure it read correctly
dim(countdata)
head(countdata)


### Input the meta data or phenotype data
# Note: Make sure PHENO_DATA is laid out in the correct order, I put my version in the Github

##  Make sure the individual names match between the count data and the metadata
coldata <-(read.table("PHENO_DATA.txt", header=TRUE, row.names=1))
#Check to make sure it read correctly
dim(coldata)
head(coldata)


#Check all sample IDs in colData are also in CountData and match their orders
#Should both be TRUE, if FALSE something is wrong
all(rownames(coldata) %in% colnames(countdata))
countdata <- countdata[, rownames(coldata)]
all(rownames(coldata) == colnames(countdata))


## Create the DESEQ dataset and define the statistical model (page 6 of the manual)
dds <- DESeqDataSetFromMatrix(countData = countdata, colData=coldata,  design = ~treatment)
#Look at it
dds


#####   Prefiltering    Manual - starting at  1.3.6 
# Here we perform a minimal pre-filtering to remove rows that have less than 20 reads mapped.
dds <- dds[ rowSums(counts(dds)) > 20, ]
#Compare with your other dds output, look at dim and make sure the row number went down 
dds

## set factors for statistical analyses
#Make sure your PHENO_DATA has the treatments labeled as control/heat_stress, or rename these to whatever you have it named. 
dds$condition <- factor(dds$treatment, levels=c("control","heat_stress"))


######     1.4 Differential expression analysis
dds <- DESeq(dds)
res <- results(dds)
#Look to ensure it worked
res


# We can order our results table by the smallest adjusted p value:
  resOrdered <- res[order(res$padj),]
  resOrdered
# We can summarize some basic tallies using the summary function the default is p<0.1.
  summary(res)
#Tells how many adjusted p-values were less than 0.1
  sum(res$padj < 0.1, na.rm=TRUE)
#If the adjusted p value will be a value other than 0.1, alpha should be set to that value
#Not really important for us 
  res05 <- results(dds, alpha=0.05)
  summary(res05)
  sum(res05$padj < 0.05, na.rm=TRUE)

  
  
  
  
###    1.5.1 MA-plot
  ##plotMA shows the log2 fold changes attributable to a given variable over the meanof normalized counts. 
  ## Points will be colored if the adjusted p value is less than 0.1. 
  ## Points which fall out of the window are plotted as open triangles pointing either up or down
  plotMA(res, main="DESeq2", ylim=c(-10,10))
  
   
##  1.5.2 Plot counts - sanity check!
  
# You can select the gene to plot by rowname or by numeric index.
#plotCounts(dds, gene="FUN_012193", intgroup="treatment")
  
# You can plot the gene with th lowest adjuated P-value
  plotCounts(dds, gene=which.min(res$padj), intgroup="treatment")
  dds
  
  ##  Write your results to a file 
  write.csv(as.data.frame(resOrdered), file="DGESeq_results.csv")  
  
  ## 2.1.2 Extracting transformed values
  rld <- rlog(dds)
  vsd <- varianceStabilizingTransformation(dds)
  head(assay(rld), 3)
  
### Heatmap of the count matrix
  #library("genefilter")
  topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 50)
  
  library("pheatmap")
  mat  <- assay(vsd)[ topVarGenes, ]
  mat  <- mat - rowMeans(mat)
  anno <- as.data.frame(colData(vsd)[, c("treatment", "type")])
  df <- as.data.frame(colData(dds)[,c("treatment","type")])
    pheatmap(mat, annotation_col = anno)
  
  
#2.2.2 Heatmap of the sample-to-sample distances
  sampleDists <- dist(t(assay(rld)))
  library("RColorBrewer")
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- paste(rld$treatment)
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  pheatmap(sampleDistMatrix,
           clustering_distance_rows=sampleDists,
           clustering_distance_cols=sampleDists,
           col=colors)
  
  
# 2.2.3 Principal component plot of the samples
  plotPCA(rld, intgroup=c("treatment"))
  

############ Preparing Data for GSEA and Cytoscape.  #############
DGEresults <- read.csv("DGESeq_results.csv", stringsAsFactors = FALSE)
summary(DGEresults)
dim(DGEresults)

## Rename first column for next steps
names(DGEresults)[1]<- "NAME" 


############################# Make ranked list for GSEA ####################

#This organizes the entries by rank
DGEresults_Rank <-  within(DGEresults, rank <- sign(log2FoldChange) * -log10(pvalue))
DGEresults_Rank 
 
#Our reference genome was already annotated so we just have to subset it to only get entries with a gene NAME

#This removes all entries containing "LOC" within the gene_id, because those genes are unnamed
DGEresults_Rank <- DGEresults_Rank[!grepl("LOC", DGEresults_Rank$NAME), ]

#This removes everything in front of the "|" symbol because the reference genome have the name twice
DGEresults_Rank$NAME <- sub("^[^|]+\\|", "", DGEresults_Rank$NAME)

#This capitalizes all the GENE names, don't know if its needed but that's what the Daphnia one had
DGEresults_Rank$NAME <- toupper(DGEresults_Rank$NAME)

#Selects the first 2 columns
DGErank = subset(DGEresults_Rank, select = c(NAME,rank) )
DGErank

#sebset the results so only Gene Name and rank
DGErank_withName <- na.omit(DGErank)
DGErank_withName
dim(DGErank_withName)

#Makes txt file with Name and Rank
write.table(as.data.frame(DGErank_withName), file="DGErank_withName.txt", quote = FALSE, row.names=FALSE, sep = "\t") 
#Makes csv file with Name and Rank
write.table(as.data.frame(DGErank_withName), file="DGErank_withName.csv", quote = FALSE, row.names=FALSE, sep = ",")


####  We also need the normalized expression DATA
nt <- normTransform(dds) # defaults to log2(x+1)
head(assay(nt))
# compare to original count data
head(countdata)
# make it a new dataframe
NormTransExp<-assay(nt)
summary(NormTransExp)
gene <-gsub("^[^-]+-", "", rownames(NormTransExp))
NormTransExpIDs  <-cbind(gene,NormTransExp)

#This removes everything in front of the "|" symbol because the reference genome has the name twice
NormTransExpIDs <- data.frame(NormTransExpIDs)
NormTransExpIDs$gene <- sub("^[^|]+\\|", "", NormTransExpIDs$gene)
head(NormTransExpIDs)

#Outputs txt file with results
write.table(as.data.frame(NormTransExpIDs), file="NormTransExpressionData.txt", row.names=FALSE, quote = FALSE, sep = "\t")  
#Outputs csv file with results
write.table(as.data.frame(NormTransExpIDs), file="NormTransExpressionData.csv", row.names=FALSE, quote = FALSE, sep = ",")
