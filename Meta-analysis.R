####  DESeq2 Script for Differential Gene Expression Analysis in 
# Functional Genomics BIOL: 6850
### Resources and Citations:
# Love et al 2016 DESeq2 GenomeBiology
# https://bioconductor.riken.jp/packages/2.14/bioc/vignettes/DESeq2/inst/doc/DESeq2.pdf
# http://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html


### You will need to set your working directory to the location you have your data.
# You can do this by using  the Session menu to set working directory To Source File Directory
#Clearenvironment
rm(list=ls())


#### Install the DESeq2 package if you have not already
## try http:// if https:// URLs are not supported
#if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")

#BiocManager::install("DESeq2")
## Load the DESeq2 library 
library(DESeq2)

## Use the Session menu to set working directory To Source File Directory

##########   1.3 Input data   ##############

### Input the count data, the gene(/transcript) count matrix and labels
### How you inport this will depend on what your final output was from the mapper/counter that you used.
## this works with output from PrepDE.py from Ballgown folder.
countdata_standard <- as.matrix(read.csv("gene_matrix_stringtie.csv", row.names="gene_id"))
countdata_stringtie_modified <- as.matrix(read.csv("gene_matrix_stringtie_modified.csv", row.names="gene_id"))
countdata_htseq <- as.matrix(read.csv("gene_matrix_htseq.csv", row.names="gene_id"))
countdata_htseq_modified <- as.matrix(read.csv("gene_matrix_htseq_modified.csv", row.names="gene_id"))
countdata_hisat250 <- as.matrix(read.csv("gene_matrix_hisat250.csv", row.names="gene_id"))
countdata_star100 <- as.matrix(read.csv("gene_matrix_star100.csv", row.names="gene_id"))
countdata_star149 <- as.matrix(read.csv("gene_matrix_star149.csv", row.names="gene_id"))

dim(countdata_standard)
dim(countdata_stringtie_modified)
dim(countdata_htseq)
dim(countdata_htseq_modified)
dim(countdata_hisat250)
dim(countdata_star100)
dim(countdata_star149)



head(countdata_standard)

#### IF necessary,depending on what program made the count matrix. Remove the unwanted row (the * and zero row)
#countdata<- countdata[-1,c(-1)]
# OR   Remove the unwanted column (length)
#countdata<- countdata[,c(-1)]
#dim(countdata)
#head(countdata)

### Input the meta data or phenotype data
# Note: The PHENO_DATA file contains information on each sample, e.g., sex or population. The exact way to import this depends on the format of the file.
##  Make sure the individual names match between the count data and the metadata
coldata <-(read.table("PHENO_DATA.txt", header=TRUE, row.names=1))
dim(coldata)
head(coldata)


#Check all sample IDs in colData are also in CountData and match their orders
all(rownames(coldata) %in% colnames(countdata_standard))
countdata_standard <- countdata_standard[, rownames(coldata)]
all(rownames(coldata) == colnames(countdata_standard))


## Create the DESEQ dataset and define the statistical model (page 6 of the manual)
dds_standard <- DESeqDataSetFromMatrix(countData = countdata_standard, colData=coldata,  design = ~treatment)
dds_stringtie_modified <- DESeqDataSetFromMatrix(countData = countdata_stringtie_modified, colData=coldata,  design = ~treatment)
dds_htseq <- DESeqDataSetFromMatrix(countData = countdata_htseq, colData=coldata,  design = ~treatment)
dds_htseq_modified <- DESeqDataSetFromMatrix(countData = countdata_htseq_modified, colData=coldata,  design = ~treatment)
dds_hisat250 <- DESeqDataSetFromMatrix(countData = countdata_hisat250, colData=coldata,  design = ~treatment)
dds_star100 <- DESeqDataSetFromMatrix(countData = countdata_star100, colData=coldata,  design = ~treatment)
dds_star149 <- DESeqDataSetFromMatrix(countData = countdata_star149, colData=coldata,  design = ~treatment)


#look at it
dds_standard
dds_stringtie_modified
dds_htseq
dds_htseq_modified
dds_hisat250 
dds_star100
dds_star149


#####   Prefiltering    Manual - starting at  1.3.6 
# Here we perform a minimal pre-filtering to remove rows that have less than 20 reads mapped.
## You can play around with this number to see how it affects your results!
dds_standard <- dds_standard[ rowSums(counts(dds_standard)) > 20, ]
dds_stringtie_modified <- dds_stringtie_modified[ rowSums(counts(dds_stringtie_modified)) > 20, ]
dds_htseq <- dds_htseq[ rowSums(counts(dds_htseq)) > 20, ]
dds_htseq_modified <- dds_htseq_modified[ rowSums(counts(dds_htseq_modified)) > 20, ]
dds_hisat250 <- dds_hisat250[ rowSums(counts(dds_hisat250)) > 20, ]
dds_star100 <- dds_star100[ rowSums(counts(dds_star100)) > 20, ]
dds_star149 <- dds_star149[ rowSums(counts(dds_star149)) > 20, ]

# look.  How many genes were filtered out?
dds_standard
dds_stringtie_modified
dds_htseq
dds_htseq_modified
dds_hisat250
dds_star100
dds_star149


## set factors for statistical analyses
###### Note you need to change condition to treatment (to match our design above)
#  and levels to our treatment names in the PHENO_DATA: Ad_lib is the control, Caloric_Restriction is the treatment group
# example:
#dds$condition <- factor(dds$condition, levels=c("untreated","treated"))
dds_standard$condition <- factor(dds_standard$treatment, levels=c("control","heat_stress"))
dds_stringtie_modified$condition <- factor(dds_stringtie_modified$treatment, levels=c("control","heat_stress"))
dds_htseq$condition <- factor(dds_htseq$treatment, levels=c("control","heat_stress"))
dds_htseq_modified$condition <- factor(dds_htseq_modified$treatment, levels=c("control","heat_stress"))
dds_hisat250$condition <- factor(dds_hisat250$treatment, levels=c("control","heat_stress"))
dds_star100$condition <- factor(dds_star100$treatment, levels=c("control","heat_stress"))
dds_star149$condition <- factor(dds_star149$treatment, levels=c("control","heat_stress"))

######     1.4 Differential expression analysis
### Question 2. Look at the manual - what is happening at this point?
dds_standard <- DESeq(dds_standard)
res_standard <- results(dds_standard)
res_standard

dds_stringtie_modified <- DESeq(dds_stringtie_modified)
res_stringtie_modified <- results(dds_stringtie_modified)
res_stringtie_modified

dds_htseq <- DESeq(dds_htseq)
res_htseq <- results(dds_htseq)
res_htseq

dds_htseq_modified <- DESeq(dds_htseq_modified)
res_htseq_modified <- results(dds_htseq_modified)
res_htseq_modified

dds_hisat250 <- DESeq(dds_hisat250)
res_hisat250 <- results(dds_hisat250)
res_hisat250

dds_star100 <- DESeq(dds_star100)
res_star100 <- results(dds_star100)
res_star100

dds_star149 <- DESeq(dds_star149)
res_star149 <- results(dds_star149)
res_star149

###  Question 3. What does each column mean?
# We can order our results table by the smallest adjusted p value:
resOrdered_standard <- res_standard[order(res_standard$padj),]
resOrdered_standard

resOrdered_stringtie_modified <- res_stringtie_modified[order(res_stringtie_modified$padj),]
resOrdered_stringtie_modified

resOrdered_htseq <- res_htseq[order(res_htseq$padj),]
resOrdered_htseq

resOrdered_htseq_modified <- res_htseq[order(res_htseq_modified$padj),]
resOrdered_htseq_modified

resOrdered_hisat250 <- res_hisat250[order(res_hisat250$padj),]
resOrdered_hisat250

resOrdered_star100 <- res_star100[order(res_star100$padj),]
resOrdered_star100

resOrdered_star149 <- res_star149[order(res_star149$padj),]
resOrdered_star149

# We can summarize some basic tallies using the summary function the default is p<0.1.
summary(res_standard)
summary(res_stringtie_modified)
summary(res_htseq)
summary(res_htseq_modified)
summary(res_hisat250)
summary(res_star100)
summary(res_star149)


#How many adjusted p-values were less than 0.1?
sum(res_standard$padj < 0.1, na.rm=TRUE)
sum(res_stringtie_modified$padj < 0.1, na.rm=TRUE)
sum(res_htseq$padj < 0.1, na.rm=TRUE)
sum(res_htseq_modified$padj < 0.1, na.rm=TRUE)
sum(res_hisat250$padj < 0.1, na.rm=TRUE)
sum(res_star100$padj < 0.1, na.rm=TRUE)
sum(res_star149$padj < 0.1, na.rm=TRUE)

#If the adjusted p value will be a value other than 0.1, alpha should be set to that value:
res05_standard <- results(dds_standard, alpha=0.05)
summary(res05_standard)
sum(res05_standard$padj < 0.05, na.rm=TRUE)

res05_stringtie_modified <- results(dds_stringtie_modified, alpha=0.05)
summary(res05_stringtie_modified)
sum(res05_stringtie_modified$padj < 0.05, na.rm=TRUE)

res05_htseq <- results(dds_htseq, alpha=0.05)
summary(res05_htseq)
sum(res05_htseq$padj < 0.05, na.rm=TRUE)

res05_htseq_modified <- results(dds_htseq_modified, alpha=0.05)
summary(res05_htseq_modified)
sum(res05_htseq_modified$padj < 0.05, na.rm=TRUE)

res05_hisat250 <- results(dds_hisat250, alpha=0.05)
summary(res05_hisat250)
sum(res05_hisat250$padj < 0.05, na.rm=TRUE)

res05_star100 <- results(dds_star100, alpha=0.05)
summary(res05_star100)
sum(res05_star100$padj < 0.05, na.rm=TRUE)

res05_star149 <- results(dds_star149, alpha=0.05)
summary(res05_star149)
sum(res05_star149$padj < 0.05, na.rm=TRUE)

#Save_file
write.csv(as.data.frame(resOrdered_standard), file="DGESeq_results_standard.csv") 
write.csv(as.data.frame(resOrdered_stringtie_modified), file="DGESeq_results_stringtie_modified.csv") 
write.csv(as.data.frame(resOrdered_htseq), file="DGESeq_results_htseq.csv")
write.csv(as.data.frame(resOrdered_htseq_modified), file="DGESeq_results_htseq_modified.csv")
write.csv(as.data.frame(resOrdered_hisat250), file="DGESeq_results_hisat250.csv")
write.csv(as.data.frame(resOrdered_star100), file="DGESeq_results_star100.csv")
write.csv(as.data.frame(resOrdered_star149), file="DGESeq_results_star149.csv")

#Start the Meta-Analysis
rm(list=ls())
#Merge tables results for each dataset - keep the genes in common ~
library(plyr)
attach(DGESeq_results_standard)
attach(DGESeq_results_stringtie_modified)
attach(DGESeq_results_htseq)
attach(DGESeq_results_htseq_modified)
attach(DGESeq_results_hisat250)
attach(DGESeq_results_star100)
attach(DGESeq_results_star149)
library(plyr)
attach(EDGER0)
attach(EDGER10)
attach(EDGER100)
attach(EDGER20)
attach(EDGER50)
attach(DGESeq_results_0)
attach(DGESeq_results_10)
attach(DGESeq_results_100)
attach(DGESeq_results_50)

# Check if "Gene" column exists in all additional data frames
sapply(additional_dfs, function(df) "Gene" %in% colnames(df))

# Check if "Gene" column uniquely identifies each row in all additional data frames
sapply(additional_dfs, function(df) length(unique(df$Gene)) == nrow(df))


# Merge DGESeq_results_standard and DGESeq_results_stringtie_modified
merged_df <- merge(DGESeq_results_standard, DGESeq_results_stringtie_modified, by = "Gene", all.x = FALSE, suffixes = c(".Std", ".StrMo"))

# Merge additional data frames
additional_dfs <- list(DGESeq_results_htseq, DGESeq_results_htseq_modified)
for (df in additional_dfs) {
  merged_df <- merge(merged_df, df, by = "Gene", all.x = FALSE, suffixes = c(".ht",".htMo",".hi250",".st100",".st149"))
}

additional_dfs<- list(DGESeq_results_hisat250, DGESeq_results_star100)
for (df in additional_dfs) {
  merged_df <- merge(merged_df, df, by = "Gene", all.x = FALSE, suffixes = c(".hi250",".st100"))
}

additional_dfs<- list(DGESeq_results_star149)
for (df in additional_dfs) {
  merged_df <- merge(merged_df, df, by = "Gene", all.x = FALSE, suffixes = c(".st149"))
}

dim(merged_df)
head(merged_df)
summary(merged_df)

#Put p-value and Fold changes results in lists
H_rawpval <- list("pvalue_Std"= merged_df[["pvalue.Std"]], "pvalueStrMo"=merged_df[["pvalue.StrMo"]],"pvalue_ht"=merged_df[["pvalue.ht"]],"pvalue_htMo"=merged_df[["pvalue.htMo"]],"pvalue_hi250"=merged_df[["pvalue.hi250"]],"pvalue_st100"=merged_df[["pvalue.st100"]],"pvalue_st149"=merged_df[["pvalue.st149"]])
summary(H_rawpval)

H_adjpval <- list("padj_Std"=merged_df[["padj.Std"]], "padj_StrMo"=merged_df[["padj.StrMo"]],"padj_ht"=merged_df[["padj.ht"]],"padj_htMo"=merged_df[["padj.htMo"]],"padj_hi250"=merged_df[["padj.hi250"]],"padj_st100"=merged_df[["padj.st100"]],"padj_st149"=merged_df[["padj.st149"]])
summary(H_adjpval)

H_FC <- list("FC_Std"=merged_df[["log2FoldChange.Std"]], "FC"=merged_df[["log2FoldChange.StrMo"]],"FC_ht"=merged_df[["log2FoldChange.ht"]],"FC_htMo"=merged_df[["log2FoldChange.htMo"]],"FC_hi250"=merged_df[["log2FoldChange.hi250"]],"FC_st100"=merged_df[["log2FoldChange.st100"]],"FC_st149"=merged_df[["log2FoldChange.st149"]])
summary(H_FC)

#make matrix: 1 for genes DE at adjpval of 0.05, and 0 for not.
DE<- mapply(H_adjpval, FUN=function(x) ifelse(x <= 0.05, 1, 0))
studies <-c("Standard", "Stringtie_Modified","Htseq","Htseq_modified","Hisat2_50","Star_100","Star_149")
colnames(DE)=paste("DE", studies,sep=".")
rownames(DE)=paste(merged_df$"Gene")
#DE<-cbind(top.table_H_intersect$"Gene", DE)
summary(DE)
head(DE)
dim(DE)


#Calculate the meta-analysis p-value
library(metaRNASeq)
# P-value combination using Fisher method 
fishcomb <- fishercomb(H_rawpval, BHth = 0.05)
summary(fishcomb)

hist(fishcomb$rawpval, breaks = 100, col = "grey", main = "Fisher Method", xlab = "Raw p-values (meta-analysis)")

# P-value combination using inverse normal method method 
invnormcomb <- invnorm(H_rawpval, nrep=c(5,5,5,5,5,5,5), BHth=0.05)
summary(invnormcomb)
head(invnormcomb)

hist(invnormcomb$rawpval, breaks=100, col="grey", main="Inverse normal method", xlab = "Raw p-values (meta-analysis)")

## Summarize the results. DE are the results from the original analyses of the datasets individually.
# add meta-analysis p-values
H_results_metastats <- data.frame(DE, "DE.fishercomb_adjP"=fishcomb$adjpval, "DE.invnorm_adj_P"=invnormcomb$adjpval)
library(data.table)

setDT(H_results_metastats, keep.rownames = "Gene")  # make the rownames a part of the dataframe with column name "Gene"
head(H_results_metastats)

# Make columns that give 1 for significant or 0 for not significant
H_results <- data.frame(DE, "DE.fishercomb"=ifelse(fishcomb$adjpval<=0.05,1,0), "DE.invnorm"=ifelse(invnormcomb$adjpval<=0.05,1,0))
head(H_results)

## Dealing with conflicts in the sign of differential expression 
signsFC <- mapply(H_FC, FUN=function(x) sign(x))
sumsigns <- apply(signsFC,1,sum)
commonsngFC <- ifelse(abs(sumsigns)==dim(signsFC)[2], sign(sumsigns), 0)
# this produced only 1155 genes
unionDE <- unique(c(fishcomb$DEindices, invnormcomb$DEindicies))
head(unionDE)


#This keeps only the ones that are significant for either fishcomb or invnormcomb and adds the Signs columns
FC.selecDE <- data.frame(H_results[unionDE,], do.call(cbind,H_FC) [unionDE,], signFC=commonsngFC[unionDE], DE=DE[unionDE])
head(FC.selecDE)

# This is the subset of genes  that don't have a conflict (1 or -1)
keepDE <- FC.selecDE[which(abs(FC.selecDE$signFC)==1),]
# these are the genes with conflicting signs between HS2008 and HS2012
conflictDE <- FC.selecDE[which(FC.selecDE$signFC==0),]
dim(FC.selecDE)
dim(keepDE)
dim(conflictDE)
head(keepDE)

table(conflictDE$DE)

head(conflictDE)
dim(conflictDE)

dim(conflictDE$DE)

table(FC.selecDE$DE.invnorm)
table(keepDE$DE.invnorm)

### Merge results and write to  a files
setDT(FC.selecDE, keep.rownames = "Gene")
head(FC.selecDE)

head(H_results_metastats)


# merge top.tables from HS2008 and 2012 - all genes that were not originally filtered otu.
dim(DGESeq_results_standard)
dim(DGESeq_results_stringtie_modified)
dim(DGESeq_results_htseq)
dim(DGESeq_results_htseq_modified)
dim(DGESeq_results_hisat250)
dim(DGESeq_results_star100)
dim(DGESeq_results_star149)

# Merge DGESeq_results_standard and DGESeq_results_stringtie_modified
merged_df <- merge(DGESeq_results_standard, DGESeq_results_stringtie_modified, by = "Gene", all.x = TRUE, suffixes = c(".Std", ".StrMo"))

# Merge additional data frames
additional_dfs <- list(DGESeq_results_htseq, DGESeq_results_htseq_modified)
for (df in additional_dfs) {
  merged_df <- merge(merged_df, df, by = "Gene", all.x = TRUE, suffixes = c(".ht",".htMo",".hi250",".st100",".st149"))
}

additional_dfs<- list(DGESeq_results_hisat250, DGESeq_results_star100)
for (df in additional_dfs) {
  merged_df <- merge(merged_df, df, by = "Gene", all.x = TRUE, suffixes = c(".hi250",".st100"))
}

additional_dfs<- list(DGESeq_results_star149, EDGER0)
for (df in additional_dfs) {
  merged_df <- merge(merged_df, df, by = "Gene", all.x = TRUE, suffixes = c(".st149",".edge0"))
}

dim(merged_df)

#merge original results with results from meta-analysis
Heat_results <- merge(merged_df, FC.selecDE, all.x = TRUE)
Heat_results2 <- merge(Heat_results, H_results_metastats, all.x = TRUE)
head(Heat_results2)

write.table(Heat_results2, file = "Meta_Heat.txt", row.names = F, sep = "\t", quote = F)

#Make summary tables
## Make a table of results
fishcomb_de <- rownames(keepDE)[which(keepDE[,"DE.fishercomb"]==1)]
invnorm_de <- rownames(keepDE)[which(keepDE[,"DE.invnorm"]==1)]
indstudy_de <- list(rownames(keepDE)[which(keepDE[,"DE.Standard"]==1)],rownames(keepDE)[which(keepDE[,"DE.Stringtie_Modified"]==1)],rownames(keepDE)[which(keepDE[,"DE.Htseq"]==1)], rownames(keepDE)[which(keepDE[,"DE.Htseq_modified"]==1)],rownames(keepDE)[which(keepDE[,"DE.Hisat2_50"]==1)],rownames(keepDE)[which(keepDE[,"DE.Star_100"]==1)],rownames(keepDE)[which(keepDE[,"DE.Star_149"]==1)])
IDD.IRR(fishcomb_de,indstudy_de)

IDD.IRR(invnorm_de,indstudy_de)

#Make Upset plot Treatment
## upload packages
library(UpSetR)
library(ggplot2)

#UpSet Plot
#Subset the top.tables to only contain the rows where "adj.P.Val <=0.05
UP_standard_top <- DGESeq_results_standard[which(DGESeq_results_standard$padj<=0.05), ]
dim(UP_standard_top)

UP_str_mod_top <- DGESeq_results_stringtie_modified[which(DGESeq_results_stringtie_modified$padj<=0.05), ]
dim(UP_str_mod_top)

UP_htseq_top <- DGESeq_results_htseq[which(DGESeq_results_htseq$padj<=0.05), ]
dim(UP_htseq_top)

UP_htseq_mod_top <- DGESeq_results_htseq_modified[which(DGESeq_results_htseq_modified$padj<=0.05), ]
dim(UP_htseq_mod_top)

UP_hisat250_top <- DGESeq_results_hisat250[which(DGESeq_results_hisat250$padj<=0.05), ]
dim(UP_hisat250_top)

UP_star100_top <- DGESeq_results_star100[which(DGESeq_results_star100$padj<=0.05), ]
dim(UP_star100_top)

UP_star149_top <- DGESeq_results_star149[which(DGESeq_results_star149$padj<=0.05), ]
dim(UP_star149_top)

UP_edge0 <- EDGER0[which(EDGER0$padj<=0.05), ]
dim(UP_edge0)

UP_edge10 <- EDGER10[which(EDGER10$padj<=0.05), ]
dim(UP_edge10)

UP_edge100 <- EDGER100[which(EDGER100$padj<=0.05), ]
dim(UP_edge100)

UP_edge20 <- EDGER20[which(EDGER20$padj<=0.05), ]
dim(UP_edge20)

UP_edge50 <- EDGER50[which(EDGER50$padj<=0.05), ]
dim(UP_edge50)

UP_dge0 <- DGESeq_results_0[which(DGESeq_results_0$padj<=0.05), ]
dim(UP_dge0)

UP_dge10 <- DGESeq_results_10[which(DGESeq_results_10$padj<=0.05), ]
dim(UP_dge10)

UP_dge100 <- DGESeq_results_100[which(DGESeq_results_100$padj<=0.05), ]
dim(UP_dge100)

UP_dge50 <- DGESeq_results_50[which(DGESeq_results_50$padj<=0.05), ]
dim(UP_dge50)



#converting column of "Gene" in dataframe column into vector by passing as index
Treatment_Standard = UP_standard_top[['Gene']]
Treatment_Stringtie_modified = UP_str_mod_top [['Gene']]
Treatment_Htseq = UP_htseq_top[['Gene']]
Treatment_Htseq_modified = UP_htseq_mod_top [['Gene']]
Treatment_Hisat_modified = UP_hisat250_top [['Gene']]
Treatment_Star_100 = UP_star100_top [['Gene']]
Treatment_Star_149 = UP_star149_top [['Gene']]
Treatment_DGE0 = UP_dge0[['Gene']]
Treatment_DGE10 = UP_dge10[['Gene']]
Treatment_DGE50 = UP_dge50[['Gene']]
Treatment_DGE100 = UP_dge100[['Gene']]
Treatment_EDGE0 = UP_edge0[['Gene']]
Treatment_EDGE10 = UP_edge10[['Gene']]
Treatment_EDGE20 = UP_edge20[['Gene']]
Treatment_EDGE50 = UP_edge50[['Gene']]
Treatment_EDGE100 = UP_edge100[['Gene']]

# Heat MetaAnalysis
## Keep the genes from the meta-analysis that have invernorm adjusted p-value <= 0.05 and the log fold change is going in the same direction in both analyses
UP_MetaA <- Heat_results2[which(Heat_results2$DE.invnorm_adj_P<=0.05), ]
UP_MetaB <- UP_MetaA[which(abs(UP_MetaA$signFC)==1),]
dim(UP_MetaB)



#converting column of "Gene" in dataframe column into vector by passing as index
Treatment_Meta = UP_MetaB[['Gene']]

listHeat <- list(Standard = Treatment_Standard, Stringtie_modified = Treatment_Stringtie_modified,  Htseq_modified =Treatment_Htseq_modified , Hisat2_modified = Treatment_Hisat_modified, Star100 = Treatment_Star_100, Htseq =  Treatment_Htseq, Star_149 = Treatment_Star_149, DGE0=Treatment_DGE0, DGE10=Treatment_DGE10, DGE50=Treatment_DGE50, DGE100=Treatment_DGE100)

listHeat <- list(Standard = Treatment_Standard, Stringtie_modified = Treatment_Stringtie_modified,  Htseq_modified =Treatment_Htseq_modified, Htseq =  Treatment_Htseq) 

# Make upset put and save as a .pdf file
pdf("OutputGraphs/Upset_Treatment.pdf")
#make plot
intersections=list(list(listHeat))
upset(fromList(listHeat), nsets=11,mainbar.y.label = "No. of Differentially Expressed Genes (adjusted p-value=0.05)", sets.x.label = "DEGs",order.by = "freq")
dev.off()

upset(fromList(listHeat), nsets=4,mainbar.y.label = "No. of Differentially Expressed Genes (adjusted p-value=0.05)", sets.x.label = "DEGs",order.by = "freq")
dev.off()

# Get the names of genes for each treatment
gene_names_Standard <- names(Treatment_Standard)
gene_names_Stringtie_modified <- names(Treatment_Stringtie_modified)
gene_names_Htseq_modified <- names(Treatment_Htseq_modified)
gene_names_Hisat2_modified <- names(Treatment_Hisat_modified)
gene_names_Star_100 <- names(Treatment_Star_100)
gene_names_Htseq <- names(Treatment_Htseq)
gene_names_Star_149 <- names(Treatment_Star_149)
gene_names_dge0 <- names(Treatment_DGE0)
gene_names_dge10 <- names(Treatment_DGE10)
gene_names_dge50 <- names(Treatment_DGE50)
gene_names_dge100 <- names(Treatment_DGE100)
gene_names_edge0<- names(Treatment_EDGE0)
gene_names_edge10<- names(Treatment_EDGE10)
gene_names_edge20<- names(Treatment_EDGE20)
gene_names_edge50<- names(Treatment_EDGE50)
gene_names_edge100<- names(Treatment_EDGE100)

# Create a list of gene names
gene_lists <- list(Treatment_Standard, Treatment_Stringtie_modified,  Treatment_Htseq_modified , Treatment_Hisat_modified, Treatment_Star_100, Treatment_Htseq,  Treatment_Star_149, Treatment_DGE0,Treatment_DGE10, Treatment_DGE50, Treatment_DGE100)

gene_lists <- list(Treatment_Htseq_modified, Treatment_Htseq, Treatment_Stringtie_modified ) 
common_genes <- Reduce(intersect, gene_lists)
# Find common gene names
common_genes <- Reduce(intersect, gene_lists)

# Print the common gene names
print(common_genes)





