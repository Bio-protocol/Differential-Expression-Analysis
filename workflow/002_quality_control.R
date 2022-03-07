# Loading R packages
library(GEOquery)
library(DESeq2)
library(edgeR)
library(limma)
library(ggplot2)
library(pheatmap)
library(Glimma)
library(dplyr)
library(readr)

## load raw_counts data, ignore this step if you has run 001_loading_and_exploring_data file.
setwd("~/path/to/project/DGE/")
load('input/raw_counts.RData')
group <- as.factor(c('mock','mock','mock','hrcc','hrcc','hrcc'))


#------------------------------------------------------------------
##                     Data preprocessing
#--------------------------------------------------------------------

# No preliminary normalization of this data is needed in DEseq2 
colData <- data.frame(row.names=colnames(raw_counts), group= group)

# Create DESeq object
dds <- DESeqDataSetFromMatrix(countData = raw_counts,
                              colData = colData,
                              design = ~ group)

# Retrieve the original count matrix
head(counts(dds))

# Normalization to generate size factors. But this step is also automatically performed by the DESeq() function below.
dds_norm <- estimateSizeFactors(dds)

# Taking a look at the normalization factors of each sample
sizeFactors(dds_norm)

# Retrieving the normalized counts matrix from dds
normalized_counts <- counts(dds_norm, normalized=TRUE)
head(normalized_counts)

### Transform counts for data visualization
rld <- rlog(dds_norm, blind=TRUE) # blind=TRUE argument is to make sure that the rlog() function does not take our sample groups into account


#------------------------------------------------------------------
##                     PCA
#--------------------------------------------------------------------

### Plot PCA 
pdf('graphs/pca.pdf')
plotPCA(rld, intgroup="group")
dev.off()

#------------------------------------------------------------------
##                    multidimensional scaling (MDS) plot
#--------------------------------------------------------------------

# An interactive R widget for generating plots is created and exported as HTML documents.
glimmaMDS(dds)

#------------------------------------------------------------------
##                    Hierarchical Clustering Heatmap
#--------------------------------------------------------------------


### Extract the rlog matrix from the object
rld_mat <- assay(rld)

### Compute pairwise correlation values
rld_cor <- cor(rld_mat)

### Plot heatmap using the correlation matrix and the metadata object
pdf('graphs/heatmap1.pdf')
pheatmap(rld_cor, annotation_col = colData)
dev.off()












