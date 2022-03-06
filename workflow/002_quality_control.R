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

#-------------------------------------------------
## Downloading and importing of RNA-seq count data
#-------------------------------------------------

# Specify URL where file is stored
url <- "http://bioinf.wehi.edu.au/edgeR/UserGuideData/arab.rds"

# Specify destination where file should be saved
setwd("~/path/to/project/DEG/")
destfile <- "~/path/to/arab.rds"

# Apply download.file function in R
download.file(url, destfile)

# Import the input r dataset
raw_counts <- readRDS(destfile)

# Check out the import raw counts matrix
head(raw_counts)

#Create group vector
group <- c('mock','mock','mock','hrcc','hrcc','hrcc')


#------------------------------------------------------------------
##                     Data preprocessing
#--------------------------------------------------------------------

# No preliminary normalization of this data is needed in DEseq2 
colData <- data.frame(row.names=colnames(raw_counts), group=group)

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
plotPCA(rld, intgroup="group")


#------------------------------------------------------------------
##                    multidimensional scaling (MDS) plot
#--------------------------------------------------------------------

glimmaMDS(dds)


#------------------------------------------------------------------
##                    Hierarchical Clustering Heatmap
#--------------------------------------------------------------------


### Extract the rlog matrix from the object
rld_mat <- assay(rld)

### Compute pairwise correlation values
rld_cor <- cor(rld_mat)

### Plot heatmap using the correlation matrix and the metadata object
pheatmap(rld_cor, annotation_col = colData)













