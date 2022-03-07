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
save(raw_counts, file='~/path/to/input/folder/raw_counts.RData')

# Check out the import raw counts matrix
head(raw_counts)

# Create group vector
group <- c('mock','mock','mock','hrcc','hrcc','hrcc')

#------------------------------------------
## Characteristics of RNA-seq count data
#-----------------------------------------

# Histogram for a single sample ('mock1')
pdf('graphs/histogram.pdf')
ggplot(data.frame(raw_counts)) +
  geom_histogram(aes(x = mock1), stat = "bin", bins = 200) +
  xlab("Raw counts") +
  ylab("Number of genes")
dev.off()

# Mean versus variance
# compute a vector of mean values
mean_counts <- apply(raw_counts[,1:3], 1, mean)

# Compute a vector of variance values
variance_counts <- apply(raw_counts[,1:3], 1, var)

# Plot the relationship between mean and variance
df <- data.frame(mean_counts, variance_counts)

pdf('graphs/mean_variance.pdf')
ggplot(df) +
  geom_point(aes(x=(mean_counts+1), y=(variance_counts+1))) + 
  scale_y_log10(limits = c(1,1e9)) +
  scale_x_log10(limits = c(1,1e9)) +
  geom_abline(intercept = 0, slope = 1, color="red")
dev.off()


