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

#-------------------------------------
##              Time course
#----------------------------------------

# Loading the data
load("input/arab_time.RData")

# Checking the data
head(arab_time)
dim(assays(arab_time)$counts)
assays(arab_time)$counts[1:5,1:5]

# Designing formula that models the strain difference
ddsTC <- DESeqDataSet(arab_time, ~ strain + hour + strain:hour)


# Performing a likelihood ratio test, where we remove the strain-specific differences over time
ddsTC <- DESeq(ddsTC, test="LRT", reduced = ~ strain + hour)

# Extracting results and printing out the top 4 significant DE genes
resTC <- results(ddsTC)
resTC$symbol <- mcols(ddsTC)$symbol
head(resTC[order(resTC$padj),], 4)

# Getting input data for line figure
arab <- plotCounts(ddsTC, which.min(resTC$padj), 
                   intgroup = c("hour","strain"), returnData = TRUE)
arab$hour <- as.numeric(as.character(arab$hour))

# Plotting the counts for the groups over time
pdf('graphs/time_course.pdf')
ggplot(arab, aes(x = hour, y = count, color = strain, group = strain)) + 
  geom_point() + stat_summary(fun.y=mean, geom="line") +
  scale_y_log10()
dev.off()

# show some results
resultsNames(ddsTC)
res <- results(ddsTC, name="hour_4_vs_0", test="Wald")
res[which.min(resTC$padj),]

# Extracting a matrix of the shrunken log2 fold changes
betas <- coef(ddsTC)
colnames(betas)

# Top 20 significant DE genes
topGenes <- head(order(resTC$padj),20)

# log2 fold change between (-3,3)
mat <- betas[topGenes, -c(1,2)]
thr <- 3 
mat[mat < -thr] <- -thr
mat[mat > thr] <- thr

# Plotting
pdf('graphs/heatmap_time.pdf')
pheatmap(mat, breaks=seq(from=-thr, to=thr, length=101),
         cluster_col=FALSE) 
dev.off()
