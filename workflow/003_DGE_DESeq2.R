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

#-------------------------------------
##              DESeq2 
#----------------------------------------

## Run pipeline for differential expression steps
dds <- DESeq(dds)

## Plot dispersion estimates
pdf('graphs/dispersion_deseq2.pdf')
plotDispEsts(dds)
dev.off()

## Extracting results
res <- results(dds, contrast=c("group","mock","hrcc"))

# View summary of results
summary(res)
summary(res, alpha = 0.05)

# Get information on each column in results
mcols(res, use.names=T)


## Extracting results
res <- results(dds, contrast=c("group","mock","hrcc"))

# View summary of results
summary(res)
summary(res, alpha = 0.05)

# Get information on each column in results
mcols(res, use.names=T)


### Filtering the data

diff<- data.frame(res)
dim(diff)

# 1. Genes with zero counts in all samples
# Filter genes by zero expression
if(sum(which(diff$baseMean == 0))== 0){
  diff<- diff
  print('no filtering for zero expression')
}else{
  diff<- diff[-which(diff$baseMean == 0),]
}

# 2. Filter genes that have an extreme outlier

if(sum(which(is.na(diff$pvalue) & is.na(diff$padj) & diff$baseMean > 0))== 0){
  diff<- diff
  print('no filtering for extreme outliers')
}else{
  diff<- diff[-which(is.na(diff$pvalue) & is.na(diff$padj) & diff$baseMean > 0),]
}

# 3. Filter genes below the low mean threshold
if(sum(which(!is.na(diff$pvalue) & is.na(diff$padj) & diff$baseMean > 0))== 0){
  diff<- diff
  print('no filtering for low mean threshold')
}else{
  diff<- diff[-which(!is.na(diff$pvalue) & is.na(diff$padj) & diff$baseMean > 0),]
}

dim(diff)


### Set thresholds
padj.cutoff <- 0.05

# Subset the tibble to keep only significant genes
sig <- diff[which(diff$padj < padj.cutoff),]
sig<- sig[order(sig$padj), ]

# Take a quick look at this
dim(diff)
dim(sig)
head(sig)

write.csv(diff, 'output/deseq2_allgenes.csv', row.names = T)
write.csv(sig, 'output/deseq2_siggenes.csv', row.names = T)

### Visualizing the results

# MA plot

pdf('graphs/maplot.pdf')
plotMA(res)
dev.off()

# Extract normalization values
normalized_counts <- data.frame(counts(dds, normalized=T))

# Significantly normalized genes
norm_sig<- normalized_counts[match(rownames(sig), rownames(normalized_counts)),]

### Run pheatmap 
pdf('graphs/heatmap2.pdf')
pheatmap(norm_sig, 
         cluster_rows = T, 
         show_rownames = F,
         annotation = colData, 
         border_color = NA, 
         fontsize = 10, 
         scale = "row", 
         fontsize_row = 10, 
         height = 20)
dev.off()

## Volcano plot

## Obtain logical vector where TRUE values denote padj values < 0.05 and fold change > 1.5 in either direction
vol<- diff %>% mutate(threshold = padj < 0.05 & abs(log2FoldChange) >= 0.5)

pdf('graphs/Volcanoplot.pdf')
ggplot(vol) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold)) +
  ggtitle("overexpression") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  #scale_y_continuous(limits = c(0,50)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))  
dev.off()



