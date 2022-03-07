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

#-------------------------------------
##              edgeR
#----------------------------------------

# Create DEGList object
d.full <- DGEList(counts= raw_counts,group= group)
head(d.full$counts)

# total gene counts per sample
apply(d.full$counts, 2, sum) 

# keeping a gene if it has a cpm of 1 or greater for at least two samples
keep <- rowSums(cpm(d.full)> 10) >= 2
d <- d.full[keep, , keep.lib.sizes=FALSE]

# Reset the library sizes
d$samples$lib.size <- colSums(d$counts)
d$samples

# Normalization: without this, the default value is 1 for all values in d$samples$norm.factors.
d <- calcNormFactors(d)
d$samples$norm.factors


# Specify the model to be fitted
design <- model.matrix(~ 0 + d$samples$group)

# Changing column and row names
rownames(design)<-colnames(d)
colnames(design) <- levels(d$samples$group)

# Fit the common dispersion
d <- estimateGLMCommonDisp(d,design)

# Fit a trended model: the default is "auto" which chooses "bin.spline" when > 200 tags and "power" otherwise.
d <- estimateGLMTrendedDisp(d, design, method="power")

# Fit the tagwise dispersion
d <- estimateGLMTagwiseDisp(d,design)

# tagwise biological coefficient of variation (square root of dispersions) against log2-CPM
pdf('graphs/plotBCV_edger.pdf')
plotBCV(d)
dev.off()

# Fitting the GLM model
fit <- glmFit(d, design)
lrt <- glmLRT(fit,  contrast=c(-1,1))

# Exploring the result for the top 10 genes
topTags(lrt, n=10)

# Extracting all DE genes
edgeR_diff<- as.data.frame(topTags(lrt, n=nrow(d)))

### Set thresholds
fdr.cutoff <- 0.05

# Subset the tibble to keep only significant genes
edgeR_sig <- edgeR_diff[which(edgeR_diff$FDR < fdr.cutoff),]
edgeR_sig<- edgeR_sig[order(edgeR_sig$FDR), ]


# Marking significant genes
de <- decideTestsDGE(lrt, adjust.method="BH", p.value = 0.05)
detags <- rownames(d)[as.logical(de)]

# Plotting the log-fold changes of all the genes
pdf('graphs/maplot_edger.pdf')
plotSmear(lrt, de.tags=detags)
abline(h = c(-2, 2), col = "blue")
dev.off()

write.csv(edgeR_diff, 'output/edgeR_allgenes.csv', row.names = T)
write.csv(edgeR_sig, 'output/edgeR_siggenes.csv', row.names = T)

