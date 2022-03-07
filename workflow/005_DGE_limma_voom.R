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
##              limma voom
#----------------------------------------

# Same as edgeR
dge.full <- DGEList(counts=raw_counts)
keep <- rowSums(cpm(dge.full)>10) >= 2
dge <- dge.full[keep, , keep.lib.sizes=FALSE]
dge$samples$lib.size <- colSums(dge$counts)
dge <- calcNormFactors(dge, method = "TMM")

design <- model.matrix(~0+factor(group))
colnames(design)=levels(factor(group))
rownames(design)=colnames(raw_counts)

### Voom transformation
pdf('graphs/voom.pdf')
v <- voom(dge, design, plot=TRUE, normalize="quantile")
dev.off()

# Fitting a linear model using weighted least squares for each gene
fit <- lmFit(v, design)

# Comparison between two groups
cont.matrix<- makeContrasts(contrasts=c('hrcc-mock'),levels = design)

# Estimating contrast for each gene
vfit<- contrasts.fit(fit,cont.matrix)

# Empirical Bayes smoothing of standard errors
efit<- eBayes(vfit)

# Extracting all DE genes: the number of top genes displayed can be specified, where n=Inf includes all genes.
dt<- decideTests(efit)
tempOutput<-  topTable(efit, coef='hrcc-mock', n=Inf)
limma_diff<- na.omit(tempOutput)
head(limma_diff)

### Set thresholds
adj.p.cutoff <- 0.05

# Subset the tibble to keep only significant genes
limma_sig <- limma_diff[which(limma_diff$adj.P.Val < adj.p.cutoff),]
limma_sig<- limma_sig[order(limma_sig$adj.P.Val), ]

write.csv(limma_diff, 'output/limma_allgenes.csv', row.names = T)
write.csv(limma_sig, 'output/limma_siggenes.csv', row.names = T)

pdf('graphs/plotMD_limma.pdf')
plotMD(efit, column=1, status=dt[,1], main=colnames(efit)[1], xlim=c(-8,13))
dev.off()
