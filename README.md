[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](http://www.gnu.org/licenses/gpl-3.0)

# Differential Expression Analysis: Simple pair, Interaction, Time-series

Identifying differentially expressed (DE) genes across specific conditions is vital in understanding phenotypic variation. The fast-growing RNA sequencing provides much information that efficiently quantifies gene expressions. Methods and tools dedicated to differential gene expression analysis from RNA-seq data also increased rapidly. More than 30 DE methods have been published; however, many comparison studies spotlight that no single method outperforms others in all circumstances. In this study, we test and compare the performances of three widely used R packages: edgeR, DESeq2, and limma voom with Arabidopsis thaliana data. Even though the standard DE analysis has been extensively used and improved over the past years, time course RNA-seq can also provide an advanced understanding of gene regulation, biological development, and identifying biologically DE genes. Therefore, we also conducted a time course analysis using another Arabidopsis time course dataset. These methods are initiated in separate R packages, then detailed R codes and explanations are constructed to help build a more convenient user experience.

To guide eBook authors having a better sense of the workflow layout, here we briefly introduce the specific purposes of the dir system. 

1. __cache__: Here, it stores R codes for preprocessing Arabidopsis raw time course data.
2. __graphs__: The graphs/figures produced during the analysis.
3. __input__: Here, we store the raw input data, including both for simple pair DGE and time course analysis . 
4. __output__: The final output results of the workflow, including all DE genes and significant DE genes of the three DGE methods.
5. __workflow__: Step by step pipeline for DGE and time course analysis. 

## Workflow
![workflow](workflow/workflow.png)

## R packages required
R version 4.1.1 (2021-08-10)

- __Required R packages and versions__: 
    - [ggplot2 3.3.3](https://cran.r-project.org/web/packages/ggplot2/index.html), [dplyr 1.0.7](https://dplyr.tidyverse.org/), [GEOquery 2.60.0](https://www.bioconductor.org/packages/release/bioc/html/GEOquery.html), [DESeq2 1.32.0](https://bioconductor.org/packages/release/bioc/html/DESeq2.html), [edgeR 3.34.0](https://bioconductor.org/packages/release/bioc/html/edgeR.html), [limma 3.48.0](https://bioconductor.org/packages/release/bioc/html/limma.html), [pheatmap 1.0.12](https://cran.r-project.org/web/packages/pheatmap/index.html), [Glimma 2.2.0](https://bioconductor.org/packages/release/bioc/html/Glimma.html), [readr 1.4.0](https://readr.tidyverse.org/)


## Input Data
- __Case study 1: A comparison of three methods for DGE analysis__

To demonstrate, here we use the Arabidopsis thaliana RNA-Seq data published by Cumbie et al., (Cumbie et al., 2011). Summarized count data is available as an R dataset, and readers can download the data from the input folder (arab.rds). In Cumbie???s experiment, they inoculate six-week-old Arabidopsis plants with the mutant of P.syringae. Control plants were inoculated with a mock pathogen. Each treatment was done as biological triplicates, with each pair of replicates done at separate times and derived from independently grown plants and bacteria.

- __Case study 2: Time course analysis__

Here we demonstrate a fundamental time course analysis with an Arabidopsis dataset containing gene counts for an RNA-seq time course. This experiment aims to see if the differentiated endodermal cells have a distinct transcriptional response to auxin treatment. They performed a time series of 10??M NAA treatment and sample at t= 0, 2, 4, 8, 16, and 24hrs after NAA treatment (Ursache et al., 2021). For the time series, they compared roots of the solitary root 1 (slr-1) mutant to the CASP1::shy2-2/slr-1 double mutant. The raw data from the NCBI database (Ursache et al., 2021) was processed and saved as a RangedSummarizeExperiment RData file. The processed data can be downloaded from the input folder (arab_time.Rdata).

## Major steps

#### Step 1: Download and explore RNA-seq data

Note: Users need to modify a few path variables before running the code.
```
#Line 20
setwd("~/path/to/project/DGE")
```
Run 001_loading_and_exploring_data code:
```
source('workflow/001_loading_and_exploring_data.R')
```

#### Step 2: Quality control
Note: Users need to modify a few path variables before running the code.
```
#Line 13
setwd("~/path/to/project/DGE")
```
Run 002_quality_control code:
```
source('workflow/002_quality_control.R')
```

#### Step 3: DGE analysis--DESeq2
Note: Users need to modify a few path variables before running the code.
```
#Line 13
setwd("~/path/to/project/DGE/")
```
Run 003_DGE_DESeq2 code:
```
source('workflow/003_DGE_DESeq2.R')
```

#### Step 4: DGE analysis--edgeR
Note: Users need to modify a few path variables before running the code.
```
#Line 13
setwd("~/path/to/project/DGE/")
```
Run 004_DGE_edgeR code:
```
source('workflow/004_DGE_edgeR.R')
```

#### Step 5: DGE analysis--limma voom
Note: Users need to modify a few path variables before running the code.
```
#Line 13
setwd("~/path/to/project/DGE/")
```
Run 005_DGE_limma_voom code:
```
source('workflow/005_DGE_limma_voom.R')
```

#### Step 6: Time course analysis
Note: Users need to modify a few path variables before running the code.
```
#Line 15
setwd("~/path/to/project/DGE/")
```
Run 006_time_course code:
```
source('workflow/006_time_course.R')
```
## Expected results

Note: Figures generated during DGE analysis (like PCA, heatmap, MA plot, e.g.) and significant gene lists are stored in the graphs and output folder. Since the graphs are saved as 'pdf' files, they are not supported to display in README. Please go through the 'graphs' folder and the protocol for detailed checking.

Visualization of DGE results using the three selected DGE tools provides valuable insights into their generated results. As seen in Figure below (A, B), the three methods detect similar numbers of DE genes, and most of them are identical. When we set the FDR cutoff to 0.05, they would behave differently. DESeq2 detects the most DE genes, while the limma voom detects the least. However, the significant genes detected by limma voom could almost be found by edgeR and DESeq2 (Figure C). Then we could infer that limma voom is more critical than the other two when calling significant DE genes.

![comparison of genes](graphs/15abc.png)

Let us look at the detected fold changes from all three methods. Here, the genes are colored differently to label which methods find them significant. Suppose that a gene is colored yellow-wish green means both methods find it. Matching the results in Figure above, the outputs of DESeq2 and edgeR are very close indeed since their fold changes correlate much better than the other two methods (Figure below). All three methods behaved well in finding the significant DE genes, while DESeq2 caught more candidates than the other two. It turns out that these genes only detected by DESeq2 have pretty low counts. The DESeq2 goes through the logic of independent filtering within results() to save from multiple test corrections on genes with no power, showing that the likelihood of a gene being significantly differentially expressed is related to how strongly it is expressed. It advocates discarding extremely lowly expressed genes because the differential expression is likely not statistically detectable, so it is unnecessary to pre-filter low count genes as recommended in the vignette. However, in edgeR and limma voom, we performed a minimal count-based pre-filtering to keep only rows with at least a CPM of ten for at least two samples total. After filtering, the candidate genes with apparently moderate fold changes detected by DESeq2 are removed in edgeR and limma voom. We suggest users try different methods with variable parameters and then choose the most suitable one.

![more comparison](graphs/abc.png)

## License
It is a free and open source software, licensed under []() [GPLv3](https://github.com/github/choosealicense.com/blob/gh-pages/_licenses/gpl-3.0.txt).
