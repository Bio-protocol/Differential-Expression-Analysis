[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](http://www.gnu.org/licenses/gpl-3.0)

# Differential Expression Analysis: Simple pair, Interaction, Time-series

Identifying differentially expressed (DE) genes across specific conditions is vital in understanding phenotypic variation. The fast-growing RNA sequencing provides much information that efficiently quantifies gene expressions. Methods and tools dedicated to differential gene expression analysis from RNA-seq data also increased rapidly. More than 30 DE methods have been published; however, many comparison studies spotlight that no single method outperforms others in all circumstances. In this study, we test and compare the performances of three widely used R packages: edgeR, DESeq2, and limma voom with Arabidopsis thaliana data. Even though the standard DE analysis has been extensively used and improved over the past years, time course RNA-seq can also provide an advanced understanding of gene regulation, biological development, and identifying biologically DE genes. Therefore, we also conducted a time course analysis using another Arabidopsis time course dataset. These methods are initiated in separate R packages, then detailed R codes and explanations are constructed to help build a more convenient user experience.

To guide eBook authors having a better sense of the workflow layout, here we briefly introduce the specific purposes of the dir system. 

1. __cache__: Here, it stores R codes for preprocessing Arabidopsis raw time course data.
2. __graphs__: The graphs/figures produced during the analysis.
3. __input__: Here, we store the raw input data, including both for simple pair DGE and time course analysis . 
4. __output__: The final output results of the workflow, including all DE genes and significant DE genes of the three DGE methods.
5. __workflow__: Step by step pipeline for DGE and time course analysis. 


## R packages required
R version 4.1.1 (2021-08-10)

- __Required R packages and versions__: 
    - [ggplot2 3.3.3](https://cran.r-project.org/web/packages/ggplot2/index.html), [dplyr 1.0.7](https://dplyr.tidyverse.org/), [GEOquery 2.60.0](https://www.bioconductor.org/packages/release/bioc/html/GEOquery.html), [DESeq2 1.32.0](https://bioconductor.org/packages/release/bioc/html/DESeq2.html), [edgeR 3.34.0](https://bioconductor.org/packages/release/bioc/html/edgeR.html), [limma 3.48.0](https://bioconductor.org/packages/release/bioc/html/limma.html), [pheatmap 1.0.12](https://cran.r-project.org/web/packages/pheatmap/index.html), [Glimma 2.2.0](https://bioconductor.org/packages/release/bioc/html/Glimma.html), [readr 1.4.0](https://readr.tidyverse.org/)


## Input Data
- __Case study 1: A comparison of three methods for DGE analysis__

To demonstrate, here we use the Arabidopsis thaliana RNA-Seq data published by Cumbie et al., (Cumbie et al., 2011). Summarized count data is available as an R dataset, and readers can download the data from the input folder (arab.rds). In Cumbie’s experiment, they inoculate six-week-old Arabidopsis plants with the mutant of P.syringae. Control plants were inoculated with a mock pathogen. Each treatment was done as biological triplicates, with each pair of replicates done at separate times and derived from independently grown plants and bacteria.

- __Case study 2: Time course analysis__

Here we demonstrate a fundamental time course analysis with an Arabidopsis dataset containing gene counts for an RNA-seq time course. This experiment aims to see if the differentiated endodermal cells have a distinct transcriptional response to auxin treatment. They performed a time series of 10µM NAA treatment and sample at t= 0, 2, 4, 8, 16, and 24hrs after NAA treatment (Ursache et al., 2021). For the time series, they compared roots of the solitary root 1 (slr-1) mutant to the CASP1::shy2-2/slr-1 double mutant. The raw data from the NCBI database (Ursache et al., 2021) was processed and saved as a RangedSummarizeExperiment RData file. The processed data can be downloaded from the input folder (arab_time.Rdata).

## Major steps

#### Step 1: Download and explore RNA-seq data

```
Rscript workflow/001_loading_and_exploring_data.R
```

#### Step 2: Quality control

```
Rscript workflow/002_quality_control.R
```

#### Step 3: DGE analysis--DESeq2

```
Rscript workflow/003_DGE_DESeq2.R
```

#### Step 4: DGE analysis--edgeR

```
Rscript workflow/004_DGE_edgeR.R
```

#### Step 5: DGE analysis--limma voom

```
Rscript workflow/005_DGE_limma_voom.R
```

#### Step 6: Time course analysis

```
Rscript workflow/006_time_course.R
```
## Expected results

Expected results are stored in output and graphs folders. A several examples of outputs:

![PCA plot](graphs/pca.pdf)
![MA plot](graphs/maplot.pdf)
![Volvano plot](graphs/Volcanoplot.pdf)
![Line grapg of time course](graphs/time_course.pdf)
![Heatmap of time course](graphs/heatmap_time.pdf)


## License
It is a free and open source software, licensed under []() [GPLv3](https://github.com/github/choosealicense.com/blob/gh-pages/_licenses/gpl-3.0.txt).
