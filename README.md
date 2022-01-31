[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](http://www.gnu.org/licenses/gpl-3.0)

# Differential Expression Analysis: Simple pair, Interaction, Time-series


To guide eBook authors having a better sense of the workflow layout, here we briefly introduce the specific purposes of the dir system. 

1. __cache__: Here, it stores R codes for preprocessing Arabidopsis raw time course data.
2. __graphs__: The graphs/figures produced during the analysis.
3. __input__: Here, we store the raw input data, including both for simple pair DGE and time course analysis . 
4. __output__: The final output results of the workflow, including all DE genes and significant DE genes of the three DGE methods.
5. __workflow__: Step by step pipeline for DGE and time course analysis. 


## R packages required

- __Required software and versions__: 
    - [FastQC v0.11.9](http://www.bioinformatics.babraham.ac.uk/projects/download.html#fastqc)
    - [multiqc](https://github.com/ewels/MultiQC)
    - [R 3.6.3](https://cran.r-project.org/) for results ploting
        - [RStudio 1.4](https://rstudio.com/), [ggplot2 3.3.3](https://cran.r-project.org/web/packages/ggplot2/index.html), [tidyr 1.1.2](https://github.com/tidyverse/tidyr)


## Input Data

The example data used here is the paired-end fastq file generated by using Illumina platform.  

- R1 FASTQ file: `input/reads1.fastq`  
- R2 FASTQ file: `input/reads2.fastq`  

Each entry in a FASTQ files consists of 4 lines:  

1. A sequence identifier with information about the sequencing run and the cluster. The exact contents of this line vary by based on the BCL to FASTQ conversion software used.  
2. The sequence (the base calls; A, C, T, G and N).  
3. A separator, which is simply a plus (+) sign.  
4. The base call quality scores. These are Phred +33 encoded, using ASCII characters to represent the numerical quality scores.  

The first entry of the input data:
```
@HWI-ST361_127_1000138:2:1101:1195:2141/1
CGTTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNGGAGGGGTTNNNNNNNNNNNNNNN
+
[[[_BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
```


## Major steps

#### Step 1: running the FastQC to conduct quality checking
- Note that you have to normalize the path in the shell script.

```
sh workflow/1_run_fastqc.sh
```

#### Step 2: aggregate results from FastQC

```
sh workflow/2_aggregate_results.sh
```

#### Step 3: view the results

- Results can be visualized by clicking `output/multiqc_report.html`.
- Alternatively, you can plot the results yourself using the below R code.

```
3_visualize_results.Rmd
```

## Expected results

![](graphs/figure1.png)

## License
It is a free and open source software, licensed under []() (choose a license from the suggested list:  [GPLv3](https://github.com/github/choosealicense.com/blob/gh-pages/_licenses/gpl-3.0.txt), [MIT](https://github.com/github/choosealicense.com/blob/gh-pages/LICENSE.md), or [CC BY 4.0](https://github.com/github/choosealicense.com/blob/gh-pages/_licenses/cc-by-4.0.txt)).
