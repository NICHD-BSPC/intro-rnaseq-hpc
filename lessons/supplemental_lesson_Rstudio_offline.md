# So, you want to work on RStudio on your local computer

### Installation Requirements

(1) Download the most recent versions of R and RStudio for your laptop:

-   [R](http://lib.stat.cmu.edu/R/CRAN/) (version 4.0.0 or above)
-   [RStudio](https://www.rstudio.com/products/rstudio/download/#download)

> **Note 1:**  When installing the following packages, if you are asked to select (a/s/n) or (y/n), please select “a” or "y" as applicable.

> **Note 2**: If you have a Mac with an M1 chip, download and install this tool before intalling your packages: <https://mac.r-project.org/tools/gfortran-12.2-universal.pkg>

(2) Install the below packages on your laptop from CRAN. You DO NOT have to go to the CRAN webpage; you can use the `install.packages` function to install them one by one.

``` r
install.packages("insert_package_name_in_quotations")
install.packages("insert_package_name_in_quotations")
& so on ...
```

``` r
# Install these packages using install.packages()
BiocManager 
tidyverse 
RColorBrewer 
pheatmap 
ggrepel 
cowplot
```

(2) Install the below packages from Bioconductor. Load BiocManager, then run BiocManager's `install()` function for the only package not already installed by default in our interactive RStudio session:

``` r
library(BiocManager)
# One example if installing packages from Bioconductor now that BiocManager is loaded
install("DEGreport")
```

Note that these package names are case sensitive!

``` r
DESeq2 
clusterProfiler 
DOSE 
org.Hs.eg.db 
pathview 
tximport 
AnnotationHub 
ensembldb
apeglm
DEGreport
```

(3) Set up your own local DEAnalysis RProject, probably in its own new directory. The instructions from [Week 5 Lesson 01](https://nichd-bspc.github.io/intro-rnaseq-hpc/lessons/wk5_lesson01_introR_Rstudio.html) are still mostly relevant, but you will need to make sure the paths and and directories are adjusted for your local operating system.

4.  Get the data minimally necessary for re-creating your downstream differential expression analyses:

    -   [Metadata](../data/mov10_AllSamples_metadata.txt)
    -   [Read Counts](../data/mov10_AllSamples_featurecounts.Rmatrix.txt)
    -   [Gene Symbols](../data/gtf_names.txt) data frame, so you don't need to mess with a whole GTF

5.  Optionally, use a script to get yourself back up to date, such as this one from Week 7:

``` r
# Gene-level differential expression analysis using DESeq2

# Setup
# Bioconductor and CRAN libraries used - already installed on Biowulf
library(tidyverse)
library(RColorBrewer)
library(DESeq2)
library(pheatmap)
library(ggrepel)


# installing DEGreport
# library(BiocManager)
# install("DEGreport")

# load DEGreport into environment
library(DEGreport)

# Load in data
data <- read.table("data/mov10_AllSamples_featurecounts.Rmatrix.txt", header=T, row.names=1)

meta <- read.table("data/mov10_AllSamples_metadata.txt", header=T, row.names=1)

# Create DESeq2Dataset object
dds <- DESeqDataSetFromMatrix(countData = data, colData = meta, design = ~ sampletype)

## This would be a good step to add to both of our scripts
# Run analysis. This does a lot!
dds <- DESeq(dds)

########## Week 7: Lesson 02 #########################################

contrast_oe <- c("sampletype", "MOV10_overexpression", "control")
res_tableOE <- results(dds, contrast=contrast_oe, alpha = 0.1)
res_tableOE_unshrunken <- res_tableOE
res_tableOE <- lfcShrink(dds, coef="sampletype_MOV10_overexpression_vs_control", type="apeglm")

############ Week 7 Lesson 03#############################
padj.cutoff <- 0.1

res_tableOE_df <- res_tableOE %>%
data.frame() %>%
rownames_to_column(var="gene")

## Next time, just read in this TSV
gtf_names <- read.table("/data/Bspc-training/shared/rnaseq_jan2025/downstream_data/gtf_names.txt", header=TRUE)

# merge gene names
res_tableOE_df <- merge(gtf_names,res_tableOE_df, by.x="ensgene", by.y="gene")

sigOE <- res_tableOE_df %>%
  dplyr::filter(padj < 0.1)
```

6.  As a bonus, here is a [R-focused workshop from Harvard HPC](https://hbctraining.github.io/Intro-to-R-flipped/) that assumes you are working locally. In case you want to pick up some more important R Skills!
