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

(3) Set up your own local DEAnalysis RProject, probably in its own new directory
