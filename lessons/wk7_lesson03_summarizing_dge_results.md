---
title: "Summarizing results from the Wald test"
author: "Harvard HPC Staff, adapted by Sally Chang at NICHD"
date: "Last Modified May 2025"
---

Approximate time: 30 minutes

## Learning Objectives

-   Evaluate the number of differentially expressed genes produced for each comparison
-   Construct R objects containing significant genes from each comparison

## **Getting Back up to Speed**

You need to have the `res_tableOE` object in your environment. If you need to be completely caught up, you can copy and paste the following into an R Script and run it. If you don't already have the files in your `/data` directory, please see [Wk 5 Lesson 01](../wk5_lesson01_introR_Rstudio.md) for instructions on where to obtain the input files.

Remember to get your HPC On Demand session going, if applicable, and open your `DEAnalysis` R project!

``` r
# Gene-level differential expression analysis using DESeq2

# Setup
# Bioconductor and CRAN libraries used - already installed on Biowulf
library(tidyverse)
library(RColorBrewer)
library(DESeq2)
library(pheatmap)
library(ggrepel)
library(DEGreport)

# Load in data
data <- read.table("data/mov10_AllSamples_featurecounts.Rmatrix.txt", header=T, row.names=1)

meta <- read.table("data/mov10_AllSamples_metadata.txt", header=T, row.names=1)

# Create DESeq2Dataset object
dds <- DESeqDataSetFromMatrix(countData = data, colData = meta, design = ~ sampletype)

# Run analysis. This does a lot!
dds <- DESeq(dds)

#Set up control vs. OE contrast
contrast_oe <- c("sampletype", "MOV10_overexpression", "control")
res_tableOE <- results(dds, contrast=contrast_oe, alpha = 0.1)

#LFC shrinking
res_tableOE_unshrunken <- res_tableOE
res_tableOE <- lfcShrink(dds, coef="sampletype_MOV10_overexpression_vs_control", type="apeglm")
```

## Summarizing results

To summarize the results table, a handy function in DESeq2 is `summary()`. Confusingly it has the same name as the function used to inspect data frames. This function when called with a DESeq results table as input, will summarize the results using a default threshold of padj \< 0.1, which is comparable to the results you would get from BSPC.

``` r
## Summarize results
summary(res_tableOE, alpha = 0.1)
```

In addition to the number of genes up- and down-regulated at the default threshold, **the function also reports the number of genes that were tested (genes with non-zero total read count), and the number of genes not included in multiple test correction due to a low mean count**.

```         
out of 43674 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)       : 3204, 7.3%
LFC < 0 (down)     : 4284, 9.8%
outliers [1]       : 0, 0%
low counts [2]     : 22701, 52%
(mean count < 5)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results
```

**Based on the contrast we set up - how would you interpret "upregulated"? In which group are those genes upregulated in?**

## Extracting significant DE genes and adding gene symbols

Let's first create variables that contain our threshold criteria. We will only be using the adjusted p-values in our criteria:

``` r
### Set thresholds
padj.cutoff <- 0.1
```

We can easily subset the results table to only include those that are significant using the `filter()` function, and merge a data frame from the

``` r
# Create a dataframe to simplify the res_tableOE object. What do you think is happening in rownames_to_columns()?
res_OE_df <- rownames_to_column(data.frame(res_tableOE), var="gene")
```

*Wouldn't it be nice if we could have gene symbols associated with each of our genes instead of just the Ensembl gene IDs?*

First, let's read in the GTF and extract just the Ensembl Gene IDs and the Gene Symbols so that we can associate those Gene Symbols to your data for better visualizations. **Run readGFF() this time, but in the future and your setup script, just read in the CSV option:**

``` r
## Convert our GTF to a large data frame.
library(rtracklayer)
gtf <- readGFF("/data/Bspc-training/shared/rnaseq_mov10/human_GRCh38/gencode.v47.primary_assembly.annotation.gtf")

# Extract only the columns with ensembl gene names and gene symbols
gtf_names <- gtf %>% dplyr::select(gene_id, gene_name) %>%
  dplyr::distinct() %>% 
  dplyr::rename(ensgene = gene_id, symbol = gene_name)
```

``` r
## Next time, just read in this tab-separated data file
gtf_names <- read.table("/data/Bspc-training/shared/rnaseq_mov10/downstream_data/gtf_names.txt", header=TRUE)
```

> **Discussion**: What are some commands we can use to preview the contents of this table and if it has the \~78k features we expected from our original GTF?

Now we merge these two data frames by their columns that each contain Ensembl Gene IDs:

```{r}
#merge gene symbols
res_OE_df <- merge(gtf_names,res_OE_df, by.x="ensgene", by.y="gene")
```

Now we can subset that table to only keep the significant genes using our pre-defined thresholds:

``` r
# Subset the dataframe to keep only significant genes 
sigOE <- res_OE_df[res_OE_df$padj < padj.cutoff,]
```

``` r
# Take a quick look at this dataframe
head(sigOE)
```

**Discussion**: Do we remember what these columns represent?

-   `baseMean`: mean of normalized counts for all samples
-   `log2FoldChange`: log2 fold change
-   `lfcSE`: standard error
-   `stat`: Wald statistic
-   `pvalue`: Wald test p-value
-   `padj`: BH adjusted p-values

## Writing out our table

Let's say we wanted to save our `res_tableOE_df` object as a table so we could read it in again later. Let's use the following command:

``` r
write.table(res_OE_df, file = "res_OE_df.tsv", sep = "\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
```

**Discussion:** What do the different arguments mean? How can we access help documentation about `write.table()`?

------------------------------------------------------------------------

## **ASSIGNMENT:**

**MOV10 Differential Expression Analysis: Control versus Knockdown**

1.  Using the same p-adjusted threshold as above (`padj.cutoff < 0.1`), subset `res_tableKD` to report the number of genes that are up- and down-regulated in Mov10_knockdown compared to control. Save it as an object called `sigKD`.
2.  How many genes are differentially expressed in the Knockdown compared to Control? How does this compare to the overexpression significant gene list (in terms of numbers)?

------------------------------------------------------------------------

## Wrap-Up

Now that we have extracted the significant results, we are ready for visualization!

In summary, to help prepare for the [next lesson](../lessons/wk7_lesson04_visualizing_results.md), you can add the following to your script along with the code required to create the `sigKD` object:

``` r
### Set thresholds
padj.cutoff <- 0.1

#formatting res_OE_df 
res_OE_df <- rownames_to_column(data.frame(res_tableOE), var="gene")

#import relevant parts of GTF file
gtf_names <- read.table("/data/Bspc-training/shared/rnaseq_jan2025/downstream_data/gtf_names.txt", header=TRUE)

#merge gene symbols
res_OE_df <- merge(gtf_names,res_OE_df, by.x="ensgene", by.y="gene")

# Subset the dataframe to keep only significant genes 
sigOE <- res_OE_df[res_OE_df$padj < padj.cutoff,]
```

------------------------------------------------------------------------

*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*

*Some materials and hands-on activities were adapted from [RNA-seq workflow](http://www.bioconductor.org/help/workflows/rnaseqGene/#de) on the Bioconductor website*

------------------------------------------------------------------------
