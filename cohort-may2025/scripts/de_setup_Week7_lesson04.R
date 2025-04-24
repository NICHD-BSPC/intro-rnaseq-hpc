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

contrast_kd <- c("sampletype", "MOV10_knockdown", "control")
res_tableKD <- results(dds, contrast=contrast_kd, alpha = 0.1)
res_tableKD <- lfcShrink(dds, coef="sampletype_MOV10_knockdown_vs_control", type="apeglm")

res_tableKD_df <- res_tableKD %>%
data.frame() %>%
rownames_to_column(var="gene")

res_tableKD_df <- merge(gtf_names,res_tableKD_df, by.x="ensgene", by.y="gene")

sigKD <- res_tableKD_df %>%
dplyr::filter(padj < padj.cutoff)
############################# Week 7 Lesson 04 ###################################

mov10_meta <- meta %>% 
  rownames_to_column(var="samplename")

# set up normalized counts
normalized_counts <- counts(dds, normalized=T) %>% 
  data.frame() %>%
  rownames_to_column(var="gene") 

normalized_counts <- merge(gtf_names, normalized_counts, by.x="ensgene", by.y="gene")
