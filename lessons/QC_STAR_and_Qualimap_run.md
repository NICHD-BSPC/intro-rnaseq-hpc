---
title: "QC with STAR and Qualimap"
author: "Harvard HPC Staff, Modified by Sally Chang at NICHD"
date: Last Modifed February 2025
---

Approximate time: 50 minutes

## Learning Objectives:

-   Brief explanation about SAM file
-   Running an alignment tool to generate BAM files
-   Running Qualimap to compute metrics on alignment files

## Aligning Reads to a Genome

**In order for us to test for differential expression we need to gain information we will gain genomic coordinate information for where each read maps** by aligning it to the reference genome to produce a BAM file. We will also **assess the quality of the mapping** by using the BAM file as input to a tool called [Qualimap](http://qualimap.bioinfo.cipf.es/doc_html/intro.html) which computes various quality metrics such as DNA or rRNA contamination, 5'-3' biases, and coverage biases.

To determine where on the human genome our reads originated from, we will align our reads to the reference genome using [STAR](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3530905/) (Spliced Transcripts Alignment to a Reference). STAR is an aligner designed to specifically address many of the challenges of RNA-seq data mapping using a strategy to account for spliced alignments.

## Alignment file format: SAM/BAM

The **Sequence Alignment Map format** (SAM) file is **a tab-delimited text file that contains all information from the FASTQ file, with additional fields containing alignment information for each read**. Specifically, we can obtain the genomic coordinates of where each read maps to in the genome and the quality of that mapping. A **BAM file is the binary, compressed version of the SAM file**. It is significantly smaller in size and is usually the file format requested for by downstream tools that require alignment data as input. The paper by [Heng Li et al](http://bioinformatics.oxfordjournals.org/content/25/16/2078.full) provides a lot more detail on the specification, and we will go into detail in the next lesson about how to

![SAM1](../img/sam_bam.png)

## In Brief: STAR Alignment Strategy

STAR is shown to have high accuracy and outperforms other aligners by more than a factor of 50 in mapping speed, but it is memory intensive. The algorithm achieves this highly efficient mapping by performing a two-step process:

1.  Seed searching
2.  Clustering, stitching, and scoring

#### Seed searching

For every read that STAR aligns, STAR will search for the longest sequence that exactly matches one or more locations on the reference genome. These longest matching sequences are called the Maximal Mappable Prefixes (MMPs):

![](../img/alignment_STAR_step1.png){alt="STAR_step1"}

The different parts of the read that are mapped separately are called 'seeds'. So the first MMP that is mapped to the genome is called *seed1*.

STAR will then search again for only the unmapped portion of the read to find the next longest sequence that exactly matches the reference genome, or the next MMP, which will be *seed2*.

![](../img/alignment_STAR_step2.png){alt="STAR_step2"}

This sequential searching of only the unmapped portions of reads underlies the efficiency of the STAR algorithm. STAR uses an uncompressed suffix array (SA) to efficiently search for the MMPs, this allows for quick searching against even the largest reference genomes. Other slower aligners use algorithms that often search for the entire read sequence before splitting reads and performing iterative rounds of mapping.

**If STAR does not find an exact matching sequence** for each part of the read due to mismatches or indels, the previous MMPs will be extended.

![STAR_step3](../img/alignment_STAR_step3.png){alt="STAR_step3"}

**If extension does not give a good alignment**, then the poor quality or adapter sequence (or other contaminating sequence) will be soft clipped.

![STAR_step4](../img/alignment_STAR_step4.png){alt="STAR_step4"}

#### Clustering, stitching, and scoring

The separate seeds are stitched together to create a complete read by first clustering the seeds together based on proximity to a set of 'anchor' seeds, or seeds that are not multi-mapping.

Then the seeds are stitched together based on the best alignment for the read (scoring based on mismatches, indels, gaps, etc.).

![STAR_step5](../img/alignment_STAR_step5.png){alt="STAR_step5"}

## Sinteractive with scratch space

#### Aligning reads

Since we already have the reference index ready, we can move on to aligning reads to the genome.

Create an output directory for our alignment files:

``` bash
$ cd ~/rnaseq/raw_data

$ mkdir ../results/STAR
```

For now, we're going to work on just one sample to set up our workflow. To start we will use the first replicate in the Mov10 over-expression group, `Mov10_oe_1.subset.fq`. Details on STAR and its functionality can be found in the [user manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf); we encourage you to peruse through to get familiar with all available options.

The basic options for aligning reads to the genome using STAR are:

-   `--runThreadN`: number of threads / cores
-   `--readFilesIn`: /path/to/FASTQ_file
-   `--genomeDir`: /path/to/genome_indices_directory
-   `--outFileNamePrefix`: prefix for all output files

Listed below are additional parameters that we will use in our command:

-   `--outSAMtype`: output filetype (SAM default), but also produce a BAM that has been sorted by genomic coordinated
-   `--outSAMunmapped`: what to do with unmapped reads

> **NOTE:** Note that "**STAR’s default parameters are optimized for mammalian genomes.** Other species may require significant modifications of some alignment parameters; in particular, the maximum and minimum intron sizes have to be reduced for organisms with smaller introns" [[1](http://bioinformatics.oxfordjournals.org/content/early/2012/10/25/bioinformatics.bts635.full.pdf+html)].

The full command is provided below for you to copy paste into your terminal. If you want to manually enter the command, it is advisable to first type out the full command in a text editor (i.e. [Sublime Text](http://www.sublimetext.com/) or [Notepad++](https://notepad-plus-plus.org/)) on your local machine and then copy paste into the terminal. This will make it easier to catch typos and make appropriate changes.

``` bash
# assumes you are in raw_data directory
$ STAR --genomeDir /data/NICHD-core0/references/human/gencode-v28/genome/star/human_gencode-v28 \
--runThreadN 6 \
--readFilesIn Mov10_oe_1.subset.fq \
--outFileNamePrefix ../results/STAR/Mov10_oe_1_ \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard 
```

Since this will take longer with a full-size human reference genome, let's submit an SBATCH script:

``` bash
#!/bin/bash

#SBATCH --partition=quick
#SBATCH --time=03:00:00       # time limit
#SBATCH --cpus-per-task=6        # number of cores
#SBATCH --mem=32g   # requested memory
#SBATCH --job-name salmon_in_serial      # Job name
#SBATCH -o %j.out           # File to which standard output will be written
#SBATCH -e %j.err       # File to which standard error will be written
#SBATCH --mail-type=BEGIN,END

# Load Salmon module
module load star/2.7.6a

# Change directory to where the data is
cd /data/NICHD-core0/test/changes/rc_training/rnaseq/raw_data

STAR --genomeDir /data/NICHD-core0/references/human/gencode-v28/genome/star/human_gencode-v28 --runThreadN 6 --readFilesIn Mov10_oe_1.subset.fq --outFileNamePrefix ../results/STAR/Mov10_oe_1_ --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard 
```

#### STAR output

After running our single FASTQ file through the STAR aligner, you should have a number of output files in the `~/rnaseq/results/STAR` directory. Let's take a quick look at some of the files that were generated and explore their content.

``` bash
$ cd ~/rnaseq/results/STAR
    
$ ls -lh
```

You should have **5 output files** plus a single `.tmp` directory for the Mov10_oe_1 sample. The contents of the 5 files are described below:

-   `Log.final.out` - a summary of mapping statistics for the sample
-   `Aligned.sortedByCoord.out.bam` - the aligned reads, sorted by coordinate, in BAM format
-   `Log.out` - a running log from STAR, with information about the run
-   `Log.progress.out` - job progress with the number of processed reads, % of mapped reads etc., updated every \~1 minute
-   `SJ.out.tab` - high confidence collapsed splice junctions in tab-delimited format. Only junctions supported by uniquely mapping reads are reported

**We are most interested in the BAM file, which will allow us to proceed with Qualimap.**

### Running Qualimap

Alignment data files frequently contain biases that are introduced by sequencing technologies and/or during sample preparation. Therefore, one of the fundamental requirement during analysis of these data is to perform quality control. In this way, we get an idea of **how well our reads align to the reference** and **how well data fit with the expected outcome**.

We will take the BAM file we generated in the previous step and use it as input to Qualimap which computes various quality metrics such as DNA or rRNA contamination, 5'-3' biases, and coverage biases.

To **run Qualimap**, change directories to the `rnaseq` folder and make a `qualimap` folder inside the `results` directory:

``` bash
$ cd ~/rnaseq

$ mkdir -p results/qualimap
```

By default, Qualimap will try to open a GUI to run Qualimap, so we need to run the `unset DISPLAY` command:

``` bash
$ unset DISPLAY
```

We also need to load the qualimap module:

``` bash
$ module load qualimap/2.2.1
```

Now we are ready to run Qualimap on our BAM file! There are different tools or modules available through Qualimap, and the [documentation website](http://qualimap.bioinfo.cipf.es/doc_html/command_line.html) details the tools and options available. We are interested in the `rnaseq` tool. To see the arguments available for this tool we can search the help:

``` bash
$ qualimap rnaseq 
```

We will be running Qualimap with the following specifications:

-   `-outdir`: output directory for html report
-   `-a`: Counting algorithm - uniquely-mapped-reads(default) or proportional (each multi-mapped read is weighted according to the number of mapped locations)
-   `-bam`: path/to/bam/file(s)
-   `-p`: Sequencing library protocol - strand-specific-forward, strand-specific-reverse or non-strand-specific (default)
-   `-gtf`: path/to/gtf/file - **needs to match the genome build and GTF used in alignment**
-   `--java-mem-size=`: set Java memory

``` bash
# assuming you are running the script from /results
qualimap rnaseq -outdir results/qualimap/Mov10_oe_1 -a proportional /
-bam STAR/Mov10_oe_1_Aligned.sortedByCoord.out.bam /
-p strand-specific-reverse /
-gtf /data/NICHD-core0/references/human/gencode-v28/annotation/human_gencode-v28.gtf --java-mem-size=8G
```

### The Qualimap report

The Qualimap report in HTML format should be present in the `results/qualimap` directory. To view this report you need a web browser, so you would need to transfer it over to your laptop. You can do so by Connecting to Server again.

However, you don't need to do that. We generated this report on a subset of data, to get a better idea of the metrics let's **take a look at the report of the full dataset for `Mov_oe_1`**, available [in this zipped folder](https://www.dropbox.com/scl/fi/t3t41hjubytbce6d9g29h/Mov10_oe_1_fulldata_qualimap.zip?rlkey=aum3nw3jym3t4f3l31p6564jn&dl=1). Please download and unzip the folder; find the HTML report within and open it in your browser.

#### **Read alignment summary**

The first few numbers listed in the report are the mapping statistics. Qualimap also computes counts by assigning reads to genes and [reports associated statistics](http://qualimap.conesalab.org/doc_html/analysis.html#id7). For example it computes the following:

-   the number of reads aligned to genes

-   number of ambiguous alignments (reads matching several genes)

-   number of alignments without any feature (intronic and intergenic)

![](../img/qualimap_read_alignment.png){width="700"}

> -   The percentage of mapped reads is a global indicator of the overall sequencing accuracy. We expect between 70-90% of reads to be mapped for the human genome.
> -   Expect a small fraction of reads to be mapping equally well to multiple regions in the genome (‘multi-mapping reads’).
> -   The count related metrics are not as relevant to us since we have quantified with Salmon at the transcript level.

#### **Genomic Origin of Reads: Exon, Intron or Intergenic Region?**

This section reports how many alignments fall into exonic, intronic and intergenic regions along with a number of intronic/intergenic alignments overlapping exons. Exonic region includes 5’UTR,protein coding region and 3’UTR region. This information is summarized in table in addition to a pie chart as shown below.

![](../img/qualimap_genomic_feature.png){width="700"}

> -   Even if you have high genomic mapping rate for all samples, check to see where the reads are mapping. Expect a high proportion of reads mapping to exonic regions (\> 60%) and lower intronic mapping rates (20 -30%).
> -   A higher intronic mapping rate is expected for rRNA removal compared to polyA selection. The intronic reads likely originate from immature transcripts which include either full-length pre-mRNA molecules or nascent transcripts where the RNA polymerase has not yet attached to the 3′ end of the gene.
> -   A roughly equal distribution of reads mapping to intronic, exonic and intergenic regions suggests that there is DNA contamination.
> -   Ribosomal RNA (rRNA) constitutes a large majority of the RNA species in any total RNA preparation. Despite depletion methods, you can never achieve complete rRNA removal. Even with Poly-A enrichment a small percentage of ribosomal RNA can stick to the enrichment beads non-specifically. Excess ribosomal content (\> 2%) will normally have to be filtered out so that differences in rRNA mapped reads across samples do not affect alignment rates and skew subsequent normalization of the data.

#### **Transcript coverage profile**

The profile provides ratios between mean coverage at the 5’ region, the 3’ region and the whole transcript. Coverage plots are generated for all genes total, and also for the 500 highest-expressed and 500 lowest-expressed genes separately.

-   **5’ bias**: the ratio between mean coverage at the 5’ region (first 100bp) and the whole transcript
-   **3’ bias**: is the ratio between mean coverage at the 3’ region (last 100bp) and the whole transcript
-   **5’-3’ bias**: is the ratio between both biases.

<p align="center">

![](../img/qualimap_transcript_coverage.png){width="700"}

</p>

<p align="center">

![](../img/qualimap_coverage_profile.png){width="700"}

</p>

> -   In a perfect sequencing experiment you would expect to see a 5'-3' bias ratio of 1 with low coverage at both ends of the transcript. This would suggest no bias is present.
> -   It is well-documented that libraries prepared with polyA selection can lead to high expression in 3' region (3' bias). At least [one study](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2014-15-6-r86) shows the reverse effect for rRNA removal.
> -   If reads primarily accumulate at the 3’ end of transcripts in poly(A)-selected samples, this might indicate low RNA quality in the starting material.

#### **Junction Analysis**

Qualimap also reports the total number of reads mapping to splice junctions and the 10 most frequent junction rates. The pie chart shows analysis of junction positions in spliced alignments.

-   Known category represents percentage of alignments where both junction sides are known.
-   Partly known represents alignments where only one junction side is known.
-   All other alignments with junctions are marked as Novel.

<p align="center">

![](../img/qualimap_junctions.png){width="700"}

</p>

> Other tools like [RNASeQC](https://software.broadinstitute.org/cancer/cga/rna-seqc) will plot figures that can help evaluate **GC content bias**. This is also an important aspect of QC, as low/high GC content regions will tend to have low coverage.

------------------------------------------------------------------------

## Summary

Taken together, these metrics give us some insight into the quality of our samples and help us in identifying any biases present in our data. The conclusions derived from these QC results may indicate that we need to correct for these biases and so you may want to go back and modify the parameters for Salmon (mapping) accordingly.

------------------------------------------------------------------------

*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*

------------------------------------------------------------------------

*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
