---
title: "Reference genomes and genome indices"
author: "Harvard HPC Staff, modified by Sally Chang @ NICHD"
date: Last modified February 2025
---

Approximate time: 50 minutes

## Learning Objectives:

-   Download reference genome files needed for mapping
-   Find an appropriate reference genome for "your" organism
-   Build a genome index that will be used in the next lesson
-   Use Bash commands to summarize aspects of a GTF file

## Intro to Read Alignment

<img src="../img/RNAseqWorkflow.png" width="400"/>

Now that we have explored the quality of our raw reads, we can move on to read alignment. We perform read alignment or mapping to determine where in the genome the reads originated from. The alignment process consists of choosing an appropriate reference genome to map our reads against and performing the read alignment using one of several splice-aware alignment tools. We are going to make use of the [STAR](http://bioinformatics.oxfordjournals.org/content/early/2012/10/25/bioinformatics.bts635) aligner.

**No matter what aligner you are using steps involved are:**

1.  Creating a genome index (after identifying the appropriate reference for your organism or use case). For STAR, the necessary input is:
    1.  The reference sequence in FASTA format
    2.  Annotations associated with the reference genome in GTF format
2.  Mapping reads to the genome, which needs:
    1.  Your quality-checked reads in FASTQ format
    2.  The output of the genome indexing step, above

## Prepping for this lesson

To get started with this lesson, start an interactive session:

``` bash
$ sinteractive
```

You should have a directory tree setup similar to that shown below. It is best practice to have all files you intend on using for your workflow present within the same directory.

``` bash
rnaseq
    ├── logs
    ├── meta
    ├── raw_data
    │   ├── Irrel_kd_1.subset.fq
    │   ├── Irrel_kd_2.subset.fq
    │   ├── Irrel_kd_3.subset.fq
    │   ├── Mov10_oe_1.subset.fq
    │   ├── Mov10_oe_2.subset.fq
    │   └── Mov10_oe_3.subset.fq
    ├── results
    └── scripts
```

## Finding Reference Genomes

As mentioned above, we need to identify the relevant reference genome and an associated annotation file. [According to NHGRI](https://www.genome.gov/genetics-glossary/Human-Genome-Reference-Sequence), a reference genome (or reference assembly) is an accepted representation of the human genome sequence that is used as a standard for comparison to DNA/RNA sequences generated other studies. Having this standard genome assembly allows authors to

> **Discussion**: As of February 2025, NCBI hosts [1,831 human genome assemblies](https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=9606). How do we pick a genome to serve as reference for an organism? What qualities make for a good reference genome? Who makes these decisions?

For humans and several other model organisms, the authority for

How do we know which one of these assemblies and annotation files is right for my organism?

Once we have identified files we are interested in, we use `wget` , which is a built-in piece of software that downloads files stored at HTTP, HTTPS, FTP and FTPS addresses.

``` bash
# DO NOT RUN - I have already downloaded this file for us
$ wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/GRCh38.primary_assembly.genome.fa.gz
```

``` bash
# DO NOT RUN - I have already downloaded this file for us
$ wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.primary_assembly.annotation.gtf.gz
```

It contains the comprehensive gene annotation on the primary assembly (chromosomes and scaffolds) sequence region but not any alternative loci haplotypes.

## The GTF Format

**Exercises**

Now that we know what type of information is inside the GTF file, let's use the commands we have learned so far to answer a simple question about our data: **how many unique exons are present on chromosome 1 using `chr1-hg19_genes.gtf`?**

To determine the number of unique exons on chromosome 1, we are going to perform a series of steps as shown below. In this exercise, you need to figure out the command line for each step.

1.  Extract only the genomic coordinates of exon features
2.  Subset dataset to only keep genomic coordinates
3.  Remove duplicate exons
4.  Count the total number of exons

Your end goal is to have a single line of code, wherein you have strung together multiple commands using the pipe operator. But, we recommend that you do it in a stepwise manner as detailed below.

#### 1. Extract only the genomic coordinates of exon features

We only want the exons (not CDS or start_codon features), so let's use `grep` to search for the word "exon". You should do sanity check on the first few lines of the output of `grep` by piping the result to the `head` command. ***Report the command you have at this stage.***

#### 2. Subset the extracted information from step 1 to only keep genomic coordinates

We will define the uniqueness of an exon by its genomic coordinates, both start and end. Therefore, from the step 1 output, we need to keep 4 columns (chr, start, stop, and strand) to find the total number of unique exons. The column numbers you want are 1, 4, 5, and 7.

You can use `cut` to extract those columns from the output of step 1. ***Report the command you have at this stage.***

At this point, the first few lines should look like this:

```         
chr1    14362   14829   -
chr1    14970   15038   -
chr1    15796   15947   -
chr1    16607   16765   -
chr1    16858   17055   -
```

#### 3. Remove duplicate exons

Now, we need to remove those exons that show up multiple times for different transcripts. We can use the `sort` command with the `-u` option. ***Report the command you have at this stage.***

Do you see a change in how the sorting has changed? By default the `sort` command will sort and what you can't see here is that it has removed the duplicates. We will use step 4 to check if this step worked.

#### 4. Count the total number of exons

First, check how many lines we would have without using `sort -u` by piping the output to `wc -l`.

Now, to count how many unique exons are on chromosome 1, we will add back the `sort -u` and pipe the output to `wc -l`. Do you observe a difference in number of lines?

***Report the command you have at this stage and the number of lines you see with and without the `sort -u`.***

<details>

<summary><b><i>Answers</i></b></summary>

<p><i>Question 1</i><br> <code>grep exon chr1-hg19_genes.gtf \| head</code><br></p>

<p><i>Question 2</i><br> <code>grep exon chr1-hg19_genes.gtf \| cut -f 1,4,5,7 \| head</code><br></p>

<p><i>Question 3</i><br> <code>grep exon chr1-hg19_genes.gtf \| cut -f 1,4,5,7 \| sort -u \| head</code><br></p>

<p><i>Question 4</i><br> <code>grep exon chr1-hg19_genes.gtf \| cut -f 1,4,5,7 \| wc -l</code><br> The output returns 37,213 lines.<br> <code>grep exon chr1-hg19_genes.gtf \| cut -f 1,4,5,7 \| sort -u \| wc -l</code><br> The output returns 22,769 lines, indicating that repetitive lines have been removed.<br></p>

</details>

## Create a genome index with STAR

To use the STAR aligner, load the module:

``` bash
$ module load star/2.7.6a 
```

For this workshop we have generated the genome indices for you, so that we don't get held up waiting on the generation of the indices (it takes a while and requires a lot of memory). The index can be found in `/data/NICHD-core0/references/human/gencode-v28/genome/star/human_gencode-v28`

The command to create an index can be found in the job submission script we have linked [here](../scripts/star_genome_index.run)

> **NOTE:** By default the latest human genome build, GRCh38, contains information about alternative alleles for various locations on the genome. If using this version of the GRCh38 genome then it is advisable to use the HISAT2 aligner as it is able to utilize this information during the alignment. There is a version of GRCh38 available that does not have these alleles represented, which is the appropriate version to use with STAR. This is because STAR does not have the functionality to appropriately deal with the presence of alternate alleles as yet.

The basic options to **generate genome indices** using STAR are as follows:

-   `--runThreadN`: number of threads
-   `--runMode`: genomeGenerate mode
-   `--genomeDir`: /path/to/store/genome_indices
-   `--genomeFastaFiles`: /path/to/FASTA_file
-   `--sjdbGTFfile`: /path/to/GTF_file
-   `--sjdbOverhang`: readlength -1

> *NOTE:* In case of reads of varying length, the ideal value for `--sjdbOverhang` is max(ReadLength)-1. In most cases, the default value of 100 will work similarly to the ideal value.

Now let's create a job submission script to generate the genome index:

``` bash
$ vim ~/rnaseq/scripts/genome_index.run
```

Within `vim` we now add our shebang line, the SLURM directives, and our STAR command.

``` bash
#!/bin/bash

#SBATCH -p short        # partition name
#SBATCH -t 0-2:00       # hours:minutes runlimit after which job will be killed
#SBATCH -n 6        # number of cores requested -- this needs to be greater than or equal to the number of cores you plan to use to run your job
#SBATCH --mem 16G
#SBATCH --job-name STAR_index       # Job name
#SBATCH -o %j.out           # File to which standard out will be written
#SBATCH -e %j.err       # File to which standard err will be written

cd /n/scratch2/username/

module load gcc/6.2.0 star/2.5.2b

STAR --runThreadN 6 \
--runMode genomeGenerate \
--genomeDir chr1_hg38_index \
--genomeFastaFiles /n/groups/hbctraining/intro_rnaseq_hpc/reference_data_ensembl38/Homo_sapiens.GRCh38.dna.chromosome.1.fa \
--sjdbGTFfile /n/groups/hbctraining/intro_rnaseq_hpc/reference_data_ensembl38/Homo_sapiens.GRCh38.92.gtf \
--sjdbOverhang 99
```

``` bash
$ sbatch ~/rnaseq/scripts/genome_index.run
```

## Reference Genomes on Biowulf

> ***UPDATE BASED ON SHARED LOCATION FOR BSPC AND/OR BIOWULF***
>
> A quick note on shared databases for human and other commonly used model organisms. The O2 cluster has a designated directory at `/n/groups/shared_databases/` in which there are files that can be accessed by any user. These files contain, but are not limited to, genome indices for various tools, reference sequences, tool specific data, and data from public databases, such as NCBI and PDB. So when using a tool that requires a reference of sorts, it is worth taking a quick look here because chances are it's already been taken care of for you.
>
> ``` bash
> $ ls -l /n/groups/shared_databases/igenome/
> ```
