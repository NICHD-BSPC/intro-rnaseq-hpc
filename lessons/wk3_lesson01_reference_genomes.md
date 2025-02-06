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

As mentioned above, we need to identify the relevant reference genome and an associated annotation file.

#### Reference Consortia

[According to NHGRI](https://www.genome.gov/genetics-glossary/Human-Genome-Reference-Sequence), a reference genome (or reference assembly) is an accepted representation of the human genome sequence that is used as a standard for comparison to DNA/RNA sequences generated other studies. Having this standard genome assembly allows researchers to "speak the same language" when it comes to genomic locations and features.

> **Discussion**: As of February 2025, NCBI hosts [1,831 human genome assemblies](https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=9606). How do we pick a genome to serve as reference for an organism? What qualities make for a good reference genome? Who makes these decisions?

For humans and several other model organisms, the research community makes use of a well-vetted and standardized assembly from the [Reference Genome Consortium](https://www.ncbi.nlm.nih.gov/grc). There are other organism-specific consortia doing annotation and assembly such as [FlyBase](https://flybase.org/) for Drosophila.

As of this writing, the most recent version of this genome is known as `GRCH38.p14`, which is the 38th major "build" (released in 2013) of the assembly, and the 14th update, of this build (released in 2022). These new builds and updates represent corrections and updates to our knowledge of the content and structure of the human genome. This [YouTube video](https://www.youtube.com/watch?v=DeZTPCOKZrg) is a nice introduction to some of the genome assembly nomenclature and compares several versions of the human reference.

#### Sources of annotation 

Briefly, genome annotation is the process of inferring the identity and location of functional elements, like genes, on an assembled genome. As you can probably imagine, this is crucial for making sense of any sequencing projects!

As is the case for assemblies, major annotation efforts are conducted by research consortia, such as [ENCODE](https://www.encodeproject.org/help/project-overview/), which was started and funded by NHGRI:

> The goal of ENCODE is to build a comprehensive parts list of functional elements in the human genome, including elements that act at the protein and RNA levels, and regulatory elements that control cells and circumstances in which a gene is active. The discovery and annotation of gene elements is accomplished primarily by sequencing a diverse range of RNA sources, comparative genomics, integrative bioinformatic methods, and human curation.

A look at the landing page for the [current GENCODE human genome](https://www.gencodegenes.org/human/) shows us a number of different GTF and FASTA files we could download - *which one is most relevant for mapping with STAR?*

#### Obtaining Reference Genome Files

Luckily for us, the developers of STAR make pretty specific recommendations for us on page 6 of the [STAR manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf):

> GENCODE: files marked with PRI (primary) strongly recommended for mouse and human

The PRI files contain the comprehensive gene annotation on the primary assembly (chromosomes and scaffolds) sequence region but not any alternative loci haplotypes, which STAR cannot make use of during mapping.

So, from that page, find the following listings and copy the links for the `PRI` GTF file and the FASTA for `Genome sequence, primary assembly`.

Once we have those URLS handy, we use `wget` , which is a built-in piece of software that downloads files stored at HTTP, HTTPS, FTP and FTPS addresses. These are the commands I used to download these files into our shared space:

``` bash
# DO NOT RUN - I have already downloaded this file for us
$ wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/GRCh38.primary_assembly.genome.fa.gz
```

``` bash
# DO NOT RUN - I have already downloaded this file for us
$ wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.primary_assembly.annotation.gtf.gz
```

**IMPORTANT NOTE:** Due to differences in formatting and gene/chromosome names between versions and consortia, it is recommended you use GTFs and FASTAs from the SAME version and SAME provider (i.e. ENCODE) for analyses.

## A note about the GTF Format

**The GTF (Gene Transfer Format) file** is a *tab-delimited* file arranged in a very specific manner to store info about gene structure. The ENCODE [data format page](https://www.gencodegenes.org/pages/data_format.html) specifies that the columns of a GTF provided by that consortia are as follows:

| Column Number | Content                                                         | Format                                                      |
|--------------|-----------------------------|-----------------------------|
| 1             | chromosome name                                                 | chr1, chr2 etc.                                             |
| 2             | annotation source                                               | ENSEMBL, HAVANA                                             |
| 3             | feature type                                                    | gene, transcript, exon etc.                                 |
| 4             | genomic start location                                          | 1-based integer                                             |
| 5             | genomic end location                                            | integer                                                     |
| 6             | score (not used)                                                | .                                                           |
| 7             | genomic strand                                                  | +, -                                                        |
| 8             | genomic phase                                                   | 0, 1, 2, .                                                  |
| 9             | additional values, stored as key pairs, separated by semicolons | For example: gene_type 'protein_coding'; gene_name "C2CD4C" |

Here is what a few lines of that look like in practice:

```         
chr19   HAVANA   gene   405438   409170   .   -   .   gene_id "ENSG00000183186.7"; gene_type "protein_coding"; gene_name "C2CD4C"; level 2; havana_gene "OTTHUMG00000180534.3";
chr19   HAVANA   transcript   405438   409170   .   -   .   gene_id "ENSG00000183186.7"; transcript_id "ENST00000332235.7"; gene_type "protein_coding"; gene_name "C2CD4C"; transcript_type "protein_coding"; transcript_name "C2CD4C-001"; level 2; protein_id "ENSP00000328677.4"; transcript_support_level "2"; tag "basic"; tag "appris_principal_1"; tag "CCDS"; ccdsid "CCDS45890.1"; havana_gene "OTTHUMG00000180534.3"; havana_transcript "OTTHUMT00000451789.3";
```

Let's see what kind of features are annotated for the study's gene of interest Mov10 Helicase (`MOV10`), which is on Chromosome 1.

``` bash
$ grep "MOV10" gencode.v47.primary_assembly.annotation.gtf | wc -l
$ 772 # Does this seem like a biologically reasonable number? 
```

Taking a peek at the results using `less` or something like `grep "MOV10" gencode.v47.primary_assembly.annotation.gtf | tail -n 5` , we see some results for a different gene with a similar name: `MOV101L` on `chr22`.

How do we

## Create a genome index with STAR

Before we set up a script, let's explore the possible versions of STAR are on Biowulf.

``` bash
$ module spider star
```

The basic options to **generate genome indices** using STAR are as follows:

-   `--runThreadN`: number of threads, should match the number of CPUS we specify in the SLURM directives
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

#SBATCH --partition=quick #quick partition
#SBATCH --time=04:00:00 # time limit
#SBATCH --cpus-per-task=8 # number of cores
#SBATCH --mem=48g # requested memory
#SBATCH --job-name grch38_star_index # Job name
#SBATCH -o %j.out # File to which standard output will be written
#SBATCH -e %j.err # File to which standard error will be written
#SBATCH --mail-type=BEGIN,END

# load most recent STAR module 
module load STAR

# Change directory to where the data is
cd /data/changes/reference_genomes/human_GRCh38

STAR --runThreadN 8 --runMode genomeGenerate --genomeDir /data/changes/reference_genomes/human_GRCh38 --genomeFastaFiles /data/changes/reference_genomes/GRCh38.primary_assembly.genome.fa --sjdbGTFfile /data/changes/reference_genomes/gencode.v47.primary_assembly.annotation.gtf --sjdbOverhang 99
```

``` bash
$ sbatch ~/rnaseq/scripts/genome_index.run
```

For this workshop we have generated the genome indices for you, so that we don't get held up waiting on the generation of the indices (it takes a while and requires a lot of memory). The index can be found in `/data/NICHD-core0/references/human/gencode-v28/genome/star/human_gencode-v28`

## Genome Indices on Biowulf

Going right to the source to download the GTF Files is

A quick note on shared databases for human and other commonly used model organisms. The O2 cluster has a designated directory at `/n/groups/shared_databases/` in which there are files that can be accessed by any user. These files contain, but are not limited to, genome indices for various tools, reference sequences, tool specific data, and data from public databases, such as NCBI and PDB. So when using a tool that requires a reference of sorts, it is worth taking a quick look here because chances are it's already been taken care of for you.

> ``` bash
> $ ls -l /n/groups/shared_databases/igenome/
> ```

## Assignment 
