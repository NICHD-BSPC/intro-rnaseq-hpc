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

## Reference Genomes

As mentioned above, we will need to find both a

> ***UPDATE BASED ON SHARED LOCATION FOR BSPC AND/OR BIOWULF***
>
> A quick note on shared databases for human and other commonly used model organisms. The O2 cluster has a designated directory at `/n/groups/shared_databases/` in which there are files that can be accessed by any user. These files contain, but are not limited to, genome indices for various tools, reference sequences, tool specific data, and data from public databases, such as NCBI and PDB. So when using a tool that requires a reference of sorts, it is worth taking a quick look here because chances are it's already been taken care of for you.
>
> ``` bash
> $ ls -l /n/groups/shared_databases/igenome/
> ```

``` bash
$ wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/GRCh38.primary_assembly.genome.fa.gz
```

``` bash
$ wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.primary_assembly.annotation.gtf.gz
```

It contains the comprehensive gene annotation on the primary assembly (chromosomes and scaffolds) sequence region but not any alternative loci haplotypes.

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

#### 
