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

[According to NHGRI](https://www.genome.gov/genetics-glossary/Human-Genome-Reference-Sequence), a reference genome (or reference assembly) is an accepted representation of the human genome sequence that is used as a standard for comparison to DNA/RNA sequences generated other studies. Having this standard genome assembly allows researchers to "speak the same language" when it comes to genomic locations and features. While this is a great idea in theory, in practice there is disagreement on chromosome nomenclature across provider.

> **Discussion**: As of February 2025, NCBI hosts [1,831 human genome assemblies](https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=9606). How do we pick a genome to serve as reference for an organism? What qualities make for a good reference genome? Who makes these decisions?

For humans and several other model organisms, the research community makes use of a well-vetted and standardized assembly from the [Reference Genome Consortium](https://www.ncbi.nlm.nih.gov/grc). There are other organism-specific consortia doing annotation and assembly such as [FlyBase](https://flybase.org/) for Drosophila.

As of this writing, the most recent version of this genome is known as `GRCH38.p14`, which is the 38th major "build" (released in 2013) of the assembly, and the 14th update, of this build (released in 2022). These new builds and updates represent corrections and updates to our knowledge of the content and structure of the human genome. This [YouTube video](https://www.youtube.com/watch?v=DeZTPCOKZrg) is a nice introduction to some of the genome assembly nomenclature and compares several versions of the human reference.

**Not all reference genomes are the same.** In general, the sequences are usually identical but very often the chromosome nomenclature will differ. For example, the UCSC Genome Brower usually has chromosomes prefixed by a `chr` (e.g., `chr1`) while Ensembl often does not (e.g., `1`). Even for human, there is a surprising amount of confusion. See [this blog post from 2017 by Heng Li](https://lh3.github.io/2017/11/13/which-human-reference-genome-to-use) (the author of bwa, samtools, seqtk, and others) for the subtleties of different assemblies. Those details matter a lot for variant calling, but for RNA-seq we don't need to be quite as careful.


#### Sources of annotation

Briefly, genome annotation is the process of inferring the identity and location of functional elements, like genes, on an assembled genome. As you can probably imagine, this is crucial for making sense of any sequencing projects!

As is the case for assemblies, major annotation efforts are conducted by research consortia, such as [ENCODE](https://www.encodeproject.org/help/project-overview/), which was started and funded by NHGRI:

> The goal of ENCODE is to build a comprehensive parts list of functional elements in the human genome, including elements that act at the protein and RNA levels, and regulatory elements that control cells and circumstances in which a gene is active. The discovery and annotation of gene elements is accomplished primarily by sequencing a diverse range of RNA sources, comparative genomics, integrative bioinformatic methods, and human curation.

A look at the landing page for the [current GENCODE human genome](https://www.gencodegenes.org/human/) shows us a number of different GTF and FASTA files we could download - *which one is most relevant for mapping with STAR?*

#### Obtaining Reference Genome Files

Luckily for us, the developers of STAR make pretty specific recommendations for us on page 6 of the [STAR manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf):

> GENCODE: files marked with PRI (primary) strongly recommended for mouse and human.

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

**OTHER IMPORTANT NOTE:** For organisms not on GENCODE: Ensembl, UCSC or an organism-specific databases, like FlyBase, are another good choice.

## A note about the GTF Format

**The [GTF](https://web.archive.org/web/20031212200757/http://genes.cse.wustl.edu/GTF2.html) (Gene Transfer Format) file** is a *tab-delimited* file arranged in a very specific manner to store info about gene structure. The ENCODE [data format page](https://www.gencodegenes.org/pages/data_format.html) specifies that the columns of a GTF provided by that consortia are as follows:

| Column Number | Content                                                         | Format                                                      |
|-------------------|---------------------------|---------------------------|
| 1             | chromosome name                                                 | chr1, chr2 etc.                                             |
| 2             | annotation source                                               | ENSEMBL, HAVANA                                             |
| 3             | feature type                                                    | gene, transcript, exon etc.                                 |
| 4             | genomic start location                                          | 1-based integer                                             |
| 5             | genomic end location                                            | integer                                                     |
| 6             | score (not used)                                                | .                                                           |
| 7             | genomic strand                                                  | +, -                                                        |
| 8             | genomic phase                                                   | 0, 1, 2, .                                                  |
| 9             | additional values, stored as key pairs, separated by semicolons. May have spaces and/or quotes | For example: gene_type 'protein_coding'; gene_name "C2CD4C" |

Here is what a few lines of that look like in practice:

```         
chr19   HAVANA   gene   405438   409170   .   -   .   gene_id "ENSG00000183186.7"; gene_type "protein_coding"; gene_name "C2CD4C"; level 2; havana_gene "OTTHUMG00000180534.3";
chr19   HAVANA   transcript   405438   409170   .   -   .   gene_id "ENSG00000183186.7"; transcript_id "ENST00000332235.7"; gene_type "protein_coding"; gene_name "C2CD4C"; transcript_type "protein_coding"; transcript_name "C2CD4C-001"; level 2; protein_id "ENSP00000328677.4"; transcript_support_level "2"; tag "basic"; tag "appris_principal_1"; tag "CCDS"; ccdsid "CCDS45890.1"; havana_gene "OTTHUMG00000180534.3"; havana_transcript "OTTHUMT00000451789.3";
```

> NOTE:
>  [GFF](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md) is a similar format, with some subtle differences. If you can only find a GFF, usually the tool will support it (though you may need an additional flag). Also, many GTF files don't exactly follow the [specification](https://web.archive.org/web/20031212200757/http://genes.cse.wustl.edu/GTF2.html). [This page](https://agat.readthedocs.io/en/latest/gxf.html) gives a well-researched historical overview of the subtle differences and how they have changed over time.

### Refresher about, `grep` `cut` and `sort`

Given our understanding of splice isoforms, we know that a given exon can be part of 2 or more different transcripts generated from the same gene. In a GTF file, this exon will be represented multiple times, once for each transcript (or splice isoform).

``` bash
$ grep "PLEKHN1" /data/Bspc-training/shared/rnaseq_jan2025/human_GRCh38/gencode.v47.primary_assembly.annotation.gtf | head -n 5
```

**`cut` is a command that extracts columns from files.**

We will use `cut` with the `-f` argument to specify which specific fields or columns from the dataset we want to extract. Let's say we want to get the 1st column (chromosome number) and the 4th column (starting genomic position) from `chr1-hg19_genes.gtf` file, we can say:

``` bash
$ cut -f 1,4 /data/Bspc-training/shared/rnaseq_jan2025/human_GRCh38/gencode.v47.primary_assembly.annotation.gtf  | head
```

```         
##description: evidence-based annotation of the human genome (GRCh38), version 47 (Ensembl 113)
##provider: GENCODE
##contact: gencode-help@ebi.ac.uk
##format: gtf
##date: 2024-07-19
chr1    11121
chr1    11121
chr1    11121
chr1    12010
chr1    12613
```

> The `cut` command assumes our data columns are separated by tabs (i.e. tab-delimited). Since GTFs are a tab-delimited file, so the default `cut` command works for us. However, data can be separated by other types of delimiters like "," or ";". If your data is not tab delimited, there is an argument you can add to your `cut` command, `-d` to specify the delimiter (e.g. `-d ","` with a .csv file).

The output of `cut` doesn't seem very useful yet as you still have the header lines (that start with \##). How do we disregard those?

#### Grep -v : Reverse grep

The `-v` option is shorthand for **`--invert-match`**, which selects all the lines that DO NOT have matches to your pattern. We can chain it together with the `cut` command using a pipe:

``` bash
$ grep -v /data/Bspc-training/shared/rnaseq_jan2025/human_GRCh38/gencode.v47.primary_assembly.annotation.gtf | cut -f 1,4 | head -n 5
```

```         
chr1	11121
chr1	11121
chr1	11121
chr1	12010
chr1	12613
```

#### Sort 

**`sort` is a command used to sort the contents of a file in a particular order.** It has arguments that let you pick which column to sort by (`-k`), what kind of sorting you want to do (numeric `n`) and also if the result of the sorting should only return unique (`-u`) values. These are just 2 of the many features of the sort command.

Let's do a quick test of how the `-u` argument returns only unique lines (and remove duplicates).

``` bash
$ grep "##" -v /data/Bspc-training/shared/rnaseq_jan2025/human_GRCh38/gencode.v47.primary_assembly.annotation.gtf | cut -f 1,4 | wc -l 
```

*How many lines are returned to you?*

<details>

<summary><b><i>Click here to check your output</i></b></summary>

<p>Your command should have returned 4117647 lines</p>
</details>

Now apply the `sort -u` command before counting the lines.

``` bash
$ grep "##" -v /data/Bspc-training/shared/rnaseq_jan2025/human_GRCh38/gencode.v47.primary_assembly.annotation.gtf | cut -f 1,4 | sort -u | wc -l 
```

*How many lines do you see now?*

<details>

<summary><b><i>Click here to check your output</i></b></summary>

<p>Your command should have returned 796,715 lines</p>

</details>

You can practice chaining together these skills in the GTF Exercise in Week 1 Lesson 04 if you haven't done so already. The assignment for this week also challenges you to put a few of these together.

------------------------------------------------------------------------

## Create a genome index with STAR

We will be using STAR for aligning reads to the genome. Aligners like STAR need an **index** in order to align reads efficiently. An index is a computationally-efficient data structure created from a FASTA file -- it turns out to be unreasonably computationally expensive to try to match reads to the text FASTA file, so this preparation step is required. It only has to be done once for each FASTA file.

Before we set up a script, let's explore the possible versions of STAR on Biowulf.

``` bash
$ module spider star
```

```         
  Versions:
        STAR/2.5.4a
        STAR/2.6.1c
        STAR/2.7.0f
        STAR/2.7.3a
        STAR/2.7.6a
        STAR/2.7.8a
        STAR/2.7.9a
        STAR/2.7.10b
        STAR/2.7.11b
```

Let's go ahead and load the default version of STAR, which happens to be the most recent:

```         
$ module load STAR
$ [+] Loading STAR  2.7.11b
```

> **IMPORTANT NOTE**: This is one of those times where you may not want to always use the most recent version. Newer versions of STAR are not backwards compatible with all genome indices created by older versions. Something to keep in mind if you are trying to re-work old analyses! We'll look at available indices on Biowulf below.

The basic options to **generate genome indices** using STAR are as follows:

-   `--runThreadN`: number of threads, should match the number of CPUS we specify in the SLURM directives
-   `--runMode`: genomeGenerate mode
-   `--genomeDir`: /path/to/store/genome_indices
-   `--genomeFastaFiles`: /path/to/FASTA_file
-   `--sjdbGTFfile`: /path/to/GTF_file
-   `--sjdbOverhang`: readlength -1

> *NOTE:* In case of reads of varying length, the ideal value for `--sjdbOverhang` is max(ReadLength)-1. In most cases, the default value of 100 will work similarly to the ideal value.

Now let's create a job submission script to generate the genome index. *We will NOT run this script, as it uses a lot of memory, will take up a bunch of disk space and will take a long time to run! This is the script I used to generate the files we will use in the next lesson.*

``` bash
# Since we are not actually going to run this, you can work on this text document from any directory
$ vim genome_index.sh
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
cd /data/Bspc-training/shared/rnaseq_jan2025/human_GRCh38

STAR --runThreadN 8 --runMode genomeGenerate --genomeDir /data/Bspc-training/shared/rnaseq_jan2025/human_GRCh38 --genomeFastaFiles /data/Bspc-training/shared/rnaseq_jan2025/human_GRCh38/GRCh38.primary_assembly.genome.fa --sjdbGTFfile /data/Bspc-training/shared/rnaseq_jan2025/human_GRCh38/gencode.v47.primary_assembly.annotation.gtf --sjdbOverhang 99
```

#### Results

For this workshop we have generated the genome indices for you, so that we don't get held up waiting on the generation of the indices (it takes a while and requires a lot of memory). To see what all the genome index result files are, take a look:

``` bash
$ ls -lh /data/Bspc-training/shared/rnaseq_jan2025/human_GRCh38
$ -rw-rw----+ 1 changes Bspc-training 1.2K Feb  6 11:32 chrLength.txt
-rw-rw----+ 1 changes Bspc-training 3.2K Feb  6 11:32 chrNameLength.txt
-rw-rw----+ 1 changes Bspc-training 2.0K Feb  6 11:32 chrName.txt
-rw-rw----+ 1 changes Bspc-training 2.1K Feb  6 11:32 chrStart.txt
-rw-rw----+ 1 changes Bspc-training  74M Feb  6 11:32 exonGeTrInfo.tab
-rw-rw----+ 1 changes Bspc-training  30M Feb  6 11:32 exonInfo.tab
-rw-rw----+ 1 changes Bspc-training 1.7G Feb  6 11:32 gencode.v47.primary_assembly.annotation.gtf #input
-rw-rw----+ 1 changes Bspc-training 3.2M Feb  6 11:32 geneInfo.tab
-rw-rw----+ 1 changes Bspc-training 3.1G Feb  6 11:32 Genome
-rw-rw----+ 1 changes Bspc-training  872 Feb  6 11:32 genomeParameters.txt
-rw-rw----+ 1 changes Bspc-training 3.0G Feb  6 11:32 GRCh38.primary_assembly.genome.fa #input
-rw-rw----+ 1 changes Bspc-training  33K Feb  6 11:32 Log.out
-rw-rw----+ 1 changes Bspc-training  24G Feb  6 11:32 SA
-rw-rw----+ 1 changes Bspc-training 1.5G Feb  6 11:32 SAindex
-rw-rw----+ 1 changes Bspc-training  15M Feb  6 11:32 sjdbInfo.txt
-rw-rw----+ 1 changes Bspc-training  17M Feb  6 11:32 sjdbList.fromGTF.out.tab
-rw-rw----+ 1 changes Bspc-training  14M Feb  6 11:32 sjdbList.out.tab
-rw-rw----+ 1 changes Bspc-training  25M Feb  6 11:32 transcriptInfo.tab
```

## Genome Indices on Biowulf

Going right to the source to download the GTF and FASTA is a great way of making sure you are using the most up-to-date version.

However, Biowulf also provides centrally-maintained [scientific reference databases](https://hpc.nih.gov/refdb/index.php) for users, which includes many pre-generated genome indices for various alignment software.

According to that page, STAR indices are all kept in a directory: `/fdb/STAR_indices`, which contains the following files and subdirectories:

```         
00init.sh  01build.sh       03check_repo.sh  2.6.1c  2.7.10b  2.7.3a  2.7.8a  refdb.yml 00lib.sh   02dedup_existing.sh  2.5.4 2.7.0f  2.7.11b  2.7.6a  2.7.9a  repo
```

Navigating further into 2.7.11b (indices built using the most recent version of STAR on Biowulf), we can eventually get to the ENCODE human genome versions:

``` bash
$ ls -lh /fdb/STAR_indices/2.7.11b/GENCODE/Gencode_human/
drwxrwxr-x 2 wresch staff 4.0K Mar 14  2024 release_19
drwxrwxr-x 2 wresch staff 4.0K Mar 14  2024 release_39
drwxrwxr-x 2 wresch staff 4.0K Mar 14  2024 release_45
```

**Looking into this `release_45` we see a few important things:**
- A number of subdirectories with names like `genes-100/`. These are the actual genome index directories, categorized by read length! Remember how we had to set that `--sjdbOverhang 99` parameter? 
- `genes.gtf` and other files used to actually create the indicies
- `slurm-21952863.out` and similar. These are logs from when Biowulf staff created these indices. You can check out what actual command they used. 

**Generally speaking, STAR indices for each version of STAR available on Biowulf are stored in:**

-   **`/fdb/STAR_indices/[STAR VERSION]`**

-   **`/fdb/STAR_current`**`→ /fdb/STAR_indices/[current default STAR version]`

Keep in mind the lack of backwards compatibility between some versions of STAR! 

## Exercise:

Note that the most recent Gencode annotations they have are for `release_45`, whereas the up-to-date version we used earlier are from Release 47.

-   Does it look like there is a `STAR 2.7.11b` index available for a recent release on Biowulf somewhere in `/fdb/STAR_indices`? You may need to check both in both the Gencode and UCSC subfolders.
-   If not - look in the `/fdb/STAR_indices/2.7.10b` directory or other, older directories. This slightly older version of STAR has more pre-prepared indices. If you find a relevant index folder - report the full path to that directory.

## Assignment

By looking at recent publications, talking to your labmates etc. - what is the most recent major genome build for your organism used by your research community? For example `GRCh38` for human - sometimes abbreviated to `hg38` by providers such as UCSC.

1.  Create a `reference_genome` subdirectory of your `/Bspc-training/$USER/rnaseq` directory. Create a `ref_notes.txt` file in there to answer the next questions.
2.  Find a **GTF file**, in a Biowulf STAR index directory, that corresponds with your chosen genome version. OR if necessary, use `wget` to download the GTF file you identify from a consortium website. *Post in the main `rnaseq_jan2025` channel if you need help with this!*
3.  In the text file, make note of: the name of your organism, the genone build you would use, and the full Biowulf directory paths or `wget` commands you used to download the GTF and FASTA files.

## BONUS PRACTICE ASSIGNMENT:

A few Fridays ago, we demonstrated how to create a frequency table of different **feature types** from a GTF file using a combination of commands like `cut`, `grep` and `uniq`. The result ended up looking like this:

```         
gene  78724
transcript 4038715
exon 3651434
CDS 1872631
UTR 664938
start_codon 99382
stop_codon 93628
Selenocysteine 130
```

Using a combination of commands like `cut`, `grep` and `uniq` - can you figure out how to create a table like this for the GTF file you found? **Store the results in another text file along with the command you used!**

Hint: The end of [Week 1 Lesson 04](https://nichd-bspc.github.io/intro-rnaseq-hpc/lessons/wk1_lesson04_searching_files.html) gives some great examples of how to build more complicated commands from those individual BASH commands.
