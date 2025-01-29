---
title: "Quality control using FASTQC"
author: Modified by Sally Chang @ NICHD
date: Last modifed in January 2025
duration: 45 minutes
---

## Learning Objectives:

-   Practice using more features of Biowulf (loading modules, running software, parallelization)
-   Describe the contents and format of a FASTQ file
-   Create a quality report using FASTQC
-   Open and look at the report(s) you generated

## Quality Control of FASTQ files

The first step in the RNA-Seq workflow is to take the FASTQ files received from the sequencing facility and assess the quality of the sequence reads.

<p align="center">

<img src="../img/rnaseq_salmon_workflow.png" width="400"/>

</p>

### Unmapped read data (FASTQ)

The [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format) file format is the defacto file format for sequence reads generated from next-generation sequencing technologies. This file format evolved from FASTA in that it contains sequence data, but also contains quality information. Similar to FASTA, the FASTQ file begins with a header line. The difference is that the FASTQ header is denoted by a `@` character. For a single record (sequence read), there are four lines, each of which are described below:

| Line | Description                                                                                               |
|------------|------------------------------------------------------------|
| 1    | Always begins with '\@', followed by information about the read                                           |
| 2    | The actual DNA sequence                                                                                   |
| 3    | Always begins with a '+', and sometimes the same info as in line 1                                        |
| 4    | Has a string of characters representing the quality scores; must have same number of characters as line 2 |

Let's use the following read as an example:

```         
@HWI-ST330:304:H045HADXX:1:1101:1111:61397
CACTTGTAAGGGCAGGCCCCCTTCACCCTCCCGCTCCTGGGGGANNNNNNNNNNANNNCGAGGCCCTGGGGTAGAGGGNNNNNNNNNNNNNNGATCTTGG
+
@?@DDDDDDHHH?GH:?FCBGGB@C?DBEGIIIIAEF;FCGGI#########################################################
```

The line 4 has characters encoding the quality of each nucleotide in the read. The legend below provides the mapping of quality scores (Phred-33) to the quality encoding characters. *Different quality encoding scales exist (differing by offset in the ASCII table), but note the most commonly used one is fastqsanger, which is the scale output by Illumina since mid-2011.*

``` Quality encoding: !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHI``|         |         |         |         | ```

`Quality score:    0........10........20........30........40`

Using the quality encoding character legend, the first nucelotide in the read (C) is called with a quality score of 31 (corresponding to encoding character `@`), and our Ns are called with a score of 2 (corresponding to encoding character `#`). **As you can tell by now, this is a bad read.**

Each quality score represents the probability that the corresponding nucleotide call is incorrect. This quality score is logarithmically based and is calculated as:

```         
Q = -10 x log10(P), where P is the probability that a base call is erroneous
```

These probability values are the results from the base calling algorithm and dependent on how much signal was captured for the base incorporation. The score values can be interpreted as follows:

| Phred Quality Score | Probability of incorrect base call | Base call accuracy |
|:------------------:|:-------------------------------:|:-----------------:|
|         10          |              1 in 10               |        90%         |
|         20          |              1 in 100              |        99%         |
|         30          |             1 in 1000              |       99.9%        |
|         40          |            1 in 10,000             |       99.99%       |

Therefore, for the first nucleotide in the read (C), there is less than a 1 in 1000 chance that the base was called incorrectly. Whereas, for the the end of the read there is greater than 50% probability that the base is called incorrectly.

## **Loading the FASTQC module**

Now that we understand what information is stored in a FASTQ file, the next step is to examine quality metrics for our data.

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) provides a simple way to do some quality checks on raw sequence data coming from high throughput sequencing pipelines. It provides a modular set of analyses, which you can use to obtain an impression of whether your data has any problems that you should be aware of before moving on to the next analysis.

FastQC does the following:

-    accepts FASTQ files (or BAM files) as input

-   generates summary graphs and tables to help assess your data

-   generates an easy-to-view HTML-based report with the graphs and tables

------------------------------------------------------------------------

> NOTE: Before we run FastQC, **you should be on a compute node** in an interactive session. We will start with the default `sinteractive` allocation which is 1 core (2 CPUs) and 768 MB/CPU (1.5 GB) of memory, which should be just fine for our purposes. See this [Biowulf page](https://hpc.nih.gov/docs/userguide.html#int) for more information.
>
> ``` bash
> $ sinteractive 
> ```
>
> ***And we should be in our rnaseq/ directory from within student directories:***
>
> ``` bash
> $ cd /data/Bspc-training/$USER/rnaseq
> ```

### Load the FastQC module

Before we start using software, we have to load the module for each tool. On Biowulf, this is done using an **LMOD** system. It enables users to access software installed on Biwoulf easily, and manages every software's dependencies. The LMOD system adds directory paths of software executables and their dependencies (if any) into the `$PATH` variable.

If we check which modules we currently have loaded, we should not see FastQC.

``` bash
$ module list
```

This is because the FastQC program is not in our \$PATH (i.e. it's not in a directory that shell will automatically check to run commands/programs).

To run the FastQC program, we first need to load the appropriate module, so it puts the program into our path. To find the FastQC module to load we need to search the versions available:

``` bash
$ module spider fastqc
```

Once we know which version we want to use (often the most recent version), we can load the FastQC module.

> **Discussion point**: What are some reasons we might NOT want to use the most updated version of a tool in certain cases?

``` bash
$ module load fastqc/0.12.1
```

Once a module for a tool is loaded, you have essentially made it directly available to you to use like any other basic shell command (for example, `ls`). Now what happens when you list the modules?

``` bash
$ module list
```

**As a reminder - some LMOD commands are listed below**

|            LMOD command            |                                    description                                    |
|:----------------------------------:|:----------------------------------:|
|          `module spider`           |                     List all possible modules on the cluster                      |
|     `module spider modulename`     |                     List all possible versions of that module                     |
|           `module avail`           |                  List available modules available on the cluster                  |
|       `module avail string`        |                   List available modules containing that string                   |
|  `module load modulename/version`  | Add the full path to the tool to `$PATH` (and modify other environment variables) |
|           `module list`            |                                List loaded modules                                |
| `module unload modulename/version` |                             Unload a specific module                              |
|           `module purge`           |                             Unload all loaded modules                             |

## **Running FASTQC on our samples**

Now, let's create a directory to store the output of FastQC inside of the `results` directory you set up last week:

``` bash
$ mkdir results/fastqc
```

We will need to specify this directory in the command to run FastQC. How do we know which argument to use?

``` bash
$ fastqc --help
```

From the help manual, we know that `-o` (or `--outdir`) will create all output files in the specified output directory.

> **You can always check out the Biowulf page for more information about running specific modules in the context of our specific cluster! See the** [fastqc page](https://hpc.nih.gov/apps/fastqc.html) **as an example.**

### **Specifying Input and Output Options**:

**Navigating to where we want to run the program**: To make things simpler, let's navigate to the directory where our input data is so we don't need to use full path names for each of the input files.

```bash
cd raw_data
```

**Specifying output location**: From this directory, how would we specify that we want our output data to end up in the directory we just created? There are a few options based on what we know about full and relative paths. For example: 

```bash
# don't run these, these are just options
fastqc -o ../results/fastqc
fastqc -o /data/Bspc-training/$USER/rnaseq/results/fastqc
```

**Specifying input files**: FastQC will accept multiple file names as input, and we could simply list them individually like so. Note that for FASTQC you don't need to specify an argument before listing input files like we did before specifying the output location.

``` bash
$ cd raw_data
$ fastqc -o ../results/fastqc file1_1.fq file1_2.fq file2_1.fq file2_2.fq
```

 But we could also use our fancy, space-saving `*.fq` wildcard instead. 

``` bash
$ cd raw_data
$ fastqc -o ../results/fastqc *.fq
```

*Did you notice how each file was processed pretty much one at a time?*

## **Using Parallelization** 

FastQC has the capability of splitting up a single process to run on multiple cores! To do this, we will need to:

-   Exit the current interactive session and start another with more additional CPUs , since we started this interactive session with only the default 2 cores. We cannot have a tool to use more cores than requested on a compute node.

-   Specify an argument in FASTQC itself (-t) to process more than one input file at once

Exit the interactive session and start a new one with 6 cores:

``` bash
$ exit  #exit the current interactive session (you will be back on a login node)
$ sinteractive --cpus-per-task=6 --mem=2G
```

Once you are on the compute node, check what job(s) you have running and what resources you are using.

``` bash
$ squeue -u $USER
```

Now that we are in a new interactive session with the appropriate resources, we will need to load the module again for this new session.

``` bash
$ module load fastqc/0.12.1  #reload the module for the new (6-core) interactive session
```

Because we are on a new compute node, `raw_data` directory (remember we are on a new compute node now):

``` bash
$ cd /data/Bspc-training/rnaseq/changes/raw_data
```

Run FastQC and use the multi-threading functionality of FastQC to run 6 jobs at once (with an additional argument `-t`): 

``` bash
$ fastqc -o ../results/fastqc -t 6 *.fq
# note the extra parameter we specified for 6 threads
```

**Discussion Points:**

*Do you notice a difference? Is there anything in the ouput that suggests this is no longer running serially?*
*This overwrote our results. What is a way that we could prevent this from happening?*

## **Viewing results from FastQC**

For each individual FASTQ file that is input to FastQC, there are **two output files that are generated**.

``` bash
$ ls -lh ../results/fastqc/
```

1.  The first is **an HTML file** which is a self-contained document with various graphs embedded into it. Each of the graphs evaluate different quality aspects of our data, we will discuss in more detail in this lesson.
2.  Alongside the HTML file is **a zip file** (with the same name as the HTML file, but with .zip added to the end). This file contains the different plots from the report as separate image files but also contains data files which are designed to be easily parsed to allow for a more detailed and automated evaluation of the raw data on which the QC report is built.

## Viewing the HTML report

We will only need to look at the HTML report for a given input file. It is not possible to view HTML files directly on the cluster from the command line. We will view the HTML result for `Mov10_oe_1.subset.fq` by [locally mounting an HPC System Directory](https://hpc.nih.gov/docs/hpcdrive.html) so you can access the HTMLs locally.

> ### What does this do?
>
> The [HPC System Directories](https://hpc.nih.gov/storage), which include /home, /data, and /scratch, can be mounted to your local workstation if you are on the NIH network or VPN, allowing you to easily drag and drop files between the two places. Note that this is most suitable for transferring small file. Users transferring large amounts of data to and from the HPC systems should continue to use scp/sftp/globus.
>
> Mounting your HPC directories to your local system is particularly userful for viewing HTML reports generated in the course of your analyses on the HPC systems. For these cases, you should be able to navigate to and select the desired html file to open them in your local system's web browser.

Follow the instructions on this Biowulf page for your operating system, and navigate to the `smb://hpcdrive.nih.gov/Bspc-training/username` directory. You will actually need to write out your username here - the `$USER` variable will not work in this context. From here, you can click through to navigate to open `Mov10_oe_1.subset_fastqc.html`. In another lesson this week you will learn more about interpreting this result! 

------------------------------------------------------------------------

*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*

-   *The materials used in this lesson was derived from work that is Copyright © Data Carpentry (<http://datacarpentry.org/>). All Data Carpentry instructional material is made available under the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0).*
