---
title: "Automating the RNA-seq workflow"
author: "Harvard HPC Staff, Sally Chang @ NICHD"
date: "Last Edited February 2025"
---

## Learning Objectives:

-   Create a reusable and efficient workflow for RNA-seq data analysis using shell scripts

## Automating the analysis path from Sequence reads to Count matrix

Once you have optimized all the tools and parameters using a single sample (likely using an interactive session), you can write a script to run the whole workflow on all the samples in parallel.

This will ensure that you run every sample with the exact same parameters, and will enable you to keep track of all the tools and their versions. In addition, the script is like a lab notebook; in the future, you (or your colleagues) can go back and check the workflow for methods, which enables efficiency and reproducibility.

### Start an interactive session

We will be working with an interactive session has enough memory, CPUs and scratch space to run all of the commands we have run in the whole course up to this point!

Remember from Week 3 Lesson 02 that the aligning step using STAR seemed to be the most memory intensive and required scratch space to be designated:

``` bash
$ sinteractive --cpus-per-task=12 --mem=48g --gres=lscratch:1
```

### More Flexibility with variables

We can write a shell script that will run on a specific file, but to make it more flexible and efficient we would prefer that it lets us give it an input fastq file when we run the script. To be able to provide an input to any shell script, we need to use **Positional Parameters**.

For example, we can refer to the components of the following command as numbered variables **within** the actual script:

``` bash
# * DO NOT RUN *
sh  run_rnaseq.sh  input.fq  input.gtf  12
```

`$0` =\> run_rnaseq.sh

`$1` =\> input.fq

`$2` =\> input.gtf

`$3` =\> 12

The variables `$1`, `$2`, `$3`,...`$9` and so on are **positional parameters** in the context of the shell script, and can be used within the script to refer to the files/number specified on the command line. Basically, the script is written with the expectation that `$1` will be a fastq file and `$2` will be a GTF file, and so on.

Note that `$1`, which you may have seen before, is actually a short form of `${1}` and we can only use `$1` when it is **not** followed by a letter, digit or an underscore but we can always use `${1}`. Using `${1}` is best practice and what we will use for the rest of this lesson.

*There can be virtually unlimited numbers of inputs to a shell script, but it is wise to only have a few inputs to avoid errors and confusion when running a script that used positional parameters.*

> [This is an example of a simple script that used the concept of positional parameters and the associated variables](http://steve-parker.org/sh/eg/var3.sh.txt). You should try this script out after the class to get a better handle on positional parameters for shell scripting. You can also learn more about positional parameters [here](https://hbctraining.github.io/Training-modules/Intermediate_shell/lessons/positional_params.html)

We will be using this concept in our automation script, wherein we will accept the full or relative path to a file as input.

## Writing the automation script: Setting up variables

We're about to do substantial editing of a file. Use vim if you're comfortable with it (e.g., `vim run_rnaseq.sh`); otherwise use another text editor on your laptop and you can paste it into vim later.

Let's begin with the shebang line and a `cd` command to make paths more convenient to write:

``` bash
#!/bin/bash/

# change directories to your course directory so that all the analysis is
# stored there and we can write all paths relative to this directory.

cd /data/Bspc-training/$USER/rnaseq
```

**We want users to input the path to the fastq file as input to the shell script.** To make this happen, we will use the `${1}` positional parameter variable within the script.

Since `${1}` will store the path to the fastq file, including the file name, we will be referring to it every time we need to specify the fastq file in any commands. We could just use the variable `${1}`, but that is not an intuitive variable name for a fastq file, is it? So we want to create a new variable called `fq` and copy the contents of `${1}` into it.

``` bash
# we're expecting the first argument to be a fastq file, so copy it to a more
# clearly-named variable

fq=${1}
```

In the rest of the script, we can now call the fastq file using `${fq}` instead of `${1}`, which will make our code easier to understand and debug.

> When we set up variables we do not use the `$` before it, but when we *use the variable*, we always have to have the `$` before it. \>
>
> For example:
>
> initializing the `fq` variable =\> `fq=${1}`
>
> using the `fq` variable =\> `fastqc ${fq}`

Next, we want to extract the name of the sample from `${fq}` which contains the full name of the file and possibly the path to the file as well. The reason to extract the sample name is so that all the output files from this workflow are appropriately named with sample identifier.

We can obtain the sample name by using the `basename` command on the `${fq}` (or `${1}`) variable, and save it in a new variable called `samplename`.

``` bash
# grab base of filename for naming outputs for your script PUT THIS IN YOUR SCRIPT

samplename=$(basename ${fq} .subset.fq)
echo "Sample name is ${samplename}"           
```

> **Intro to `basename`**
>
> 1. the `basename` command: this command takes a path or a name and trims away all the information before the last `/`. If you also specify the string to clear away at the end, it will do that as well. In this case, if the variable `${fq}` contains the path `*/data/Bspc-training/changes/rnaseq/raw_data/Mov10_oe_1.subset.fq*`, then `basename ${fq} .subset.fq` will output `Mov10_oe_1`.
> 2. We encapsulate the `basename...` command in `$(...)` which is called command substitution. It places the results of the command into the variable.

``` bash
# basename demo. Run on the command line to inspect; do not put in your script
$ fq=/data/Bspc-training/changes/rnaseq/raw_data/Mov10_oe_1.subset.fq
$ echo $fq
$ samplename=$(basename ${fq} .subset.fq)
$ echo "Sample name is ${samplename}"
```

Next we want to specify how many cores the script should use to run the analysis. This provides us with an easy way to modify the script to run with more or fewer cores without have to replace the number within all commands where cores are specified.

Here, we take advantage of the fact that [Slurm creates environment variables](https://slurm.schedmd.com/sbatch.html#SECTION_OUTPUT-ENVIRONMENT-VARIABLES) for us when a job starts; we'll use the one that indicates how many CPUs we asked for.

``` bash
# specify the number of cores to use. For an initial value, match the number to the interactive node we just started.

cores=${SLURM_CPUS_PER_TASK}
```

Next we'll initialize 2 more variables named `genome` and `gtf`, these will contain the paths to where the reference files are stored. This makes it easier to modify the script for when you want to use a different genome, i.e. you'll just have to change the contents of these variables at the beginning of the script.

``` bash
# directory with the genome index files + name of the gene annotation file

genome=/data/Bspc-training/shared/rnaseq_jan2025/human_GRCh38/
gtf=/data/Bspc-training/shared/rnaseq_jan2025/human_GRCh38/gencode.v47.primary_assembly.annotation.gtf
```

We'll create output directories, but with the `-p` option. This will make sure that `mkdir` will create the directory only if it does not exist, it won't throw an error if it does exist, and even if we make a nested directory

``` bash
# make all of the output directories
# The -p option means mkdir will create the whole path if it does not exist and refrain from complaining if it does exist
# Let's add the suffix "auto" to the names to differentiate it from the other output folders we already have

mkdir -p results/fastqc_auto/
mkdir -p results/STAR_auto/
mkdir -p results/qualimap_auto/
```

Now that we have already created our output directories, we can now specify variables with the path to those directories both for convenience but also to make it easier to see what is going on in a long command.

``` bash
# set up output filenames and locations. 

fastqc_out=results/fastqc_auto/
align_out=results/STAR_auto/${samplename}_
align_out_bam=results/STAR_auto/${samplename}_Aligned.sortedByCoord.out.bam
qualimap_out=results/qualimap_auto/${samplename}.qualimap
```

## Keeping track of tool versions

All of our variables are now staged. Next, let's make sure all the modules are loaded. This is also a good way to keep track of the versions of tools that you are using in the script:

``` bash
# set up the software environment (use version numbers)

module load fastqc/0.12.1
module load STAR/2.7.11b
module load qualimap/2.2.1
unset DISPLAY #so qualimap doesn't try to render results on the command line
```

### Preparing for future debugging

In the script, it is a good idea to use `echo` for debugging. `echo` basically displays the string of characters specified within the quotations. When you have strategically place `echo` commands specifying what stage of the analysis is next, in case of failure you can determine the last `echo` statement displayed to troubleshoot the script. For example:

``` bash
echo "Processing file ${fq}"
```

> You can also use `set -x`:
>
> `set -x` is a debugging tool that will make bash display the command before executing it. In case of an issue with the commands in the shell script, this type of debugging lets you quickly pinpoint the step that is throwing an error. Often, tools will display the error that caused the program to stop running, so keep this in mind for times when you are running into issues where this is not available. You can turn this functionality off by saying `set +x`

## Running the tools using variables

Now that we have established all of our variables, now we can use them in the context of running the three tools: FASTQC, STAR and Qualimap.

**Running FASTQC:**

Remember that FASTQC needs the following arguments, which we can now fill with variables that we set up above: Output directory, number of threads to use, list of input file(s).

``` bash
echo "Starting FASTQC: ${samplename}"

# Run FastQC and move output to the appropriate folder
fastqc -o ${fastqc_out} -t ${cores} ${fq} 

echo "Finished FASTQC: ${samplename}"
```

**STAR**:

``` bash
echo "Mapping: ${samplename}"

# Run STAR
STAR \
  --runThreadN ${cores} \
  --genomeDir ${genome} \
  --readFilesIn ${fq} \
  --outFileNamePrefix ${align_out} \
  --outSAMtype BAM SortedByCoordinate \
  --outSAMunmapped Within \
  --outSAMattributes Standard \
  --outTmpDir=/lscratch/$SLURM_JOB_ID/STARtmp
```

``` bash
echo "Running Qualimap: ${samplename}"

# Run Qualimap
qualimap rnaseq \
  -outdir ${qualimap_out} \
  -a proportional \
  -bam ${align_out_bam} \
  -p strand-specific-reverse \
  -gtf ${gtf} \
  --java-mem-size=8G
```

### Last addition to the script

It is best practice to have the script **usage** specified at the top any script. This should have information such that when your future self, or a co-worker, uses the script they know what it will do and what input(s) are needed. For our script, we should have the following lines of comments right at the top after `#!/bin/bash/`:

``` bash
# This script takes a fastq file of RNA-seq data, runs FastQC, STAR and Qualimap
# USAGE: sh rnaseq_analysis_on_input_file.sh <name of fastq file>
```

It is okay to specify this after everything else is set up, since you will have most clarity about the script only once it is fully done.

Additional helpful notes might be what the user should expect in terms of output file and directories.

### Saving and running script

To transfer the contents of the script from your laptop to Biowulf, you can copy and paste the contents into a new file called `rnaseq_analysis_on_input_file.sh` using `vim`.

``` bash
$ cd /data/rnaseq/scripts/

$ vim rnaseq_analysis_on_input_file.sh 
```

> *Alternatively, you can save the script on your computer and transfer it to your `/rnaseq/scripts/` directory using the mounted directory system or `scp`*

We should all have an interactive session with 12 cores, so we can run the script as follows from your `/rnaseq/` directory:

``` bash
$ sh scripts/rnaseq_analysis_on_input_file.sh /data/Bspc-training/$USER/rnaseq/raw_data/Mov10_oe_1.subset.fq
```

## Running the script to submit jobs in parallel to the Slurm scheduler

The above script will run in an interactive session **one file at a time**. But the whole point of writing this script was to run it on all files at once. How do you think we can do this?

To run the above script **"in serial"** for all of the files on a worker node via the job scheduler, we can create a separate submission script that will need 2 components:

1.  **Slurm directives** at the **beginning** of the script. This is so that the scheduler knows what resources we need in order to run our job on the compute node(s) so we don't need to specify them from the command line.
2.  a **`for`** loop that iterates through and runs the above script for all the fastq files.

Below is what this second script (`rnaseq_analysis_on_allfiles.slurm`) would look like BUT DON'T RUN THIS :

``` bash
#!/bin/bash
#SBATCH --gres=lscratch: #designate lscratch space based on size of input FASTQ file
#SBATCH --partition=quick
#SBATCH --time=3:00:00     # time limit
#SBATCH --cpus-per-task=       # number of cores
#SBATCH --mem=   # requested memory
#SBATCH --job-name salmon_in_serial      # Job name
#SBATCH -o .out           # File to which standard output will be written
#SBATCH -e .err       # File to which standard error will be written
#SBATCH --mail-type=BEGIN,END

# this `for` loop, will take the fastq files as input and run the script for all of them one after the other. 
for fq in ~/rnaseq/raw_data/*.fq
do
  echo "running analysis on ${fq}"
  sh ~/rnaseq/scripts/rnaseq_analysis_on_input_file.sh ${fq}
done
```

**But we don't want to run the analysis on these 6 samples one after the other!** We want to run them "in parallel" as 6 separate jobs.

**Note:** If you create and run the above script, or something similar to it, i.e. with Slurm directives at the top, you should give the script name `.run` or `.slurm` as the extension. This will make it obvious that it is meant to submit jobs to the Slurm scheduler.

------------------------------------------------------------------------

**Exercise**

How would you run `rnaseq_analysis_on_allfiles.slurm`, i.e. the above script?

------------------------------------------------------------------------

## Parallelizing the analysis for efficiency - convert to SWARM

Parallelization will save you a lot of time with real (large) datasets. To parallelize our analysis, we will still need to write a second script that will call the script we just wrote that takes a fastq file as input (rnaseq_analysis_on_input_file.sh). We will still use a `for` loop, but we will be creating a regular shell script and we will be specifying the Slurm directives differently.

> Alternatively, this could also be done using a ***Slurm array***, which lets you submit a collection of similar jobs easily and quickly. You can learn more about Slurm arrays [here](https://hbctraining.github.io/Training-modules/Intermediate_shell/lessons/arrays_in_slurm.html).

Use `vim` to start a new shell script called `rnaseq_analysis_on_allfiles-for_slurm.sh`:

``` bash
$ vim rnaseq_analysis_on_allfiles_for-slurm.sh
```

This script loops through the same files as in the previous (demo) script, but the command being submitted within the `for` loop is `sbatch` with Slurm directives specified on the same line:

``` bash
#! /bin/bash

for fq in ~/rnaseq/raw_data/*.fq
do

sbatch -c 6 --job-name rnaseq-workflow --mem 8G --wrap="sh ~/rnaseq/scripts/rnaseq_analysis_on_input_file.sh ${fq}"
sleep 1 # wait 1 second between each job submission
  
done
```

> Please note that after the `sbatch` directives the command `sh ~/rnaseq/scripts/rnaseq_analysis_on_input_file.sh ${fq}` is in quotes.

``` bash
$ cd ~/rnaseq/scripts/
$ sh rnaseq_analysis_on_allfiles_for-slurm.sh
```

What you should see on the output of your screen would be the jobIDs that are returned from the scheduler for each of the jobs that your script submitted.

You can use `O2sacct` to check progress. And we can check if there are any additional files in our analysis folder.

``` bash
$ O2sacct

$ ls -l /n/scratch/users/r/$USER/rnaseq_hbc-workshop/
```

Don't forget about the `scancel` command, should something go wrong and you need to cancel your jobs.

> **NOTE:** All job schedulers are similar, but not the same. Once you understand how one works, you can transition to another one without too much trouble. They all have their pros and cons which are considered by the system administrators when picking one for a given HPC environment. Some examples of other job schedulers are LSF, SGE, PBS/Torque.

------------------------------------------------------------------------

*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
