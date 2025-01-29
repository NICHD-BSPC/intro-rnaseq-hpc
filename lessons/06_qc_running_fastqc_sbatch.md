---
title: "Quality control using FASTQC - script running"
author: "Harvard HPC Staff, Modified by Sally Chang @ NICHD"
date: Last edited January 2025
duration: 45 minutes
---

## Learning Objectives:

-   Modify a template SLURM job submission script
-   Run a SLURM job submission script to perform quality assessment for a full-sized FASTQ file

## Preparing SLURM Directives

### Performing quality assessment using job submission scripts

So far in our FASTQC analysis, we have been directly submitting commands to Biowulf using an interactive session (ie. `fastqc -o ../results/fastqc/ -t 6 *.fq`). However, we can submit a command or series of commands to these partitions using job submission scripts. This is extremely useful when we want to run a job that will take longer than we want to wait around for in an interactive job or will take more memory that we can use during an interactive job. *In this case, we will be running the analysis on the full-sized FASTQ file for one of our samples, which is 10x larger than the subset files we have thus far been working with.*

**Job submission scripts** for Biowulf are regular shell (command line) scripts, but contain the Slurm **options/directives** for our job submission. These directives define the various resources we are requesting for our job (i.e *number of cores, name of partition, runtime limit* ).

Submission of the script using the `sbatch` command allows Slurm to run your job when its your turn. Let's create a job submission script to automate what we have done in the previous lesson.

Our script will do the following:

1.  Specify important details to the SLURM scheduler
2.  Load the FastQC module
3.  Run FastQC on the full-sized `Mov10_oe_1` FASTQ file

Let's first change the directory to `/rnaseq/scripts`, and copy a template that already has some SLURM parameters listed and re-name this file something informative:

``` bash
$ cd rnaseq/scripts
$ cp /data/Bspc-training/shared/rnaseq_jan2025/slurm_template.sh .
$ mv slurm_template.sh mov10_fastqc_full.sh
```

Open the renamed file with Vim and take a look at the components.The first thing needed in a SLURM script is the **shebang line**. This tells the command line to interpret this file as a BASH script:

``` bash
#!/bin/bash
```

Following the shebang line are the Slurm directives. These are directives that your instructor finds particularly informative. If you are using these instructions without access to the course directory, you could also copy and paste this into a text editor :

``` bash
#SBATCH --job-name=
#SBATCH --partition=
#SBATCH --mail-type=                # Mail events 
#SBATCH --ntasks=                  # Run a single task     
#SBATCH --cpus-per-task=            # Number of CPU cores per task
#SBATCH --mem=                    # Job memory request
#SBATCH --time=            # Time limit hrs:min:sec
#SBATCH --output=          # Standard output log
#SBATCH --error=           #error log
```

-   `--job-name=` name of the job (ex. mov10_fastqc_full)
-   `--partition=` : name of compute partition (ex. norm, largemem, quick) [Biowulf partitions](https://hpc.nih.gov/docs/userguide.html#partitions)
-   `--time=` : how much time to allocate to job (ex. 08:00:00 for 8 hours) [Walltime limits for partitions](https://hpc.nih.gov/docs/userguide.html#wall)
-   `--cpus-per-task=`: Number of CPUs required (ex. '8' for 8 CPUS that can be used for multithreading.
-   `--mem=`: maximum memory, i.e. 8g (8 gigabytes - note the `g`)
-   `--mail-type=`: send an email to your NIH account when job starts, ends or quits with an error (ex. END, ALL)
-   `--output=`: The name of a log file to send output that would normally be printed to screen (ex. mov10_fastqc.log)
-   `--error=`: The name of a log file specifically to capture error messages

More about these options and many others can found in the [Biowulf User Guide](https://hpc.nih.gov/docs/userguide.html).

Use Vim in insert mode (`i`) to modify your copied file to use 1 task that makes use of 6 CPUs, which will run on the `quick` partition for one hour and use 6 gigabytes of memory. Make sure to also give your job and output/error logs informative names!

``` bash
#!/bin/bash
#SBATCH --job-name=mov10_oe1_full_fastqc
#SBATCH --partition=quick
#SBATCH --mail-type=ALL                # Mail events 
#SBATCH --ntasks=1                  # Run a single task     
#SBATCH --cpus-per-task=6           # Number of CPU cores per task
#SBATCH --mem=6g                    # Job memory request
#SBATCH --time=01:00:00            # Time limit hrs:min:sec
#SBATCH --output=%j.out          # Use job name as a variable
#SBATCH --error=%j.err           # Use job name as a variable
```

## Modify the body of the script

Now in the body of the script, we can include any commands we want to run, specifying that our input file is now a file in the shared class directory. In this case, it will be the following:

``` bash
## Load modules required for script commands
module load fastqc/0.12.1

## Run FASTQC on a full-sized file in the class shared directory
fastqc -o /data/Bspc-training/$USER/rnaseq/results/fastqc /data/Bspc-training/shared/rnaseq_jan2025/Mov10_oe_1.fq
```

> **NOTE:** These are the same commands we used when running FASTQC in the interactive session. Since we are writing them in a script, the `tab` completion function will **not work**, so please make sure you don't have any typos when writing the script!

Once done with your script, click `esc` to exit the INSERT mode. Then save and quit the script by typing `:wq`. You may double check your script by typing `less mov10_fastqc.run`.

## Submit the job script

Now, if everything looks good submit the job!

``` bash
$ sbatch mov10_fastqc.run
```

You should immediately see a prompt saying `Submitted batch job JobID`. Your job is assigned with that unique identifier `JobID`. You can check on the status of your job with:

``` bash
$ squeue -u $USER
```

Look for the row that corresponds to your `JobID`. The third column indicates the state of your job. Possible states include `PENDING (PD)`, `RUNNING (R)`, `COMPLETING (CG)`. For this example, once your job is `RUNNING`, you should expect it to finish in less than two minutes (although of course this will not always be the case for longer analyses).

> **NOTE:** Here is a [table of Job State Codes](https://hpc.nih.gov/docs/userguide.html#states) from the Biowulf User Guide.

Check out the output files in your directory:

``` bash
$ ls -lh ../results/fastqc/
```

There should also be one standard error (`.err`) and one standard out (`.out`) files from the job listed in `~/rnaseq/scripts`. You can move these over to your `logs` directory and give them more intuitive names:

``` bash
$ mv *.err ../logs/fastqc.err
$ mv *.out ../logs/fastqc.out
```

> **NOTE:** The `.err` and `.out` files store log information during the script running. They are helpful resources, especially when your script does not run as expected and you need to troubleshoot the script.

------------------------------------------------------------------------

**Exercise**

1\. Take a look at what's inside the `.err` and `.out` files. What do you observe? Do you remember where you see those information when using the interactive session?

2\. How would you change your script to analyze the 6 files we used in the last episode?

------------------------------------------------------------------------

*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*

-   *The materials used in this lesson was derived from work that is Copyright © Data Carpentry (<http://datacarpentry.org/>). All Data Carpentry instructional material is made available under the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0).*
