---
title: "Quality control using FASTQC - script running"
author: "Harvard HPC Staff, Modified by Sally Chang @ NICHD"
date: Last edited January 2025
duration: 45 minutes
---

## Learning Objectives:

-   Modify a template SLURM job submission script
-   Run a SLURM job submission script to perform quality assessment for a full-sized FASTQ file

## Quality Control of FASTQ files

### Performing quality assessment using job submission scripts

So far in our FASTQC analysis, we have been directly submitting commands to Biowulf using an interactive session (ie. `fastqc -o ../results/fastqc/ -t 6 *.fq`). However, we can submit a command or series of commands to these partitions using job submission scripts. This is extremely useful when we want to run a job that will take longer than we want to wait around for in an interactive job or will take more memory that we can use during an interactive job. *In this case, we will be running the analysis on the full-sized FASTQ file for one of our samples, which is 10x larger than the subset files we have thus far been working with.*

**Job submission scripts** for Biowulf are regular shell (command line) scripts, but contain the Slurm **options/directives** for our job submission. These directives define the various resources we are requesting for our job (i.e *number of cores, name of partition, runtime limit* ).

Submission of the script using the `sbatch` command allows Slurm to run your job when its your turn. Let's create a job submission script to automate what we have done in the previous lesson.

Our script will do the following:

1.  Specify important details about our
2.  Load the FastQC module
3.  Run FastQC on the full-sized Mov10_1_oe FASTQ file

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

Following the shebang line are the Slurm directives. For the script to run, we need to include options for **queue/partition (--partition) and runtime limit (--time)**. To specify our options, we precede the option with `#SBATCH`. Some key resources to specify are:

| Resource  |      Flag      |                                                        Description                                                         |
|:----------------------:|:----------------------:|:----------------------:|
| partition |  --partition   |                                                       partition name                                                       |
|   time    |     --time     |                                hours:minutes run limit, after which the job will be killed                                 |
|   core    | -cpus-per-task | number of cores requested -- this needs to be greater than or equal to the number of cores you plan to use to run your job |
|  memory   |     --mem      |                                         memory limit per compute node for the job                                          |

Let's specify those options as follows:

``` bash
#SBATCH --partition=quick        # Based on running these analyses before, we can get away with running this on a quick (priority!) node
#SBATCH --time=02:00:00       # time limit
#SBATCH --cpus-per-task=6        # number of cores
#SBATCH --mem=6g   # requested memory
#SBATCH --job-name rnaseq_mov10_fastqc      # Job name
#SBATCH -o %j.out           # File to which standard output will be written
#SBATCH -e %j.err       # File to which standard error will be written
```

Now in the body of the script, we can include any commands we want to run. In this case, it will be the following:

``` bash
## Change directories to where the fastq files are located. For now I will use absolute paths until we find a good solution for a prefix. 
cd /data/NICHD-core0/test/changes/rc_training/rnaseq/raw_data

## Load modules required for script commands
module load fastqc/0.12.1

## Run FASTQC
fastqc -o /data/NICHD-core0/test/changes/rc_training/rnaseq/fastqc/ -t 6 *.fq
```

> **NOTE:** These are the same commands we used when running FASTQC in the interactive session. Since we are writing them in a script, the `tab` completion function will **not work**, so please make sure you don't have any typos when writing the script!

Once done with your script, click `esc` to exit the INSERT mode. Then save and quit the script by typing `:wq`. You may double check your script by typing `less mov10_fastqc.run`. If everything looks good submit the job!

``` bash
$ sbatch mov10_fastqc.run
```

You should immediately see a prompt saying `Submitted batch job JobID`. Your job is assigned with that unique identifier `JobID`. You can check on the status of your job with:

``` bash
$ squeue -u changes
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

**Exercise** 1. Take a look at what's inside the `.err` and `.out` files. What do you observe? Do you remember where you see those information when using the interactive session? 2. How would you change the `mov10_fastqc.run` script if you had 9 fastq files you wanted to run in parallel?

------------------------------------------------------------------------

*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*

-   *The materials used in this lesson was derived from work that is Copyright © Data Carpentry (<http://datacarpentry.org/>). All Data Carpentry instructional material is made available under the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0).*
