---
title: Working on HPC
author: Harvard HPC Staff, Modified by E. Sally Chang @ NICHD
date: Last Edited January 2025
---

# Working in an HPC environment

## Learning Objectives

-   Define the basic components of a computing environment
-   Explain parallelization in computing
-   Describe the salient features of the SLURM job scheduling system

## Clarifying some terminology

### Cores

CPUs (Central Processing Units) are composed of a single or multiple cores. If I have a processor with 4 "cores", i.e. a "quad-core" processor. This means that my processor can do 4 distinct computations at the same time.

### Data Storage

Data storage is measured in bytes, usually Gigabytes (GB), Terabytes (TB), Petabytes (PB).

### Memory

Memory, in an HPC setting, refers to volatile or temporary information used by running processes - typically this refers to RAM (random access memory). This is distinct from data storage we discussed in the previous point. It is an essential component of any computer, as this is where data is stored *when some computation is being performed* on it. If you have ever used R, all the objects in your environment are usually stored in memory until you save your environment. My laptop has 16 GB of memory available to it, for example.

## Why use the cluster or an HPC environment?

<p align="center">

<img src="https://github.com/hbctraining/Intro-to-shell-flipped/blob/master/img/compute_cluster.png?raw=true" width="450"/>

</p>

1.  A lot of software is designed to work with the resources on an HPC environment and is either unavailable for, or unusable on, a personal computer.
2.  If you are performing analysis on large data files (e.g. high-throughput sequencing data), you should work on the cluster to avoid issues with memory and to get the analysis done a lot faster with the superior processing capacity. Essentially, Biowulf has:
    -   96,000 processor cores
    -   40+ Petabytes of storage!
    -   Anywhere from 128 GB - 3TB of memory depending on the node!

Check out

### Parallelization

Point 2 in the last section brings us to the idea of **parallelization** or parallel computing that enables us to efficiently use the resources available on the cluster.

What if we had 3 input files that we wanted to analyze? Well, we could process these files **in serial**, i.e. use the same core(s) over and over again, with or without multithreading, as shown in the image below.

<p align="center">

<img src="../img/serial_hpc_3samples.png" width="450"/>

</p>

This is great, but it is not as efficient as multithreading each analysis, and using a set of 8 cores for each of the three input samples. With this type of parallelization, several samples can be analysed at the same time!

<p align="center">

<img src="../img/multithreaded_hpc_3samples.png" width="650"/>

</p>

## Refresher: Connect to a *login* node on Biouwulf

Let's get started with the hands-on component by typing in the following command to log in to O2:

``` bash
ssh username@biowulf.nih.gov
```

You will receive a prompt for your password, and you should type in your associated password; **note that the cursor will *not move* as you type in your password**.

A warning might pop up the first time you try to connect to a remote machine, type "Yes" or "Y".

## Specify options for a *compute* node on Biowulf

You can access compute nodes in 2 ways.

1\. Directly using an interactive session (`sinteractive` command): This is essentially a way for us to do work on the compute node directly from the terminal. **If the connectivity to the cluster is lost in the middle of a command being run that work will be lost in an interactive session.**

2\. By starting a "batch" job (Slurm's `sbatch` command): The `sbatch` command with a few mandatory parameters + a specialized shell script will result in the script being run on a compute node. This "job" will not be accessible directly from the Terminal and will run in the background. **Users do not need to remain connected to the cluster when such a "batch job" is running.**

For now let's start an interactive session, but specify more more particulars:

``` bash
# This is an example, you don't need to run this
$ sinteractive --time=02:00:00 --mem==1g 
```

In the above command the parameters we are using are requesting specific resources:

`--time=02:00:00` - time needed for this work: 0 days, 2 hours, 0 minutes.

`--mem=1G` - memory needed - 1 gigibyte (GiB)

> These parameters are used for `sbatch` as well, but they are listed differently within the script used to submit a batch job. We will be reviewing this later in this lesson.

Let's check how many jobs we have running currently, and what resources they are using.

``` bash
$ squeue -u $USER
```

## More about Slurm

-   Slurm = Simple Linux Utility for Resource Management
-   Fairly allocates access to resources (computer nodes) to users for some duration of time so they can perform work
-   Provides a framework for starting, executing, and monitoring batch jobs
-   Manages a queue of pending jobs; ensures that no single user or core monopolizes the cluster

### Requesting resources from Slurm

Below are some of the arguments you can specify when requesting resources from Slurm for both `srun` and `sbatch`:

* `--partition=` : name of compute partition (ex. norm (default) largemem, quick) [Biowulf partitions](https://hpc.nih.gov/docs/userguide.html#partitions)        
* `--time=` : how much time to allocate to job (ex. 08:00:00 for 8 hours) [Walltime limits for partitions](https://hpc.nih.gov/docs/userguide.html#wall) 
* `--cpus-per-task=`: Number of CPUs required (ex. '8' for 8 CPUS)
*  `--mem=`: maximum memory, i.e. 8G (8 gigabytes)
* `--job-name=` name of the job (ex. Fastqc_run, rnaseq_workflow_mov10) 
* `--mail-type=`: send an email to your NIH account when job starts, ends or quits with an error (ex. END, ALL)
* `--output=`: The name of a log file to send any output that would normally be printed to screen (ex. mov10_fastqc.log)

More about these options and many others can found in the [Biowulf User Guide](https://hpc.nih.gov/docs/userguide.html).                                           

### `sbatch` job submission script: Example

An `sbatch` job submission script is essentially a normal shell script with the Slurm resource request specified at the top (Slurm directives) preceded by `#SBATCH`. Below is an example of an sbatch shell script that is requesting the following:

-   The "quick" partition: For jobs \<4 hours, which are scheduled at higher priority.
-   2 hours
-   4 CPUS
-   Using 400 megabytes (100MiB for each core)

***DO NOT RUN - you'll practice writing one of these yourself soon.***

```         
#! /bin/sh

#SBATCH --partition=quick
#SBATCH –-time=02:00:00
#SBATCH –-cpus-per-task=4
#SBATCH --mem=400M
#SBATCH –o %j.out #where to send output messages
#SBATCH –e %j.err #where to send error messages
#SBATCH -J fastqc_run 
#SBATCH --mail-type=ALL
#SBATCH –-mail-user=user@nih.gov

## Load the fastqc module
module load fastqc/0.11.5

# run fastqc (multithreaded)
fastqc -t 4 file1_1.fq file1_2.fq file2_1.fq file2_2.fq
```

## Using software on Biowulf

### Brief Intro to \$PATH and other Environmental Variables

Environment variables are, in short, variables that describe the environment in which programs run, and they are predefined for a given computer or cluster that you are on.

Here are some of the most

-   `$USER` Recall how that we can use `$USER` as a variable instead of actually writing out our username every time? That's because this is a built-in environmental variable on

-   \$HOME defines the full path for the home directory of a given user. Try typing `echo $HOME` to confirm!

-   `$PATH` defines a list of directories to search in when looking for a command/program to execute.

In this lesson, we are going to focus on \$PATH. If you want to see what is on your \$PATH, you can use the echo command again. As a beginning Biowulf user, you probably don't have as many entries as I do!

``` bash
$ echo $PATH

/gpfs/gsfs11/users/changes/mambaforge/condabin:/usr/local/slurm/bin:/usr/local/bin:/usr/X11R6/bin:/usr/local/jdk/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/usr/local/mysql/bin:/home/changes/opt/bin:/data/changes/mambaforge/condabin:/data/NICHD-core1/bin
```

This output is a lot more complex than the output from \$HOME! When you look closely at the output of `echo $PATH`, you should a list of full paths separated from each other by a ":".

### What are all these paths? And what do they represent?

These are the directories that the shell will look through (in the same order as they are listed) for any given command or executable file that you type on the command prompt.

For example, we have been using the `ls` command to list contents in a directory. When we type `ls` at the command prompt, the shell searches through each path in `$PATH` until it finds an executable file called `ls`. So which of those paths contain that executable file?

For any command you execute on the command prompt, you can find out where the executable file is located using the `which` command.

What path was returned to you? Does it match with any of the paths stored in `$PATH`?

Try it on a few of the basic commands we have learned so far:

```         
$ which <your favorite command>
$ which <your favorite command>
```

Check the path `/usr/bin/` and see what other executable files you recognize. (Note that executable files will be listed as green text or have the `*` after their name).

```         
$ ls -lF /usr/bin/
```

The path `/usr/bin` is often where executables for commonly used commands are stored, but things get a little more complicated when you start working with things like Conda environments.

> As pointed out earlier, a lot of the folders listed in the `$PATH` variable are called `bin`. This is because of a convention in Unix to call directories that contain all the commands (in ***binary*** format) **`bin`**.

### LMOD system

In the next example, we want to run the FastQC tool to look at sequence quality. However, before we use the `fastqc` command, we have to make sure it is in our `$PATH` - luckily there is an easy way of doing this on Biowulf!

``` bash
$ module load fastqc 
```

This `module load` command is part of the LMOD system available on Biowulf. It enables users to access software installed on Biowulf easily, and manages every software's dependency. The LMOD system adds directory paths of software executables and their dependencies (if any) into the `$PATH` variable.

Some key LMOD commands are listed below:
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

------------------------------------------------------------------------

#### Exercises

1.  What are the contents of the `$PATH` environment variable?
2.  Try running the `multiqc` command. What do you get as output?
3.  Check if the `multiqc` tool is available as a module. How many versions of `multiqc` are available?
4.  Load `multiqc/1.9`. Did you have to load any additional modules?
5.  List all the modules loaded in your environment
6.  Try running the `multiqc` command again. What do you get as output?
7.  Use the `which` command to find out where the `multiqc` tool is in the file system.
8.  Check the contents of the `$PATH` environment variable again. Any changes compared to before?

------------------------------------------------------------------------

## Filesystems on Biowulf

-   Storage on HPC systems is organized differently than on your personal machine.
-   Each node on the cluster does not have storage; instead, it is on disks bundled together externally.
-   Storage filesystems can be quite complex, with large spaces dedicated to a pre-defined purpose.
-   Filesystems are accessed over the internal network by all the nodes on the cluster.
-   It is often best practice to use temporary storage spaces (called "scratch" spaces, see below) while running analyses. We'll use scratch space in a few exercises during this course.
-   **Most importantly: Disk space on the NIH HPC should never be used as archival storage. *Before you even get your project data in the first place, think about where you'll store it after the project is done.***
-   So much more information about storage best practices can be found on the [Biowulf Storage page](https://hpc.nih.gov/storage/index.html).

### Data Management and Checking your Quota

Use the `checkquota` command to determine how much disk space you are using.

``` bash
$ checkquota

Mount                   Used      Quota  Percent    Files    Limit  Percent
/data:                5.2 GB   500.0 GB    1.04%       23 31457280    0.00%
/data(SharedDir):     4.1 TB    10.5 TB   38.67%   103357 31876696    0.32%
/home:               10.3 GB    16.0 GB   64.62%       73      n/a    0.00%
```

You can also see your usage history in the [User Dashboard](https://hpcnihapps.cit.nih.gov/auth/dashboard), under the Disk Usage tab. In fact, I would recommend getting used to looking at what is available in your user dashboard!

**Data Management**

Biowulf disk storage is intended for active data and cannot be used for longterm archiving. To help users and groups manage their data on Biowulf, the [HPC User Dashboard](https://hpc.nih.gov/dashboard) provides [Krona](https://github.com/marbl/Krona/wiki/KronaTools)-based visualizations of their data directories that show the distribution and age of data. These interactive hierarchical maps of data and shared directories will assist users in finding unused or infrequently accessed old data that could be archived.
