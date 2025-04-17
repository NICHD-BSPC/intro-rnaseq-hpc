# Introduction to RNA-seq using high-performance computing (January 2025)

### Description

This repository has teaching materials for a 9-week Introduction to RNA-sequencing data analysis workshop held at NICHD from Jan-March 2025. This workshop focuses on teaching basic computational skills to enable the effective use of an high-performance computing environment to implement an RNA-seq data analysis workflow. The first four and a half weeks includes an introduction to shell (bash) and shell scripting, running the RNA-seq workflow from FASTQ files to count data and covers best practice guidelines for RNA-seq experimental design and data organization/management. The latter half covers differential expression analyses using DEseq2 and downstream analyses of those results. 

### Learning Objectives

* Gain practical knowledge about analyzing RNAseq from experimental design through functional enrichment analysis
* Learn broadly applicable bioinformatics skills such as command line and R programming
* Work with real data sets and real bioinformatics environments on NIH’s high-performance compute cluster (Biowulf)
* Apply what you learn to your own bulk RNAseq data

### Weekly Course Components

*	**In-person lecture and hands on demonstrations**: Tuesdays, 10am-12pm starting January 28th
*	**In-person practice sessions**: Fridays at 10am-12pm starting on January 31st
    - This is a chance to review material from Tuesday and work on the weekly assignments (see below) with help from Dr. Chang, BSPC staff as well as fellow students. You are strongly encouraged to attend these sessions to get the most out of this course.
* **Asynchronous Practice**: 3-4 hours of asynchronous practice and reading. This will prepare you for the following week and allow you to apply knowledge to new examples. You will submit small weekly assignment to me so I can check your progress and adjust content if needed. 


### Topics and Links to weekly materials

Note: Links to weekly lessons will be provided several days in advance in the table below. All course materials are living documents - some changes may be made after publishing to reflect our progress through the material. 


| Link to Materials                      | Topic                                                     |
|-----------------|-------------------------------------------------------|
| [Week 1](schedule/links-to-lessons.md#week-1) | Introduction to the command line and logging into Biowulf |
| [Week 2](schedule/links-to-lessons.md#week-2) | Scientific software on Biowulf, quality control of sequence data, experimental design|
| [Week 3](schedule/links-to-lessons.md#week-3) | Reference genomes, theory and practice of mapping RNAseq reads to a reference, Mapping QC|
| [Week 4](schedule/links-to-lessons.md#week-4) | Theory and practice of counting RNAseq reads|
| [Week 5](schedule/links-to-lessons.md#week-5) | Automation of the RNAseq workflow, transition to RStudio for differential expression analyses | 
| [Week 6](schedule/links-to-lessons.md#week-6) | Overview and prep for DEseq2 analysis pipeline, assessing sample quality, design formulas |
| [Week 7](schedule/links-to-lessons.md#week-7) | Hypothesis testing in DESeq2, summarizing and visualizing results |
| [Week 8](schedule/links-to-lessons.md#week-8) | Likelihood ratio test results, overenrichment analyses, functional class scoring | 
| [Week 9](schedule/links-to-lessons.md#week-9) | Any material left from Week 8, course wrap-up, review of requested topics | 

### Course communication: 
* For any content questions: We will be using the `rnaseq_jan2025` Slack channel under the Bioinformatics @NICHD workspace. Posting your questions in this channel is great, as multiple staff and fellow students will be able to respond. Also, you probably aren't the only one with that question!
* For scheduling/personal concerns: Contact Dr. Chang by e-mail. 


### Software Requirements

***Mac users:***

-   Plain text editor: TextEdit should be installed by default on Macs. 

***Windows users:***

-   [GitBash](https://git-scm.com/download/win)
-   Plain text editor: Microsoft Notepad should be installed by default on Windows.

**Note about text plain text editors:** A plain text editor is a program to edit text files such as a script that doesn’t interfere with formatting like a full word processor (like Word) would. The built-in text editors for each operating system are listed above. 
 
As we progress in the course and your research, you may find that you want a "fancier" text editor that has more features for coding efficiently. 

The following options are approved by NICHD IT but may need a license after a trial period: 
* For Macs: [BBEdit](https://www.barebones.com/products/bbedit/index.html), which has a lot more features for coding but needs a license to keep using all of those features after the trial period.
* For PCs: [NotePad++](https://notepad-plus-plus.org/), which likewise has many more features and is FREE!

  

------------------------------------------------------------------------

### Citation

These materials were modified by E. Sally Chang at NICHD from the following citation:

> Mary E. Piper, Meeta Mistry, Jihe Liu, William J. Gammerdinger, & Radhika S. Khetani. (2022, January 10). hbctraining/Intro-to-rnaseq-hpc-salmon-flipped: Introduction to RNA-seq using Salmon Lessons from HCBC (first release). Zenodo. <https://doi.org/10.5281/zenodo.5833880>. RRID:SCR_025373.

------------------------------------------------------------------------

*The original materials developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*

-   *Some materials used in these lessons were derived from work that is Copyright © Data Carpentry (<http://datacarpentry.org/>). All Data Carpentry instructional material is made available under the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0).*
