# Introduction to RNA-seq using high-performance computing at NICHD: May 2025
### Description

This repository contains materials for an introductory bulk RNA sequencing analysis course, held at NICHD from April 30th-July 9th, 2025. 

This workshop focuses on teaching basic computational skills to enable the effective use of an high-performance computing environment to implement an RNA-seq data analysis workflow. The first half of the course includes an introduction to shell (bash) and shell scripting, running the RNA-seq workflow from FASTQ files to count data and covers best practice guidelines for RNA-seq experimental design and data organization/management. The latter half covers differential expression analyses using the DEseq2 R package and downstream analyses of those results. 

### Learning Objectives

* Gain practical knowledge about analyzing RNAseq from experimental design through functional enrichment analysis
* Learn broadly applicable bioinformatics skills such as command line and R programming
* Work with real data sets and real bioinformatics environments on NIH’s high-performance compute cluster (Biowulf)
* Apply what you learn to your own bulk RNAseq data


### Topics and Links to weekly materials

Please note that this schedule is tentative and subject to change. Materials for each week will be made available prior to the Wednesday session for each week. 

| Link to Materials                      | Topic                                                     |
|-----------------|-------------------------------------------------------|
| Week 1 | Introduction to the command line and logging into Biowulf |
| Week 2 | Scientific software on Biowulf, quality control of sequence data, experimental design|
| Week 3 | Reference genomes, theory and practice of mapping RNAseq reads to a reference |
| Week 4 | Theory and practice of counting RNAseq reads, automation of RNAseq workflow | 
| Week 5 | Catch up as needed, intro to R/RStudio| 
| Week 6 | Overview and prep for DEseq2 analysis pipeline, assessing sample quality|
| Week 7 | Design formulas and hypothesis testing in DESeq2 (Wald and Likelihood Ratio tests) |
| Week 8 | Summarizing and visualizing results, annotation databases | 
| Week 9 | Overenrichment analyses, functional class scoring, other functional analyses | 
| Week 10 (spread over July 2nd and 9th) | Make-up time, Course wrap-up, review of requested topics | 


### Course Scheduling
*	In-person sessions: Wednesdays and Fridays, 10am-12pm starting April 30th. These are composed of lectures, live demonstrations and in-person practice time
*	Asynchronous Practice: Approximately two hours of asynchronous practice and reading. This will prepare you for the following week and allow you to apply knowledge to new examples. You will submit small weekly assignment to me so I can check your progress and adjust content if needed. 
*	Optional office hours: 11am Monday on Teams (link will be distributed in class) starting on May 5th. These are reserved for troubleshooting that absolutely can't be done via Slack (see below). 

### Course Communication
*	For content questions: We will be using the `rnaseq_may2025` Slack channel under the Bioinformatics @NICHD workspace. You will be invited shortly before the course starts. **You are STRONGLY encouraged to ask any and all of your course-related questions in this Slack channel, where other students can see the answers and more people can chime in with help**. You are never the only person with a particular question!
*	For scheduling/personal concerns: Contact Dr. Chang by e-mail. 

### Weekly Assignments/Asynchronous Work 

There will be short assignments embedded in the weekly course materials, and will be **due at 12pm Tuesday of the following week**. This will give me enough time to give you feedback and adjust upcoming materials as needed before the course starts. For example, any assignments that are part of the Wed/Fri Week 1 lessons will be due Tuesday of Week 2. 

You will submit these assignments by copying requisite materials in a specific location in your course Biowulf directory. This will be demonstrated during Week 1. 

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
