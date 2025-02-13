---
title: "SAM/BAM file and assessing quality "
author: Harvard HPC Staff, Adapted by Sally Chang @ NICHD
date: Last modified February 2025
---

## Learning objectives

-   Understanding the standard alignment file (SAM/BAM) structure
-   Using `samtools` to evaluate alignment quality

## Assessing alignment quality

After running our single **subset** FASTQ file through the STAR aligner, you should have noticed a number of output files in the `rnaseq/results/STAR` directory. Let's take a quick look at some of the files that were generated and explore the content of some of them.

``` bash
$ sinteractive --mem=4g --cpus-per-task=2 #make sure you are on an interactive node with a bit of memory
```

```         
$ cd /data/$USER/rnaseq/results/STAR 

$ ls -lh
```

What you should see, is that for each FASTQ file you have **5 output files** and a single tmp directory. Briefly, these files are described below:

-   `Log.final.out` - a summary of mapping statistics for the sample
-   `Aligned.sortedByCoord.out.bam` - the aligned reads, sorted by coordinate, in BAM format. **We are going to be focusing on this file for this lesson!**
-   `Log.out` - a running log from STAR, with information about the run
-   `Log.progress.out` - job progress with the number of processed reads, % of mapped reads etc., updated every \~1 minute
-   `SJ.out.tab` - high confidence collapsed splice junctions in tab-delimited format. Only junctions supported by uniquely mapping reads are reported

## Alignment file format: SAM/BAM

The output we requested from the STAR aligner (using the appropriate parameters) is a BAM file. By default STAR will return a file in SAM format. BAM is a binary, compressed version of the SAM file, also known as **Sequence Alignment Map format**. The SAM file, introduced is a tab-delimited text file that contains information for each individual read and its alignment to the genome. While we will go into some features of the SAM format, the paper by [Heng Li et al](http://bioinformatics.oxfordjournals.org/content/25/16/2078.full) provides a lot more detail on the specification.

The file begins with a **header**, which is optional. The header is used to describe source of data, reference sequence, method of alignment, etc., this will change depending on the aligner being used. Each section begins with character ‘\@’ followed by a two-letter record type code. These are followed by two-letter tags and values. Example of some common sections are provided below:

```         
@HD  The header line
VN: format version
SO: Sorting order of alignments

@SQ  Reference sequence dictionary
SN: reference sequence name
LN: reference sequence length
SP: species

@PG  Program
PN: program name
VN: program version
```

Following the header is the **alignment section**. Each line that follows corresponds to alignment information for a single read. Each alignment line has **11 mandatory fields for essential mapping information** and a variable number of other fields for aligner specific information.

![](../img/sam_bam.png)

An example read mapping is displayed above. *Note that the example above spans two lines, but in the file it is a single line.* Let's go through the fields one at a time. First, you have the read name (`QNAME`), followed by a `FLAG`

> DISCUSSION: Is the first BAM read ID the same as the first input FASTQ read ID? Why or why not?

The `FLAG` value that is displayed can be translated into information about the mapping.

| Flag |                Description                |
|-----:|:-----------------------------------------:|
|    1 |              read is paired               |
|    2 |       read is mapped in proper pair       |
|    4 |             read is unmapped              |
|    8 |             mate is unmapped              |
|   16 |            read reverse strand            |
|   32 |            mate reverse strand            |
|   64 |               first in pair               |
|  128 |              second in pair               |
|  256 |           not primary alignment           |
|  512 | read fails platform/vendor quality checks |
| 1024 |     read is PCR or optical duplicate      |
| 2048 |          supplementary alignment          |

-   For a given alignment, each of these flags are either **on or off** indicating the condition is **true or false**.
-   The `FLAG` is a combination of all of the individual flags (from the table above) that are true for the alignment
-   The beauty of the flag values is that **any given combination of flags can only result in one sum**.

In our example we have a number that exist in the table, making it relatively easy to translate. But suppose our read alignment has a flag of 163 -- what does this translate to? It is the sum of 4 different flags:

`163 = 1 + 2 + 32 + 128`

Which tells us that:

1.  the read is mapped
2.  the read is mapped as part of a pair
3.  this is the mate reverse strand
4.  this read is the second of the pair

**There are tools that help you translate the bitwise flag, for example [this one from Picard](https://broadinstitute.github.io/picard/explain-flags.html)**

Technically, bitwise flags are binary numbers (here, 11 digits for the 12 items in the table) where each position represents yes/no for the respective item.

```         
00000000101
       ^^^^
       |||| read is paired (add 1 if true)
       ||| read is mapped as part of pair (add 2 if true)
       || read is unmapped (add 4 if true)
       | mate is unmapped (add 8 if true)
       ... and so on
```

When [reading binary](https://www.lifewire.com/how-to-read-binary-4692830), `00000000101` is `5` -- there's no other combination of true/false (1s and 0s) that could give us a value of 5, so this must be a paired read that is unmapped. Storing a digit in decimal character (like "5") is more space-efficient than storing the 12-character full binary representation on every line.

Moving along the fields of the SAM file, we then have `RNAME` which is the reference sequence name. The example read is from chromosome 1 which explains why we see 'chr1'. `POS` refers to the 1-based leftmost position of the alignment. `MAPQ` is giving us the alignment quality, the scale of which will depend on the aligner being used.

`CIGAR` is a sequence of letters and numbers that represent the *edits or operations* required to match the read to the reference. The letters are operations that are used to indicate which bases align to the reference (i.e. match, mismatch, deletion, insertion), and the numbers indicate the associated base lengths for each 'operation'.

| Operation |            Description            |
|----------:|:---------------------------------:|
|         M |    sequence match or mismatch     |
|         I |    insertion to the reference     |
|         D |      deletion from reference      |
|         N | skipped region from the reference |

Suppose our read has a CIGAR string of `50M3I80M2D` which translates to:

-   50 matches or mismatches
-   3 bp insertion
-   80 matches/mismatches
-   2 bp deletion

Now to the remaning fields in our SAM file:

![SAM1](../img/sam_bam3.png)

The next three fields are more pertinent to paired-end data. `MRNM` is the mate reference name. `MPOS` is the mate position (1-based, leftmost). `ISIZE` is the inferred insert size.

Finally, you have the data from the original FASTQ file stored for each read. That is the raw sequence (`SEQ`) and the associated quality values for each position in the read (`QUAL`).

## `samtools`

[SAMtools](http://samtools.sourceforge.net/) is a tool that provides alot of functionality in dealing with SAM files. SAMtools utilities include, but are not limited to, viewing, sorting, filtering, merging, and indexing alignments in the SAM format. In this lesson we will explore a few of these utilities on our alignment files. To use this we need to load the module.

``` bash
# loaded default version 1.21 
$ module load samtools
```

### Viewing the SAM file

**Note**: My output BAM file run on the subset of data is `Mov10_oe_1_subsetAligned.sortedByCoord.out.bam` . Yours might be slightly different if you used a different prefix while mapping, so make sure to double check what your output is called and adjust the commands below if neede

Now that we have learned so much about the SAM file format, let's use `samtools` to take a quick peek at our own files. The output we had requested from STAR was a BAM file. The problem is the BAM file is binary and not human-readable. Using the `view` command within `samtools` we can easily convert the BAM into something that we can understand -- and pipe to other command line tools! You will be returned to screen the entire SAM file, and so we can either write to file, or pipe this to the `less` command so we can scroll through it.

We are assuming you are in `/data/Bspc-training/$USER/rnaseq/results/STAR` for the following commands.

```         
$ samtools view -h Mov10_oe_1_subsetAligned.sortedByCoord.out.bam > Mov10_oe_1_subsetAligned.sortedByCoord.out.sam
```

``` bash
$ head Mov10_oe_1_subsetAligned.sortedByCoord.out.sam

@HD VN:1.4  SO:coordinate
@SQ SN:chr1 LN:248956422
@SQ SN:chr2 LN:242193529
@SQ SN:chr3 LN:198295559
@SQ SN:chr4 LN:190214555
@SQ SN:chr5 LN:181538259
@SQ SN:chr6 LN:170805979
@SQ SN:chr7 LN:159345973
@SQ SN:chr8 LN:145138636
@SQ SN:chr9 LN:138394717
```

> NOTE: what did `-h` do? How would you find out?

So, these first 10 lines aren't super interesting, as it is just listing the names of the input sequences (in this case, chromosomes of the input genome) and their lengths. Let's try `tail` for what some human-readable lines of actual SAM-formatted data look like:

``` bash
HWI-ST330:304:H045HADXX:2:2212:4026:94975	4	*	0	0	*	*	0	0	GTTTGGTGATGAACTGGGGGTCATAGCCATTGGGGCCCTGCTTATACAAGGAGTTGTAGGCGAGCAGGCGCTCTAGCAGAGAGTAGCCAAGTCCATGCTT	@@@DDDDDDHDDD<EB>E@CC?F?D>3?3?F**067?6C#############################################################	NH:i:0	HI:i:0	AS:i:18	nM:i:2	uT:A:1
HWI-ST330:304:H045HADXX:2:2212:6059:95943	4	*	0	0	*	*	0	0	CCCGACTCCCCGGGTCCCAGGAGGTCCCAGAGGACAGTGGGTGGGAAGTAACCCACAATGCTGGTCTTACAGTAGATGTGGAGTTCATAGCTTTCACCGG	??;B?)@AFDDFD128??CCDGD0)?99?*8?(?##################################################################	NH:i:0	HI:i:0	AS:i:21	nM:i:7	uT:A:1
HWI-ST330:304:H045HADXX:2:2212:7648:96498	4	*	0	0	*	*	0	0	GCCTCATCGATGAAGATGTGTGTGAAGTGATCGATGGGAAACTGGGCTGACACCAACCTGCTGGCAGCGATGAGGGTGGCAATTAAGACCCGATATTGCT	;@@DDDADDFF<D<:AA<A@EB3A43<2<?+CF@;?D3D*):?*?FG9:9B3818=BF##########################################	NH:i:0	HI:i:0	AS:i:53	nM:i:3	uT:A:1
HWI-ST330:304:H045HADXX:2:2212:15028:97450	4	*	0	0	*	*	0	0	GGGGCACCCGGTATATCCCCCGTTTTCTTTACAGAACTCCAGGAACGTTTTCCAGTCTGGGTCGTGGCCTAGGAGGAGGGGGTTGCCCACTACGATGAGC	;<@DFFDDFD>FDGHG@CFF6?GHGH4B*:*??9*98B##############################################################	NH:i:0	HI:i:0	AS:i:38	nM:i:3	uT:A:1
HWI-ST330:304:H045HADXX:2:2212:18192:97569	4	*	0	0	*	*	0	0	CCAGATTTGCGATCTTCATTCCATACAGCATGGAGGAGAAGCCAGGGGCAGGGGTCCCATAGCTGGGCTTGAAGTCGTGGTTGTAAATCGTCCGCAGCCG	=@@BDFFFF<FF8BEFE<FHGCFHE?@<++1:8:8CB13BBD;008(('70/';54?9=.7.;@?@==35>B@###########################	NH:i:0	HI:i:0	AS:i:57	nM:i:3	uT:A:1
HWI-ST330:304:H045HADXX:2:2212:2976:98041	4	*	0	0	*	*	0	0	CACTTGTAAGGGCAGGCCCCCTTCACCCTCCCGCTCCTGGGGGAGATTCTGGGGACGGCGGGGGCCCGAGGGTGGGGGGGTGAGTTTCCTCCGACCTTTG	?;;DDD:BD4CDBEE1EBACEEEE3:CCCDDD@6).0/B=B@##########################################################	NH:i:0	HI:i:0	AS:i:36	nM:i:4	uT:A:1
HWI-ST330:304:H045HADXX:2:2212:8908:98487	4	*	0	0	*	*	0	0	CCGGCTGCCGAGGGTCTCCTGCCAGTGCCAGCTGACATCCGGCATTGACCGGTTACTTGCCAACCATCAGTCCTGCTGTGTCCACCACACTCTCAGGCAC	@<@AA@<12@?FA8CEFE3CFE3**:0000B69B##################################################################	NH:i:0	HI:i:0	AS:i:30	nM:i:7	uT:A:1
HWI-ST330:304:H045HADXX:2:2212:11799:98450	4	*	0	0	*	*	0	0	CGCTCTCCACATGCTCCCTCTGATGCCAGACTTTTTTTTTTTTTTTTTTAAGGTCTTATCATTTTAATCGGTTTTAATGATTATTTCTCATCTGCAAGGT	?@7?B?D3ADFHDHEEEGEFHIBHGG<CGCEHHIGIADFHGFDCB#######################################################	NH:i:0	HI:i:0	AS:i:46	nM:i:2	uT:A:1
HWI-ST330:304:H045HADXX:2:2212:5156:99689	4	*	0	0	*	*	0	0	ATAGCTTTCACCGGGGCCCAGTGGGCAGGGCAGGTCCTGTTCTCCATGGTAGAAGACAAACTGGGGCGGCCAGCACAGTGGGAATAGGTGAGTGAAGGAG	<@@?DAD<F<DDF1)11811:DD9)):0)?;?####################################################################	NH:i:0	HI:i:0	AS:i:22	nM:i:4	uT:A:1
HWI-ST330:304:H045HADXX:2:2212:15671:100398	4	*	0	0	*	*	0	0	TGTAAGGGCAGGCCCCCTTCACCCTCCCGCTCCTGGGGGGGATTCTGGTGACGGCGGGGCCCCGCGGCGGAGGGGCTGGGGTTCCTCCGACCCTGGGGCA	?@@DDDDDFC:DCFGA@?+C3?EHE)?D0@)8?###################################################################	NH:i:0	HI:i:0	AS:i:35	nM:i:9	uT:A:1
```

### Filtering the SAM file

Now we know that we have all of this information for each of the reads -- wouldn't it be useful to summarize and filter based on selected criteria? Suppose we wanted to set a **threshold on mapping quality**. For example, we want to know how many reads aligned with a quality score higher than 30. To do this, we can combine the `view` command with additional flags `-q 30` and `-c` (to count):

```         
$ samtools view -q 30 -c Mov10_oe_1_subsetAligned.sortedByCoord.out.sam

275951
```

> NOTE: mapping quality is the aligner's estimate on how accurate the alignment is. Low MAPQ indicates multimappers. However, it is not always clear what to use as a threhsold, as [this post describes](http://www.acgt.me/blog/2014/12/16/understanding-mapq-scores-in-sam-files-does-37-42).

*How many of reads have a mapping quality of 30 or higher? What are the different quality scores you observe in the file?*

We can also **apply filters to select reads based on where they fall within the `FLAG` categories**. Remember that the bitwise flags are like boolean values. If the flag exists, the statement is true. Similar to when filtering by quality we need to use the `samtools view` command, however this time use the `-F` or `-f` flags.

-   `-f` - to find the reads that agree with the flag statement
-   `-F` - to find the reads that do not agree with the flag statement

```         
## This will tell us how many reads are unmapped
$ samtools view -f 4 -c Mov10_oe_1_subsetAligned.sortedByCoord.out.bam


## This should give us the remaining reads that do not have this flag set (i.e reads that are mapped)
$ samtools view -F 4 -c Mov10_oe_1_subsetAligned.sortedByCoord.out.bam
```

### Indexing the BAM file

To perform some functions (i.e. subsetting, visualization) on the BAM file, an index is required. Think of an index located at the back of a textbook. When you are interested in a particular subject area you look for the keyword in the index and identify the pages that contain the relevant information. Similarly, indexing the BAM file aims to achieve fast retrieval of alignments overlapping a specified region without going through the whole alignment file. In order to index a BAM file, it must first be sorted by the reference ID and then the leftmost coordinate, which can also be done with `samtools`. However, in our case we had included a parameter in our STAR alignment run so we know our BAM files are already sorted.

To index the BAM file we use the `index` command. **This may take a while and actually require you use a SBATCH script, so you only need to do this when you plan to do the bonus visualization step.**

``` bash
$ samtools index Mov10_oe_1_subsetAligned.sortedByCoord.out.bam
```

This will create an index in the same directory as the BAM file, which will be identical to the input file in name but with an added extension of `.bai`. By convention, tools that need BAM index expect the index to be named after its respective BAM file but with a `.bai` extension.

------------------------------------------------------------------------

## **ASSIGNMENT:**

The STAR log file for `Mov10_oe_1` indicated that there were a certain number of reads mapping to multiple locations. When this happens, one of these alignments is considered primary and all the other alignments have the secondary alignment flag set in the SAM records. **Use `samtools` and your knowledge of [bitwise flags](https://github.com/hbc/NGS_Data_Analysis_Course/blob/master/sessionII/lessons/03_alignment_quality.md#bitwise-flags-explained) to find count how many secondary reads there are for `Mov10_oe_1`.**

In a text file in your `/results/STAR directory` called `secondary_reads.txt`, report the number you found above AND the command you used calculate it.

------------------------------------------------------------------------

## BONUS: Visualization

Another method for assessing the quality of your alignment is to visualize the alignment using a genome browser. One recommended piece of software is the [Integrative Genomics Viewer (IGV)](https://www.broadinstitute.org/igv/) from the Broad Institute. IGV is an interactive tool which allows exploration of large, integrated genomic datasets. It supports a wide variety of data types, including array-based and next-generation sequence data, and genomic annotations, which facilitates invaluable comparisons.

To view these files you can either

1.  Download the IGV program and open it on your local computer. This should NOT require help from NICHD IT. You will then need to transfer the following files from your `/results/STAR` directory to your local computer using either `scp` or a locally mounted server:

`Mov10_oe_1_Aligned.sortedByCoord.out.bam`, `Mov10_oe_1_Aligned.sortedByCoord.out.bam.bai`

2.  Use Biowulf's [HPC On Demand](https://hpcondemand.nih.gov/) tool which provides convenient web interfaces to your interactive Biowulf applications, while allowing you access to all of your files on Biowulf. We will be making use of this for RStudio in weeks 5-8.

    You may need to copy over your `bam` and `bai` files over to your `/data/$USER` area on Biowulf for this to work.

### Visualize

-   Start [IGV](https://www.broadinstitute.org/software/igv/download), either locally or on HPC On Demand
-   Load the Human genome (hg19) into IGV using the dropdown menu at the top left of your screen. *Note: there is also an option to "Load Genomes from File..." under the "Genomes" pull-down menu - this is useful when working with non-model organisms*
-   Load the .bam file using the **"Load from File..."** option under the **"File"** pull-down menu. *IGV requires the .bai file to be in the same location as the .bam file that is loaded into IGV, but there is no direct use for that file.*

![IGV screenshot](../img/IGV_mov10.png)

------------------------------------------------------------------------

**Exercise**

Take a look at a few other genes by typing into the search bar. For example, PPM1J and PTPN22. How do these genes look?

------------------------------------------------------------------------

*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
