---
title: "AnnotationDbi"
author: "Harvard HPC Staff, Adapted by Sally Chang at NICHD"
date: "Last Modified March 2025"
---

## AnnotationDbi

AnnotationDbi is an R package that provides an interface for connecting and querying various annotation databases using SQLite data storage (a convenient and performant flat-file database format). There are AnnotationDbi packages that contain information on gene IDs, Gene Ontology terms, and more. There is helpful [documentation](https://bioconductor.org/packages/release/bioc/vignettes/AnnotationDbi/inst/doc/IntroToAnnotationPackages.pdf) available to reference when extracting data from any of these databases.

AnnotationDbi is just the interface; you need to get an appropriate annotation package that you can use with AnnotationDbi.

There are a plethora of organism-specific *orgDb* packages, such as `org.Hs.eg.db` for human and `org.Mm.eg.db` for mouse. A list of organism databases can be found [here](https://www.bioconductor.org/packages/release/BiocViews.html#___OrgDb). These databases are best for converting gene IDs or obtaining GO information for current genome builds, but not for older genome builds. These packages provide the current builds corresponding to the release date of the package, and update every 6 months. If a package is not available for your organism of interest, advanced users can create their own witn [AnnotationForge](https://bioconductor.org/packages/release/bioc/html/AnnotationForge.html). 

Humans and common model organisms have packages that can be installed, like `org.Hs.eg.db` for human. However, only a subset of available annotations are available as packages. For anything else, use the [AnnotationHub](https://bioconductor.org/packages/release/bioc/html/AnnotationHub.html) package. Model organism annotations are also available through AnnotationHub, so this is the most generic approach.


### org.Hs.eg.db

For convenience, here we will use the already-installed `org.Hs.eg.db` OrgDb:

``` r
# Load libraries
library(org.Hs.eg.db)
library(AnnotationDbi)

# Check object metadata
org.Hs.eg.db
```

We can see the metadata for the database by just typing the name of the database, including the species, last updates for the different source information, and the source urls. Note the KEGG data from this database was last updated in 2011. KEGG is under a subscription model so BioConductor can not redistribute the content since this earlier freeze date.

```         
OrgDb object:
| DBSCHEMAVERSION: 2.1
| Db type: OrgDb
| Supporting package: AnnotationDbi
| DBSCHEMA: HUMAN_DB
| ORGANISM: Homo sapiens
| SPECIES: Human
| EGSOURCEDATE: 2024-Sep20
| EGSOURCENAME: Entrez Gene
| EGSOURCEURL: ftp://ftp.ncbi.nlm.nih.gov/gene/DATA
| CENTRALID: EG
| TAXID: 9606
| GOSOURCENAME: 
| GOSOURCEURL: 
| GOSOURCEDATE: 
| GOEGSOURCEDATE: 2024-Sep20
| GOEGSOURCENAME: Entrez Gene
| GOEGSOURCEURL: ftp://ftp.ncbi.nlm.nih.gov/gene/DATA
| KEGGSOURCENAME: KEGG GENOME
| KEGGSOURCEURL: ftp://ftp.genome.jp/pub/kegg/genomes
| KEGGSOURCEDATE: 2011-Mar15
| GPSOURCENAME: UCSC Genome Bioinformatics (Homo sapiens)
| GPSOURCEURL: ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/database
| GPSOURCEDATE: 2024-Sep22
| ENSOURCEDATE: 2024-May14
| ENSOURCENAME: Ensembl
| ENSOURCEURL: ftp://ftp.ensembl.org/pub/current_fasta
| UPSOURCENAME: Uniprot
| UPSOURCEURL: http://www.UniProt.org/
| UPSOURCEDATE: Mon Sep 23 15:46:45 2024
```

We can extract information from this database using *AnnotationDbi* with the methods: `columns`, `keys`, `keytypes`, and `select`. Here were are using our `org.Hs.eg.db` database to acquire information, but the same methods work for the *TxDb*, *Go.db*, *EnsDb*, and *BioMart* annotations.

Let's inspect what we have available in the orgdb:

```r
columns(org.Hs.eg.db)
```

```
[1] "ACCNUM"       "ALIAS"        "ENSEMBL"      "ENSEMBLPROT"  "ENSEMBLTRANS"
[6] "ENTREZID"     "ENZYME"       "EVIDENCE"     "EVIDENCEALL"  "GENENAME"
[11] "GENETYPE"     "GO"           "GOALL"        "IPI"          "MAP"
[16] "OMIM"         "ONTOLOGY"     "ONTOLOGYALL"  "PATH"         "PFAM"
[21] "PMID"         "PROSITE"      "REFSEQ"       "SYMBOL"       "UCSCKG"
```


Ideally, all databases would use the same gene identifier, but this is unfortunately not the case. Recall that our featureCounts file that we imported into R has Ensembl identifiers. Let's check the orgdb's Ensemble identifiers:

```r
keys(org.Hs.eg.db, "ENSEMBL") %>% head()
```

```
[1] "ENSG00000121410" "ENSG00000175899" "ENSG00000291190" "ENSG00000171428"
[5] "ENSG00000156006" "ENSG00000196136"
```

Recall that our Ensembl IDs had a dotted version number, which will not match anything in this orgdb. So we need to clean those up, e.g., by turning `ENSG00000000003.16` into just `ENSG00000000003`. This uses a [regular expression](https://r4ds.had.co.nz/strings.html#matching-patterns-with-regular-expressions) for the pattern matching:

``` r
# use gsub to replace column with stripped gene symbols
res_tableOE_df$ensgene <- gsub("\\..*","",res_tableOE_df$ensgene)
```


Now we can use that cleaned set of identifiers to select items from the orgdb.


``` r
# Return the Ensembl IDs for a set of genes
annotations_orgDb <- AnnotationDbi::select(org.Hs.eg.db, # database
                                     keys = res_tableOE_df$ensgene,  # data to use for retrieval
                                     columns = c("SYMBOL", "ENTREZID","GENENAME"), # information to retreive for given data
                                     keytype = "ENSEMBL") # type of data given in 'keys' argument
```

We started from at about 79k in our results table, but the results from this are quite a bit longer (83k). This is because there are Ensembl IDs that have multiple entries in other columns. Let's inspect that:

```r
multi <- table(annotations_orgDb$ENSEMBL) %>%
  sort() %>%
  tail(25)
multi
```

```
ENSG00000258992 ENSG00000199270 ENSG00000199334 ENSG00000199352 ENSG00000199396 
             36              91              91              91              91 
ENSG00000199910 ENSG00000200343 ENSG00000200370 ENSG00000200381 ENSG00000200624 
             91              91              91              91              91 
ENSG00000201355 ENSG00000201588 ENSG00000201925 ENSG00000202257 ENSG00000202521 
             91              91              91              91              91 
ENSG00000202526 ENSG00000199337 ENSG00000273730 ENSG00000275757 ENSG00000274917 
             91              92             207             207             208 
ENSG00000275215 ENSG00000276700 ENSG00000277739 ENSG00000278189 ENSG00000278233 
            208             208             208             208             208 
```

What genes are those!?

```r
annotations_orgDb %>%  filter(ENSEMBL %in% names(multi)) %>% head
```

```
          ENSEMBL       SYMBOL  ENTREZID             GENENAME
1 ENSG00000199270      RNA5S12 100169763 RNA, 5S ribosomal 12
2 ENSG00000199270 LOC124905422 124905422     5S ribosomal RNA
3 ENSG00000199270 LOC124905424 124905424     5S ribosomal RNA
4 ENSG00000199270 LOC124905425 124905425     5S ribosomal RNA
5 ENSG00000199270 LOC124905426 124905426     5S ribosomal RNA
6 ENSG00000199270 LOC124905427 124905427     5S ribosomal RNA
```

So it looks like a single Ensembl ID could refer to multiple different symbols, here for 5S ribosomal RNA.

Let's take a peek to see if we actually returned annotations for each individual Ensembl gene ID that went in to the query:

``` r
length(which(is.na(annotations_orgDb$SYMBOL)))
```

Looks like more than half of the input genes did not return any annotations. This is because the OrgDb family of database are primarily based on mapping using Entrez Gene identifiers. If you look at some of the Ensembl IDs from our query that returned NA, these map to pseudogenes (i.e [ENSG00000265439](https://useast.ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000265439;r=6:44209766-44210063;t=ENST00000580735)) or non-coding RNAs (i.e. [ENSG00000265425](http://useast.ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000265425;r=18:68427030-68436918;t=ENST00000577835)). The difference is due to the fact that each database implements different computational approaches for generating the gene builds. And some databases (notably, Gene Ontology) only have functional information on proteins, which in turn only come from coding genes. Our intention is to use Gene Ontology, so while it's interesting to know about non-coding genes that are differentially expressed, let's get rid of those genes from our data that do not have corresponding entries in the orgdb:

``` r
# Determine the indices for the non-NA genes
non_na_idx <- which(is.na(annotations_orgDb$SYMBOL) == FALSE)

# Return only the genes with annotations using indices
annotations_orgDb <- annotations_orgDb[non_na_idx, ]
```

You may have also noted the *warning* returned: *'select()' returned 1:many mapping between keys and columns*. This is always going to happen with converting between different gene IDs (i.e. one geneID can map to more than one identifier in another databse). This is an unfortunate practical issue with gene nomenclature. Unless we would like to keep multiple mappings for a single gene, then we probably want to de-duplicate our data before using it.

``` r
# Determine the indices for the non-duplicated genes
non_duplicates_idx <- which(duplicated(annotations_orgDb$SYMBOL) == FALSE)

# Return only the non-duplicated genes using indices
annotations_orgDb <- annotations_orgDb[non_duplicates_idx, ]
```

### EnsDb.Hsapiens.v86

To generate the Ensembl annotations, the *EnsDb* database can also be easily queried using AnnotationDbi. You will need to decide the release of Ensembl you would like to query. We know that our data is for GRCh38, and the most current *EnsDb* release for GRCh38 in Bioconductor is release 86, so we can install this database. All Ensembl releases are listed [here](http://useast.ensembl.org/info/website/archives/index.html). **NOTE: this is not the most current release of GRCh38 in the Ensembl database, but it's as current as we can obtain through AnnotationDbi.**

Since we are using *AnnotationDbi* to query the database, we can use the same functions that we used previously:

``` r
# Load the library
library(EnsDb.Hsapiens.v86)

# Check object metadata
EnsDb.Hsapiens.v86

# Explore the fields that can be used as keys
keytypes(EnsDb.Hsapiens.v86)

# Explore columns of data that can be used for other analyses
columns(EnsDb.Hapiens.v86)
```

Now we can return all gene IDs for our gene list:

``` r
# Return the Ensembl IDs for a set of genes
annotations_edb <- AnnotationDbi::select(EnsDb.Hsapiens.v86,
                                           keys = res_tableOE_tb$gene,
                                           columns = c("SYMBOL", "ENTREZID","GENEBIOTYPE"),
                                           keytype = "GENEID")
```

We can check for NA entries, and find that there are none:

``` r
length(which(is.na(annotations_edb$SYMBOL) == FALSE))
```

Then we can again deduplicate, to remove the gene symbols which appear more than once:

``` r
# Determine the indices for the non-duplicated genes
non_duplicates_idx <- which(duplicated(annotations_edb$SYMBOL) == FALSE)

# Return only the non-duplicated genes using indices
annotations_edb <- annotations_edb[non_duplicates_idx, ]
```

> **NOTE:** In this case we used the same build but a slightly older release, and we found little discrepancy. If your analysis was conducted using an older genome build (i.e hg19), but used a newer build for annotation some genes may be found to be not annotated (NA). Some of the genes have changed names in between versions (due to updates and patches), so may not be present in the newer version of the database.

------------------------------------------------------------------------

*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
