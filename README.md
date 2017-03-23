# immunoSeqR

## Description

An R package for analyzing TCR sequencing data from Adaptive Biotechnology's ImmunoSeq platform.
It is primarily a set of tools for importing the data into R, but it also includes a number of
scripts that can be used as examples for analyzing complex experiments.

The data is imported into a single `data.frame` containing all of the samples. Rows represent
unique T cell clones and columns represent samples. When clones are not shared between samples,
the table is completed by adding zeros. This results in a larger than necessary table, but saves
time by not having to join samples on the fly. 

`scripts/example.R` contains a full example using some small sample files.

## Importing Data 
Data is imported using `iseqr_merge()`, which takes a path to the directory containing Adaptive
tsv files as its argument.

```R

path <- '/data/ex_tsv/' # change to location of example files
setwd(path)
all_files <- list.files(pattern=".tsv")

# construct the dataset with iseqr_merge
ds <- iseqr_merge(all_files)
``` 

## Collapsing Data 
Many nucleic acid sequences can encode the same CDR3, so many analyses may require aggregated
data (data in which synonymous nucleotide sequences are combined into a single amino acid level
representation). Metrics such as Clonality and Richness should typically be computed from aggregated
data.  
An imported dataset can be aggregated using `iseqr_aggregate()`

```R
# Remove Empty Sequences
length(which(ds$aa==''))
ds <- ds[ds$aa!='',]

# Remove Sequences with stop codons
length(grep('\\*',ds$aa))
ds <- ds[grep('\\*',ds$aa,invert=TRUE),]

# aggregate the data
# this collapses synonymous nucleotide sequences
ds_agg <- iseqr_aggregate(ds,inc_nt=FALSE)
```  
