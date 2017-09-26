# immunoSeqR

## Description

An R package for analyzing TCR sequencing data from Adaptive Biotechnology's ImmunoSeq platform.
It is primarily a set of tools for importing the data into R, but it also includes a number of
scripts that can be used as examples for analyzing complex experiments.

The data is imported into a single `tcr` object, containing all of the samples and associated metadata.
the `tcr` class extends `SummarizedExperiment` and has a number of custom
methods. The data is represented as a table with rows representing
unique T cell clones and columns representing samples. When clones are not shared between samples,
the table is completed by adding zeros. This results in a larger than necessary table, but saves
time by not having to join samples on the fly. 

Metadata for the samples, as well as metrics summarizing them are stored in the
colData slot of the `tcr` object. 

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

## Building a Metadata Dictionary
`immunoSeqR` can help make plots of various metrics, but it needs to be aware of the metadata
associated with each sample. This is accomplished with a dictionary, which is simply a
`data.frame` in which each row corresponds to a column in the dataset. A simple example might
look like 

fn | patient | type | response
---|---------|------|---------- 
  samplept1post|       1|  Post|        R
   samplept1pre|       1|   Pre|        R
 samplept1tumor|       1| Tumor|        R
  samplept2post|       2|  Post|       NR
   samplept2pre|       2|   Pre|       NR
 samplept2tumor|       2| Tumor|       NR

In this example, `fn` refers to the original filename of the tsv (which is brought in as the
column name in the dataset), `patient` is a patient/subject number, `type` is a sample type (in
this case pre and post treatment as well as tumor) and `response` indicates if a patient was a
responder or non-responder.  

##Creating the `tcr` object
Next, the data is combined with any metadata available to create the `tcr`
object. To match the samples to their metadata, the metadata must contain a
column called `fn` which matches the filenames of the sample tsv file
(specifically, it must match the `colnames` of the `ds` object at this stage.
With this metadata loaded as a `data.frame`, the `tcr` object can be constructed
using

```R
#load the dictionary
dict <- readRDS('dict.Rds')

#make tcr object
ds <- iseqr_make_tcr(ds,dict)
```

## Collapsing Data 
Many nucleic acid sequences can encode the same CDR3, so many analyses may require aggregated
data (data in which synonymous nucleotide sequences are combined into a single amino acid level
representation). Metrics such as Clonality and Richness should typically be computed from aggregated
data.  
This step also removes any sequences with stop codons, or sequences with no
translation. 

An imported dataset can be aggregated using `iseqr_aggregate()`

```R
# aggregate the data
# this collapses synonymous nucleotide sequences
ds_agg <- iseqr_aggregate(ds,inc_nt=FALSE)
```  



## Metrics 
A variety of metrics are available in `immunoSeqR`, including `clonality`, `overlap`, and `morisita`. 
These typically operate on a single sample, and therefore need to be used with `apply` to
calculate the metric for all samples

```R
sapply(ds,clonality)
```

Other metrics include wrappers that calculate the metric in an intelligent way using the
metadata. For example, `iseqr_morisita()` can calculate the Morisita Index for any pairs of
samples in a dataset, and `iseqr_exp_cl()` can identify the T cells which expand after treatment.



