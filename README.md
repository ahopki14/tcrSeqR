# immunoSeqR
An R package for analyzing TCR sequencing data from Adaptive Biotechnology's ImmunoSeq platform.
It is primarily a set of tools for importing the data into R, but it also includes a number of
scripts that can be used as examples for analyzing complex experiments.

The data is imported into a single `data.frame` containing all of the samples. Rows represent
unique T cell clones and columns represent samples. When clones are not shared between samples,
the table is completed by adding zeros. This results in a larger than necessary table, but saves
time by not having to join samples on the fly. 

Get started by checking out the example.R script in scripts/
