library(immunoSeqR)
library(dplyr)
library(ggplot2)


#set path to adaptive tsv files
path <- '~/Documents/emj/ImmunoseqResults/immunoSeqR/data/ex_tsv/tsv'
setwd(path)
all_files <- list.files(pattern=".tsv")

# construct the dataset with iseqr_merge
ds <- iseqr_merge(all_files)



# Remove Empty Sequences
length(which(ds$aa==''))
ds <- ds[ds$aa!='',]

# Remove Sequences with stop codons
length(grep('\\*',ds$aa))
ds <- ds[grep('\\*',ds$aa,invert=TRUE),]

# aggregate the data
# this collapses synonymous nucleotide sequences
ds_agg <- iseqr_aggregate(ds,inc_nt=FALSE)

#moving the aa and nt columns to the end makes life easier
ds <- ds[,c(3:66,1,2)]
ds_agg <- ds_agg[,c(3:66,1,2)]

#load the dictionary
dict <- readRDS('dict.Rds')
#put it in the correct order
dict <- iseqr_order(dict, ds_agg)

#going forward, use on the aggregated dataset
ds <- ds_agg

# rows of dictionary correspond to cols of ds
names(ds)
names(ds) == dict$fn

# so you can do this:
w <- which(dict$type=='PRE') # restrict to Pre-vac sample types
head(ds[,w])

# or this:
w <- which(dict$patient=='8.001') # restrict to patient 1
head(ds[,w])

# Clonality of first sample
clonality(ds[,1])

# quickly calculate for all samples
sapply(ds[,-c(65,66)],clonality)


# What about richness?
length(which(ds[,1]>0))


# overlap
overlap(ds[,1],ds[,2])
morisita(ds[,1],ds[,2])

# make plot_ds with script using
comps <- list(c('PRE','POST1'), c('PRE','POST3'))

#now plot whatever you want
iseqr_plot_factor(plot_ds, 'Clonality', 'response', type='PRE')
iseqr_plot_factor(plot_ds, 'Log2.Fold.Change.in.Clonality', 'arm', type='POST3')
iseqr_plot_metrics(plot_ds, 'Total.Sequences', 'ALC', type='PRE')


