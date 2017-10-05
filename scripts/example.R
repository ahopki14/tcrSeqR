#set path to adaptive tsv files
path <- '~/Documents/emj/ImmunoseqResults/immunoSeqR/data/ex_tsv/tsv'
setwd(path)
all_files <- list.files(pattern=".tsv")

# construct the dataset with iseqr_merge
ds <- iseqr_merge(all_files)

#load the dictionary
dict <- readRDS('dict.Rds')

#make tcr object
ds <- iseqr_make_tcr(ds,dict)

# aggregate the data
# this collapses synonymous nucleotide sequences
ds_agg <- iseqr_aggregate(ds)


#going forward, use on the aggregated dataset
ds <- ds_agg


###############################################################

# The data object can be filtered using plyr commands
require(dplyr)
filter(ds, type=='PRE')
# or
filter(ds, patient=='8.001')


#metrics can be calculated for the samples
clonality(ds, merge=F)

# and these can be merged directly into the object's metadata
ds <- clonality(ds)
ds <- richness(ds)
ds <- total(ds)
ds$Clonality
ds$Richness

# To calculate overlaps, as well as changes in metrics we will need to define
# the comparisons that need to be made
comps <- list(c('PRE','POST1'), c('PRE','POST3'))

# now we can add the morisita index (similarity) to the object
ds <- iseqr_morisita(ds, comps)

# And finally, compute the fold change in metrics (using the same comparisons)
ds <- delta_stats(ds,comps,'Clonality')
ds <- delta_stats(ds,comps,'Richness')


#now plot whatever you want
iseqr_plot_factor(ds, 'Clonality', 'response', type='PRE')
iseqr_plot_factor(ds, 'Log2.Fold.Change.in.Clonality', 'arm', type='POST3')
iseqr_plot_metrics(ds, 'Total.Sequences', 'ALC', type='PRE')


