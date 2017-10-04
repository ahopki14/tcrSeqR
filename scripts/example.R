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

# The data object can be filtered using pyr commands
filter(ds, type=='PRE')
# or
filter(ds, patient=='8.001')


#metrics can be calculated for the samples
clonality(ds, merge=F)

# and these can be merged directly into the object's metadata
ds <- clonality(ds)
ds <- richness(ds)
ds$clonality

#############OLD#################


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


