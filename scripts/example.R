#set path to adaptive tsv files
path <- '~/Documents/emj/ImmunoseqResults/immunoSeqR/data/ex_tsv'
setwd(path)
all_files <- list.files(pattern=".tsv")

# construct the dataset with iseqr_merge
ds <- iseqr_merge(all_files)

#load the dictionary
# This can be done with read.csv, as long as factors are created rather than
# characters for metadata of interest
# here, we can use a pre-made dictionary 
dict <- readRDS('dict.Rds')

#make tcr object
ds_nt <- iseqr_make_tcr(ds,dict)

#This is the optimal time to calculate expanded clones (before you aggregate the
#data) This requires a list of comparisons
comps <- list(c('PRE','POST1'), c('PRE','POST3'))
ds_nt <- iseqr_exp_cl(ds_nt, comps, output='tcr')

# aggregate the data
# this collapses synonymous nucleotide sequences
ds_agg <- iseqr_aggregate(ds_nt)


#going forward, use on the aggregated dataset
ds <- ds_agg


###############################################################

# The data object can be filtered using plyr commands
filter(ds, type=='PRE')
# or
filter(ds, patient=='1')


#metrics can be calculated for the samples
clonality(ds, merge=F)

# and these can be merged directly into the object's metadata
ds <- clonality(ds)
ds <- richness(ds)
ds <- total(ds)
ds$Clonality
ds$Richness

# Like expanded clones (see above), overlaps and fold changes require a list of
# comparisons
comps <- list(c('PRE','POST1'), c('PRE','POST3'))

# now we can add the morisita index (similarity) to the object
ds <- iseqr_morisita(ds, comps)

# And finally, compute the fold change in metrics (using the same comparisons)
ds <- delta_stats(ds,comps,'Clonality')
ds <- delta_stats(ds,comps,'Richness')


#now plot whatever you want
iseqr_plot_factor(ds, metric='Clonality', by='response', type='PRE')
iseqr_plot_factor(ds, 'Log2.Fold.Change.in.Clonality', 'arm', type='POST3')
iseqr_plot_factor(ds, 'Total.Sequences', 'sex', type='PRE')
iseqr_plot_metrics(ds, 'Number.of.Expanded.Clones', 'age', type='POST1')


