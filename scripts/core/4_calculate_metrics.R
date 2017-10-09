library(immunoSeqR)
#set path to where ds_agg is saved
setwd('/home/student/ahopkins/emj/immunoseq_data/adjuvant')
ds <- readRDS('ds_agg.Rds')

#calculate basics 
ds <- clonality(ds)
ds <- richness(ds)
ds <- total(ds)

#specify comparisons
comps <- list(c('PRE','POST1'), c('PRE','POST3'))

# now we can add the morisita index (similarity) to the object
ds <- iseqr_morisita(ds, comps)

# And finally, compute the fold change in metrics (using the same comparisons)
ds <- delta_stats(ds,comps,'Clonality')
ds <- delta_stats(ds,comps,'Richness')

#save
saveRDS(ds, file='ds_agg_final.Rds')
