library(immunoSeqR)
#set path to where ds_nt is saved
setwd('/home/student/ahopkins/emj/immunoseq_data/adjuvant')
ds <- readRDS('ds.Rds')

#Load the dictionary
dict <- readRDS('dict.Rds')

#make tcr object
ds_nt <- iseqr_make_tcr(ds,dict)

#Calculate the expanded clones
comps <- list(c('PRE','POST1'), c('PRE','POST3'))
ds_nt <- iseqr_exp_cl(ds_nt, comps, output='tcr')

saveRDS(ds_nt, file='ds_nt.Rds')

# Aggregate the data
ds_agg <- iseqr_aggregate(ds_nt)

#save
saveRDS(ds_agg,file='ds_agg.Rds')
