library(immunoSeqR)
#set path to where ds_nt is saved
setwd('/home/student/ahopkins/emj/immunoseq_data/adjuvant')
ds_nt <- readRDS('ds_nt.Rds')

# Aggregate the data
ds_agg <- iseqr_aggregate(ds_nt)

#save
saveRDS(ds_agg,file='ds_agg.Rds')
