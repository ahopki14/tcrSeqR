all_files <- list.files(pattern=".tsv")
ds <- iseqr_merge(all_files)
# can load a 2 sample data set from example_data.Rda
ds <- ds[ds$aa!='',] # not necessary as aggregate does this first...
ds_agg <- iseqr_aggregate(ds)
save(ds,file='exp1_iseqr_agg.Rda')

