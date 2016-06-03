all_files <- list.files(pattern=".tsv")
ds <- iseqr_merge(all_files)
# can load a 2 sample data set from example_data.Rda
ds <- ds[ds$aa!='',] # not necessary as aggregate does this first...
ds <- ds[!grepl('*',ds$aa,fixed=TRUE),] # remove any with stop codons
ds_agg <- iseqr_aggregate(ds)
save(ds,file='ds_agg.Rda')
#
# find clones with high synonimity
tmp <- ds_agg[ds_agg$syn>50,]

# Make descriptors from ds
describe <- function(ds){
nm <- names(ds)
is_meta <-  nm=='syn' | nm=='aa'
til <- grepl('TIL',nm)
tumor <- grepl('PDAC',nm)
pre <- grepl('Pre',nm)
post <- grepl('V1D',nm)
type <- nm
type[is_meta] <- 'NA'
type[til] <- 'TIL'
type[tumor] <- 'Tumor'
type[pre] <- 'Pre'
type[post] <- 'Post'
type
}
