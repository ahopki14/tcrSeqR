#library(immunoSeqR)
setwd('/home/student/ahopkins/emj/immunoseq_data/adjuvant')
load('ds.Rda')
ds_agg <- iseqr_aggregate(ds,inc_nt=FALSE)
ds <- ds_agg

# Remove Empty Sequences
length(which(ds$aa==''))
ds <- ds[ds$aa!='',]

# Remove Sequences with stop codons
length(grep('\\*',ds$aa))
ds <- ds[grep('\\*',ds$aa,invert=TRUE),]

# Re-factor
ds$aa <- as.factor(as.character(ds$aa))


save(ds,file='ds_agg.Rda')
