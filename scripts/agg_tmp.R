# A minimal example data set
ds <- data.frame(nt= c('j','m','h','o','k','l','i','n'),aa=c('c','c','b','','a','d','b','b'),data=c(4,5,14,0,7,3,8,4),data2=c(6,2,1,5,6,1,1,9))
# or load the 2 sample set
#load('example_data.Rda')
ds <- ds[ds$aa!='',] # restrict to productive
ds$aa <- as.factor(as.character(ds$aa)) # necessary to remove un-productive seq from factor levels
#
# 
ind <- split(seq_len(nrow(ds)),ds$aa)
syn <- sapply(ind,length)
#
dc <- names(ds)[names(ds)!='nt' & names(ds)!='aa']
#
sum_rows <- function(x){sapply(ds[x,dc],FUN=sum)}
a <- sapply(ind,FUN=sum_rows)
b <- data.frame(t(a))

paste_names <- function(x){toString(ds[x,'nt'])}
nt <- sapply(ind,FUN=paste_names)
nt <- data.frame(n)

rn <- data.frame(aa=rownames(b))
rownames(b) <- seq_len(nrow(b))
syn <- data.frame(syn)
ds_out <- cbind(syn,nt,rn,b)
rownames(ds_out) <- seq_len(nrow(ds_out))
