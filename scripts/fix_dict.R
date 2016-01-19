colnames(ds) <- gsub('sample','',colnames(ds))
dict$fn <- gsub('\\.','',dict$fn)

ord <- match(dict$fn,names(ds))
ord[1] <- 1
ord[2] <- 2
names(ds)[ord]==dict$fn
ds <- ds[,ord]
names(ds)==dict$fn

# save
