#colnames(ds) <- gsub('sample','',colnames(ds))

ord <- match(dict$fn,names(ds))
ord[1] <- 1
ord[2] <- 2
names(ds)[ord]==dict$fn
ds <- ds[,ord]
