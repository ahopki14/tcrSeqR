#colnames(ds) <- gsub('sample','',colnames(ds))
ord <- match(colnames(ds),gsub('\\.','',as.character(dict$fn)))
ord[1] <- 1
ord[2] <- 2
ds <- ds[,ord]
