## not tested!
ntchar <-ds$nt ## assuming this is a character vector
ntfactor <- factor(ntchar, levels=unique(ntchar))
## if ntfactor is cumbersome because of size, coerce to integer 
ntint <- as.integer(ntfactor)
ntfactor2 <- factor(ntint, levels=unique(ntint))

M <- as.matrix(ds$data, ds$data2)  

indices <- split(seq_len(nrow(M)), ntfactor2)
L <- sapply(indices, length) # this is wrong
## For sequences that have no match, there is nothing to do
dat1 <- ds[indices[L == 1], ]

## Sequences that have one or more matches
indices2 <- indices[L > 1]
sum1 <- sapply(indices2, function(i, x) sum(x[i]), x=M[,1])
sum2 <- sapply(indices2, function(i, x) sum(x[i]), x=M[,2])

## could be slow?  Biostrings?
catstring <- sapply(indices2, function(i, x, y) paste(x[i], y[i], sep=“,”), x=ntchar, y=ntchar)

dat2 <- data.frame(nt=catstring, data=sum1, data2=sum2)
newdat <- rbind(dat1, dat2)
