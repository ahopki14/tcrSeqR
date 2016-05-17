sum.ds <- list()
sum.ds$nsamp <- ncol(ds)-2
sum.ds$ntypes <- length(levels(dict$type))
sum.ds$unique.tcr <- nrow(ds)
sum.ds$total.tcr <- sum(ds$syn)
sum.ds$avg.syn <- mean(ds$syn)
sum.ds$max.syn <- max(ds$syn)
sum.ds$max.syn.aa <- as.character(ds$aa[which.max(ds$syn)])

