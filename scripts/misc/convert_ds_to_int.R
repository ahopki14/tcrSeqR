setwd('~/Documents/emj/ImmunoseqResults/mega/')
files <- list.files(pattern="ds.*Rds")
for(b in files){
	ds <- readRDS(b)
	if(class(ds[,1])!="numeric"){stop("ds out of order")}
	for(a in seq(ncol(ds)-2)){
		ds[,a] <- as.integer(ds[,a])
	}
	saveRDS(ds,file=b)
}
