library(immunoSeqR)
load('~/Documents/emj/ImmunoseqResults/sampleExport.2014-07-31_10-10-24/ds_agg_reorder.Rda')
load('~/Documents/emj/ImmunoseqResults/sampleExport.2014-07-31_10-10-24/dict_reorder.Rda')
olScan <- function(x,y,byn=10){
	ds <- data.frame(x=x,y=y)
	ds <- ds[order(-ds$x), ]
	byv <- seq(10,nrow(ds),byn) 
	ol <- mapply(byv,FUN=function(n){overlap(ds$x[1:n],ds$y[1:n])})
	ol
}
ds <- ds_agg_reorder
dict <- dict_reorder

out <- data.frame()
for(a in seq(length(levels(dict$patient)))){
	name <- levels(dict$patient)[a]
	resp <- as.character(dict$response[which(dict$patient==name)][1])
	samples <- c( # include the aa data
		which(dict$patient==name & dict$type=='Pre'),
		which(dict$patient==name & dict$type=='Post')
	)
tmp <- olScan(ds[,samples[1]],ds[,samples[2]])
save(tmp,file=paste0('~/Documents/emj/ImmunoseqResults/olscan/',name,'.Rda'))
#out[,a] <- tmp
#names(out)[a] <- levels[dict$patient][a]
}
