load('~/Documents/emj/ImmunoseqResults/sampleExport.2014-07-31_10-10-24/ds_agg_reorder.Rda')
load('~/Documents/emj/ImmunoseqResults/sampleExport.2014-07-31_10-10-24/dict_reorder.Rda')
out <- data.frame()
for(a in seq(length(levels(dict$patient)))){
	name <- levels(dict$patient)[a]
	resp <- as.character(dict$response[which(dict$patient==name)][1])
	samples <- c( # include the aa data
		which(dict$patient==name & dict$type=='Pre'),
		which(dict$patient==name & dict$type=='Post')
	)
tmp <- olScan(ds[,samples[1]],ds[,samples[2]])
out[,a] <- tmp
names(out)[a] <- levels[dict$patient][a]
}
