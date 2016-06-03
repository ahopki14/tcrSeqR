
bmm <- mm
bolm <- olm
path <- '~/Documents/emj/ImmunoseqResults/mega/'
for(a in levels(dict$experiment)){
	w <- which(dict$experiment==a)
#	tds <- ds[,c(w,1,2)] # put syn and aa on the end, then dictionary does not need to be padded
#	tds <- iseqr_clean_ds(tds)
#	tdict <- dict[w,] # will result in an unpadded dictionary, but will work with new ds
#	tdict <- refactor(tdict)
#	saveRDS(tds,paste0(path,'ds_',a,'.Rds'))
#	saveRDS(tdict,paste0(path,'dict_',a,'.Rds'))
	tstats <- stats[w,]
	olm <- bolm[w,w]
	mm <- bmm[w,w]
	#saveRDS(tstats,paste0(path,'stats_',a,'.Rds'))
	save(olm,mm,file=paste0(path,'olm_',a,'.Rda'))
}


#testing
load('olm_Neoadjuvant.Rda')
stats <- readRDS('stats_Neoadjuvant.Rds')
dict <- readRDS('dict_Neoadjuvant.Rds')

all(dict$fn == stats$fn)
all(dict$fn == rownames(olm))
all(dict$fn == rownames(mm))
