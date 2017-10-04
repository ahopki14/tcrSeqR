iseqr_make_tcr <- function(ds, dict){
	#reoder the dictionary to match the ds
	dict <- iseqr_order(dict, ds, reorder=T)
	#make an ordering variable
	dict$ord <- seq(nrow(dict))
	#get metadata and data column locations
	w_m <- grep('aa|nt|syn', names(ds))
	w_d <- grep('aa|nt|syn', names(ds), invert=T)
	# make the SE object with the data
	if(!any(grepl('syn',names(ds)))){
		#if nt data
		ds_se <- SummarizedExperiment(assays=list(tcr_nt=as.matrix(ds[,w_d], rownames=F)))
				}else{
		#if agg data
		ds_se <- SummarizedExperiment(assays=list(tcr=as.matrix(ds[,w_d], rownames=F)))
	}
	# add some metadata
	rownames(ds_se) <- ds$aa
	rowData(ds_se) <- DataFrame(ds[,w_m])
	stopifnot(all(colnames(ds_se)==dict$fn))
	colData(ds_se) <- DataFrame(dict)
	colnames(ds_se) <- colnames(ds[,w_d])
	#return the SE object
	tcr(ds_se)
}
