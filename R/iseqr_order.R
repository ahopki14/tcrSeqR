iseqr_order <- function(plot_ds,ds, reorder=T){
	ord <- match(names(ds), plot_ds$fn)
	ord <- ord[!is.na(ord)]
	w <- grep('syn|aa|nt', names(ds))
	stopifnot(names(ds)[-w] == plot_ds$fn[ord])
	if(reorder){
		plot_ds[ord,]
	}else{ord}
}	
