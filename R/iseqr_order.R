#' iseqr_order
#' 
#' Checks that the dictionary and data are in the same order and, if necessary,
#' reorders the dictionary to fit the data. This is used in the iseqr_merge
#' function
#'
#' @param plot_ds A dictionary or plot_ds object
#' @param ds An object created with iseqr_merge
#' @param reorder Logical indicating if the plot_ds should be reordered
#' @return If reorder=T, a plot_ds object in the same order as the ds. If
#' reorder=F, the order in which it should be placed to acheive this.
#' @author Alexander Hopkins
iseqr_order <- function(plot_ds,ds, reorder=T){
	ord <- match(names(ds), plot_ds$fn)
	ord <- ord[!is.na(ord)]
	w <- grep('syn|aa|nt', names(ds))
	stopifnot(names(ds)[-w] == plot_ds$fn[ord])
	if(reorder){
		plot_ds[ord,]
	}else{ord}
}	
