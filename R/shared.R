#' shared
#' 
#' Computes the number of clones that are present in two samples.
#' 
#'
#' @param x,y Numeric vectors of the counts of clones in each sample
#' @return The number of clones that are present in both samples
#' @author Alexander Hopkins
#' @export
shared <- function(x,y){
	if(length(x)!=length(y)){stop('Vectors must be same length')}
	w <- which(x>0)
	out <- length(which(y[w]>0))
	out
}
