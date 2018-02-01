#' Total sequences
#' 
#' Computes the sum of a vector
#' @param x A numeric vector or tcr object
#' @param merge If x is a tcr object, the output will be a tcr object with the
#' richness merged into the metadata. 
#' @return If x is a vector, the function returns the sum. If x is a tcr
#' object, and merge=T, then the function returns a tcr object with
#' Total.Sequences included in the metadata for all samples.  
#' @author Alexander Hopkins
#' @export
total <- function(x,...){
	sum(x)
}
setMethod("total", "tcr",
	  function(x, merge=T){
		tot <- apply(assay(x),FUN=total, MARGIN=2)
	  	df <- data.frame(Total.Sequences=tot,fn=names(tot))
		if(merge){
		tmp <- DataFrame(
					merge(colData(x),df,by='fn', all.x=T, sort=F)
					)
		# This is necessary to maintain colnames(assay(ds))
		# for some reason...
		rownames(tmp) <- tmp$fn
		stopifnot(tmp$fn==colnames(assay(x))) #make sure they line up
		colData(x) <- tmp
		return(x)
		}
		else{tot}
	  }
  )
