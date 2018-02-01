#' richness
#' 
#' Computes the richness (total number of unique elements) in a sample
#' @param x A numeric vector or tcr object
#' @param merge If x is a tcr object, the output will be a tcr object with the
#' richness merged into the metadata. 
#' @return If x is a vector, the function returns the richness If x is a tcr
#' object, and merge=T, then the function returns a tcr object with richness
#' included in the metadata for all samples.  
#' @author Alexander Hopkins
#' @export
richness <- function(x,...){
	length(which(x>0))
}
setMethod("richness", "tcr",
	  function(x, merge=T){
		r <- apply(assay(x),FUN=richness, MARGIN=2)
	  	df <- data.frame(Richness=r,fn=names(r))
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
		else{r}
	  }
  )
