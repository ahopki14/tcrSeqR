#' clonality
#' 
#' Computes the clonality (normalized Shannon Entropy) of vector or all samples
#' in a tcr object
#'
#' @param x A numeric vector or tcr object
#' @param entropy Logical indicating if output should be Entropy (TRUE) or
#' clonality (FALSE), given by 1 - Entropy.
#' @param merge If x is a tcr object, the output will be a tcr object with the
#' clonality merged into the metadata. 
#' @return If x is a vector, the function returns the clonality. If x is a tcr
#' object, and merge=T, then the function returns a tcr object with clonality
#' included in the metadata for all samples.  
#' @author Alexander Hopkins
#' @export

clonality <- function(x,entropy=FALSE, ...){
    x <- x[x != 0]
    x <- x/sum(x)
    cl <- -sum(x*log(x))/log(length(x)) # Entropy
    if(!entropy){
    cl <- 1-cl
    }
    cl
}
setMethod("clonality", "tcr",
	  function(x, entropy=FALSE, merge=T){
		cl <- apply(assay(x),FUN=clonality, MARGIN=2, entropy=entropy)
	  	df <- data.frame(Clonality=cl,fn=names(cl))
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
		else{cl}
	  }
  )
