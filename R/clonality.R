# clonality
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
					merge(colData(x),df,by='fn', all.x=T)
					)
		# This is necessary to maintain colnames(assay(ds))
		# for some reason...
		rownames(tmp) <- tmp$fn
		colData(x) <- tmp
		return(x)
		}
		else{cl}
	  }
  )
