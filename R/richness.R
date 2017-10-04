#richness
richness <- function(x,...){
	length(which(x>0))
}
setMethod("richness", "tcr",
	  function(x, merge=T){
		r <- apply(assay(x),FUN=richness, MARGIN=2)
	  	df <- data.frame(Richness=r,fn=names(r))
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
		else{r}
	  }
  )
