#' delta_stats
#' 
#' Computes the change in metrics and (optionally) merges the results back into
#' a tcr object
#'
#' @param ds A tcr object
#' @param comps A list of comparison types
#' @param metric A character specifying the name of the field to be used for the
#' change calculation
#' @param merge If x is a tcr object, the output will be a tcr object with the
#' clonality merged into the metadata. 
#' @return If merge=T, a tcr object with 'Log2.Fold.Change.in[metric]' included
#' in the metadata. If merge=F, a data.frame containing the fold changes.
#' @author Alexander Hopkins
#' @export
delta_stats <- function(ds,comps, metric, merge=T){
	md <- as.data.frame(colData(ds))
	patients <- unique(md$patient)
	name <- paste0('Log2.Fold.Change.in.',metric) 
	out <- data.frame(patient=character(), type=character(), name=numeric())
	for(a in seq(length(patients))){
		for(b in comps){
			fc <- filter(md, type==b[[2]] & patient==patients[a])[,metric] /
				filter(md, type==b[[1]] & patient==patients[a])[,metric]
			lfc <- log2(fc)
			if(length(lfc)>0){
				out <- rbind(out, data.frame(patient=patients[a],type=b[[2]],name=lfc))
		}
		}
	}
	names(out)[3] <- name
	if(merge){
		tmp <- merge(md, out, by=c('patient','type'), all=T)
		o <-order(tmp$ord)
		tmp <- tmp[o,]
		#make sure the metadata and data are in the same order
		stopifnot(all(tmp$fn == colnames(assay(ds))))
		tmp <- DataFrame(tmp)
		rownames(tmp) <- tmp$fn
		colData(ds) <- tmp
		ds
	}else{out}
}
