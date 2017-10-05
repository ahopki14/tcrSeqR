#' iseqr_morisita
#' 
#' Computes the Morisita index (similarity) across multiple comparisons in a tcr
#' object 
#'
#' @param ds A tcr object
#' @param comps A list of comparisons to be made
#' @param merge Logical indicating if the results should be returned alone, or
#' merged into the tcr object
#'
#' @return A data.frame containing the Morisita indicies (if merge=F) or a tcr
#' object with the Morisita indicies included in the metadata
#' @author Alexander Hopkins
#' @export
iseqr_morisita <- function(ds,comps, merge=T){
	#make sure comps is a list with 2 elements in each
	stopifnot(all(unlist(lapply(comps, FUN=length))==2))
	#make sure none of the comparators are duplicated 
	stopifnot(!any(duplicated(sapply(comps, function(x){x[2]}))))
	#pull out the necessaries from the tcr object
	plot_ds <- colData(ds)
	patients <- unique(as.character(plot_ds$patient))
	out <- data.frame(patient=character(), type=character(), Morisita=numeric())
	for(a in seq(length(patients))){
		for(b in comps){
			m <- numeric()
			if(length(which(plot_ds$patient==patients[a] & (plot_ds$type==b[1]|plot_ds$type==b[2]))) == 2){
				m <-morisita(assay(ds)[,which(plot_ds$patient==patients[a] & plot_ds$type==b[[1]])],
					     assay(ds)[,which(plot_ds$patient==patients[a] & plot_ds$type==b[[2]])])
			}
			if(length(m)>0){
				out <- rbind(out, data.frame(patient=patients[a],type=b[[2]],Morisita=m))
			}
			
		}
	}
	if(merge){
		tmp <- merge(colData(ds), out, by=c('patient','type'), all.x=T, sort=F)
		o <-order(tmp$ord)
		tmp <- tmp[o,]
		#make sure the metadata and data are in the same order
		stopifnot(all(tmp$fn == colnames(assay(ds))))
		rownames(tmp) <- tmp$fn
		colData(ds) <- DataFrame(tmp)
		ds
	}else{out}
}
