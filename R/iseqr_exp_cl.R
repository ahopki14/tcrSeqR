#' iseqr_exp_cl
#' 
#' Calculates the number of clones that expand relative to a baseline sample in
#' a tcr experiment
#' 
#' @param ds A tcr object
#' @param comps A list of comparison types
#' @param category The field to look for 'comps'. Default: type
#' @param alpha The threshold over which clones with BH adjusted p-values will
#' be considered 'expanded'. Default: 0.05
#' @param min.count The minimum number of clones (sum of both samples) to be
#' included in the analysis, Default: 5
#; @param output One of 'tcr' (default), 'num' or 'p_vals'
#' @return If output=tcr, a tcr object with the number of expanded clones merge
#' into the dataset. If output='nums', a data.frame with the number of expanded
#' clones for each comparison. If output='p_vals', a matrix (ncol(assay(ds)) x
#' length(patients)) with the p-values for each clone in each comparison.
#' @note It is not recommended (but not STRICTLY enforced) that you use the same
#' baseline sample for comparisons (i.e. the first elements of all items in
#' comps should be the same). If you REALLY know what you are doing, then you
#' could use multiple baseline samples. 
#' @note If output='p_vals', then only one comparison can be made at a time
#' (length(comps) should be 1). 
#' @author Alexander Hopkins
#' @export
iseqr_exp_cl <- function(ds,comps, category='type', alpha=0.05,min.count=5,output='tcr'){
	by='patient'
	stopifnot(class(ds)=='tcr')
	#check that we are using nt data
       if(names(ds@assays)!='tcr_nt'){
	       warning('Expanded Clones should be calculated from nucleotide data (NOT aggregated)')
       }	
	#check if all comps are compared to the same baseline samples
	if(length(unique(sapply(comps, function(x){x[[1]]})))>1){
		warning('Comparisons to different baselines NOT recommended')}
	#check that only one comp is given if output='p_vals'
	if(output=='p_vals'){
		if(length(comps)!=1){
			stop('Only one comparison can be made at a time if output=="p_vals"')
		}
	}
	# extract metadata from tcr
	md <- as.data.frame(colData(ds))
	#initialize variables
	patients <- unique(md$patient)
	out <- data.frame(p_adj=numeric(),ind=numeric(),patient=character())
	p <- matrix(nrow=nrow(ds),ncol=length(patients))
	num_exp <- data.frame(patient=character(), type=character(), out=numeric())
	#loop through comparisons first
	for(a in comps){
		#then loop through patients
		for(b in seq_along(patients)){
			w <- c(which(md[,by]==patients[b] & md[,category]==a[[1]]),
			       which(md[,by]==patients[b] & md[,category]==a[[2]])
			       )
			if(length(w)>2){stop(paste('Multiple matches for:',patients[b]))
			}else if(length(w)==2){
				#do math
				  mat <- as.data.frame(assay(ds)[,w])
				  rownames(mat) <- seq(nrow(mat))
				  mat$p_adj <- rep(1,nrow(mat))
				  mat <- as.matrix(mat)
				  ind <- which((mat[,1]+mat[,2])>min.count) # at least 5 in both
				  ind <- as.numeric(ind)
				  s <- apply(mat[ind,],MARGIN=2,FUN=sum)
				  tm <- proc.time()
				  p_vals <- exp_clone(mat[ind,1],mat[ind,2])
				  el<- proc.time()-tm
				  p_adj <- p.adjust(p_vals,method='BH')
				  mat[ind,3] <- p_adj
				  n_exp <- length(which(mat[,3]<alpha))
				  num_exp <- rbind(num_exp, data.frame(patient=patients[b],
								       type=a[[2]],
								       out=n_exp))
				  p[,b] <- mat[,3]
			}else if(length(w)<2){
				#nothing found for that category
				cat('No sample found for',by,patients[b],category,a[[2]],'\n')
			}
		}
	}
	#change the names of the output (in case by!='patient', etc)
	names(num_exp)[1] <- by
	names(num_exp)[2] <- category
	names(num_exp)[3] <- 'Number.of.Expanded.Clones'
	if(output=='tcr'){
		#return the tcr object 
		tmp <- merge(md, num_exp, by=c('patient','type'), all=T)
		o <-order(tmp$ord)
		tmp <- tmp[o,]
		#make sure the metadata and data are in the same order
		stopifnot(all(tmp$fn == colnames(assay(ds))))
		tmp <- DataFrame(tmp)
		rownames(tmp) <- tmp$fn
		colData(ds) <- tmp
		ds

	}else if(output=='p_vals'){
		#or return the p_values 
		p
	}else if(output=='num'){
		#or the number of expanded clones only
		num_exp
	}else{stop("output must be 'tcr', 'num' or 'p_vals'")}
}
