#' iseqr_aggregate
#' 
#' Collapse synonymouse nucleotide sequences, summing the values in each sample
#'
#' @param ds A TCR object made with iseqr_make_tcr.
#' 
#' @return the aggregated dataset
#' @author Alexander Hopkins
#' @export
iseqr_aggregate <- function(ds){
	stopifnot(class(ds)=='tcr')
	stopifnot(assayNames(ds)=='tcr_nt')
	start <- proc.time()
	# restrict to productive and sanitize factor
	w_remove <- which(rownames(ds)=='')
	w_keep <- which(rownames(ds)!='')
	ds <- ds[w_keep,]
	#ds$aa <- as.factor(as.character(ds$aa))
	# make indicies and synonymity
	ind <- split(seq_len(nrow(ds)),rownames(ds))
	syn <- sapply(ind,length)
	#split out columns to be fixed
	f <- syn>1 # id clones which need to be collapsed
	to_fix <- as.vector(unlist(ind[f]))
	ok <- ds[-to_fix,]
	stopifnot(length(rownames(ok))==length(unique(rownames(ok))))
	to_fix <- ds[to_fix,]
	# do math
	gc()
	ind_to_fix <- split(seq_len(nrow(to_fix)),rownames(to_fix))
	sum_rows <- function(x){colSums(assay(to_fix)[unlist(ind_to_fix[x]),])}
	fixed_mat <- sapply(seq(length(ind_to_fix)),FUN=sum_rows)
	fixed_mat <- t(fixed_mat)
	rownames(fixed_mat) <- names(ind_to_fix)
	#assemble a tcr of the fixed sequences
	fixed <- tcr(SummarizedExperiment(assays=list(tcr_nt=fixed_mat)))
	rowData(fixed) <- DataFrame(nt='NA',aa=rownames(fixed_mat))
	ds_out <- BiocGenerics::rbind(ok, fixed)
	names(ds_out@assays$data) <- 'tcr'
	# Print sanity check
	t <- round((proc.time()-start)[['elapsed']]/60,2)
	cat(paste0('Collapsed ',sum(syn[syn>1]),' TCRs to ',length(syn[syn>1]),'\n',
		   'Left ',length(syn[syn==1]),' TCRs\n'))
	cat(paste0('Removed ',length(w_remove),' Non-productive TCRs\n'))
	stopifnot(nrow(ds_out) + sum(syn[syn>1]) - length(syn[syn>1]) == nrow(ds))
	cat(paste0('Completed in ',t,' minutes\n'))
	#
	ds_out
}