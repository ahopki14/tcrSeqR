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
	# make indicies and synonymity
	ind <- split(seq_len(nrow(ds)),rownames(ds)) # This is a little slow
	syn <- sapply(ind,length)
	#split out columns to be fixed
	f <- syn>1 # id clones which need to be collapsed
	to_fix <- as.vector(unlist(ind[f]))
	ok <- ds[-to_fix,]
	stopifnot(length(rownames(ok))==length(unique(rownames(ok))))
	to_fix <- ds[to_fix,]
	cat('Found', nrow(to_fix),'synonymous sequences\n')
	# do math
	gc()
	ind_to_fix <- split(seq_len(nrow(to_fix)),rownames(to_fix))
	to_fix_m <- assay(to_fix)
	sum_rows <- function(x){colSums(to_fix_m[unlist(ind_to_fix[x]),])}
	fixed_mat <- sapply(seq(length(ind_to_fix)),FUN=sum_rows) # This is the slow part
	fixed_mat <- t(fixed_mat)
	rownames(fixed_mat) <- names(ind_to_fix)
	#assemble a tcr of the fixed sequences
	fixed <- tcr(SummarizedExperiment(assays=list(tcr_nt=fixed_mat)))
	rowData(fixed) <- DataFrame(aa=rownames(fixed_mat),nt='NA')
	ds_out <- BiocGenerics::rbind(ok, fixed)
	names(ds_out@assays$data) <- 'tcr'
	# Print sanity check
	t <- round((proc.time()-start)[['elapsed']]/60,2)
	cat(paste0('Collapsed ',sum(syn[syn>1]),' TCRs to ',length(syn[syn>1]),'\n',
		   'Left ',length(syn[syn==1]),' TCRs\n'))
	stopifnot(nrow(ds_out) + sum(syn[syn>1]) - length(syn[syn>1]) == nrow(ds))
	cat(paste0('Completed in ',t,' minutes\n'))
	#
	ds_out
}
