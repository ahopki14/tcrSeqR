#' tcr class objects and basic methods
#' 
#' This is the constructor function for a tcr class object, which extends
#' SummarizedExperiment
#'
#' @return A tcr class object
#' @export tcr
#' @exportClass tcr
tcr <- setClass("tcr",contains='SummarizedExperiment',
	 prototype=SummarizedExperiment(assays=list(tcr_null=matrix()))
	 )
setMethod("show", "tcr",
	function(object){
		cat(class(object), "Experiment with", nrow(object@assays), 
		    "receptors and", ncol(object@assays), "samples\n")
		if(names(object@assays)=='tcr_nt'){
			cat('Nucleotide Data (NOT Aggregated)\n')}else 
				if(names(object@assays)=='tcr'){
				cat('Aggregated Data\n')}else{
					cat('Empty Data Set\n')
		}
		#
		cat('Metadata Available (', ncol(colData(object)),'):\n', sep='')
		if(ncol(colData(object))>15){
			cat(head(names(colData(object)),3), '...', tail(names(colData(object)),3), '\n')
		}else{
			cat(colnames(colData(object)), '\n')
		}
	}
	# to do:
)
