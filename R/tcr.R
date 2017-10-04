setClass("tcr",contains='SummarizedExperiment')

tcr <- function(se=SummarizedExperiment(assays=list(tcr_null=matrix()))){
	new('tcr',se)
}


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
		cat('Metadata Available:\n')
		cat(head(names(colData(object)),15), '\n')
	}
	# to do:
)
