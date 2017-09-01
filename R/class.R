setClass("tcr",representation(
			      seq="matrix",
			      receptors="character",
			      colData="data.frame",
			      rowData="data.frame"
			      )
# to do:
)



data <- new("tcr", seq=assay(ds),
	    receptors=as.character(rowData(ds)$aa),
	    colData=dict, rowData=as.data.frame(rowData(ds)))

setMethod("show", "tcr",
	function(object){
		cat(class(object), "Experiment with", nrow(object@seq), 
		    "receptors and", ncol(object@seq), "samples\n")
		if(!any(grepl('syn',names(object@rowData)))){
			cat('Nucleotide Data (NOT Aggregated)\n')}else{
				cat('Aggregated Data\n')
		}
		#
		cat('Metadata:\n')
		print(head(object@colData))
		cat(rep(' ', 15), '.\n')
		cat(rep(' ', 15), '.\n')
		cat(rep(' ', 15), '.\n')
	}
	# to do:
)
