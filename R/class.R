setClass("tcr",contains='SummarizedExperiment')

tcr <- function(se=SummarizedExperiment(assays=list(tcr_null=matrix()))){
	new('tcr',se)
}

data <- new("tcr", seq=assay(ds),
	    receptors=as.character(rowData(ds)$aa),
	    colData=dict, rowData=as.data.frame(rowData(ds)))

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

#this makes filter work with DataFrame
filter.DataFrame <- function(.data,...) {
	    DataFrame(filter(as.data.frame(.data),...))
}

# filter for tcr object
filter.tcr <- function(.data, ...){
	colData.filter <- filter(colData(.data), ...)
	if(nrow(colData.filter)>0){
		w <- grep(paste0(as.character(colData.filter$fn), collapse='|'), colData(.data)$fn)
		.data[,w]	
	}else{.data[,0]}
}
