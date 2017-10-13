#' filter.tcr
#' 
#' Method for using dplyr filter on tcr class objects
#' @author Alexander Hopkins
#' @export
#' @method dplyr::filter tcr
filter.tcr <- function(.data, ...){
	colData.filter <- filter(colData(.data), ...)
	if(nrow(colData.filter)>0){
		w <- grep(paste0(as.character(colData.filter$fn), collapse='|'), colData(.data)$fn)
		.data[,w]	
	}else{.data[,0]}
}
