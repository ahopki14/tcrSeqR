#' filter.DataFrame
#' 
#' Method for using dplyr filter on DataFrames
#' @author Alexander Hopkins
filter.DataFrame <- function(.data,...) {
	    DataFrame(filter(as.data.frame(.data),...))
}
