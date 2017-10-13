#' filter.DataFrame
#' 
#' Method for using dplyr filter on DataFrames
#' @author Alexander Hopkins
#' @export
#' @method filter DataFrame
filter = function(x,...)
{
	    UseMethod("filter")
}
filter.DataFrame <- function(.data,...) {
	    DataFrame(filter(as.data.frame(.data),...))
}
