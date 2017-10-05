#' size
#' 
#' A simple function to print the size of an R object, formatted in MB
#' @param x An R object
#' @return The size of x in MB
#' @author Alexander Hopkins
#' @export
size <- function(x){format(object.size(x), units='MB')}
