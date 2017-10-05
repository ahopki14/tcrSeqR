#' simpson
#' 
#' Computes the Simpson Index (a diversity measure) of a distribution. This is
#' used in the calculation of the Morisita Index.
#'
#' @param x An integer vector indicating the counts
#' @return The Simpson Index
#' @author Alexander Hopkins
simpson <- function(x){
    x <- x/sum(x)
    x <- x^2
    l <- sum(x)
    l
}
