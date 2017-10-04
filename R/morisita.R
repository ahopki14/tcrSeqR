#' morisita
#' 
#' Compute the Morisita distance between two vectors. This is a metric of
#' similarity between two populations. 
#'
#' @param x,y Count vectors (integers) to be compared
#' @return The Morisita distance between the vectors
#' @author Alexander Hopkins
#' @export
morisita <- function(x,y){
    stopifnot(length(x) == length(y))
    prod <- x*y
    m <- (2*(sum(prod))) / ((simpson(x) + simpson(y))*sum(x)*sum(y))
    m
}
