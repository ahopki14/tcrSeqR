#' overlap
#' 
#' Compute the overlap between two vectors. This is a metric of
#' similarity between two populations. 
#'
#' @param x,y Count vectors (integers) to be compared
#' @return The overlap between the vectors
#' @author Alexander Hopkins
#' @export
# Overlap
# A basic overlap function
overlap <- function(x,y){
    shared_x <- x[x>0 & y >0]
    shared_y <- y[x>0 & y >0]
    shared_sum <- sum(c(shared_x,shared_y))
    total <- sum(c(x,y))
    ol <- shared_sum/(total) # Adaptive has a + 1 in the denom...
#
    ol
}
