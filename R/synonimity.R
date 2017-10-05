#' synonymity
#' 
#' Computes the synonymity (how many times each unique element appears). In
#' hindsight, this does exactly what 'table' does.
#'
#' @param x A character list with (probably) repeats
#' @return A named list, the length of unique(x) counting each occurrence
#' @author Alexander Hopkins
#' @export
synonymity <- function(x){
    x <- as.character(x)
    ind <- split(seq_len(length(x)),x)
    syn <- sapply(ind,length)
#
    syn
}
