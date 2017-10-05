#' exp_clone
#' 
#' Calculates p values for each clone using a fisher test. 
#'
#' @param x,y Vectors representing the receptor counts from two samples
#'
#' @return A vector of fisher test p values for each clone, the same length as x and y
#'
#' @author Alexander Hopkins
exp_clone <- function(x,y){ # x and y are vectors of counts in the samples compared
    mat <- as.matrix(cbind(x,y))
    s <- apply(mat,MARGIN=2,FUN=sum)
    p_vals <- apply(mat,MARGIN=1,FUN=fisher,s=s)
    p_vals
}

