exp_clone <- function(x,y){ # x and y are vectors of counts in the samples compared
    mat <- as.matrix(cbind(x,y))
    s <- apply(mat,MARGIN=2,FUN=sum)
    p_vals <- apply(mat,MARGIN=1,FUN=fisher,s=s)
    p_vals
}

