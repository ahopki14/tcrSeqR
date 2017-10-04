tert <- function(x){
	q <- quantile(x,c(1/3,2/3))
	w1 <- which(x<=q[[1]])
	w2 <- which(x<=q[[2]] & x>q[[1]])
	w3 <- which(x>q[[2]])
	out <- character(length=length(x))
	out[w1] <- 1
	out[w2] <- 2
	out[w3] <- 3
	out <- factor(out,levels=c('1','2','3'))
	out
}
