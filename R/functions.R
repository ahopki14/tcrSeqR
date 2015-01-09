ss2n <- function(x,y){
	mean_x <- apply(x,1,mean)
	sd_x <- apply(x,1,sd)
	mean_y <- apply(y,1,mean)
	sd_y <- apply(y,1,sd)
	s2n <- (mean_x - mean_y) / (sd_x + sd_y)
	s2n[is.na(s2n)==TRUE] <- 0 # enforce sd > 0 
	s2n[is.infinite(s2n)==TRUE] <- 0
	s2n
}
