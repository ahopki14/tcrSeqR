#' mean_sd
#' 
#' Computes the mean +/- standard deviation. For use with stat_summary function
#' in ggplot 
#'
#' @param x data to be summarized
#' @author Alexander Hopkins
mean_sd <- function(x){
	m <- mean(x,na.rm=TRUE)
	s <- sd(x,na.rm=TRUE)
	c(y=m,ymin=m-s,ymax=m+s)
}

#' mean_95
#' 
#' Computes the mean and 95% confidence intercal. For use with stat_summary function
#' in ggplot 
#'
#' @param x data to be summarized
#' @author Alexander Hopkins
mean_95 <- function(x){
	m <- mean(x,na.rm=TRUE)
	n <- length(x)
	i <- qt(0.975,df=n-1)*(sd(x)/sqrt(n))
	c(y=m,ymin=m-i,ymax=m+i)
}

#' mean_only
#' 
#' Computes the mean and outputs it for use with stat_summary function
#' in ggplot 
#'
#' @param x data to be summarized
#' @author Alexander Hopkins
mean_only <- function(x){
	m <- mean(x,na.rm=TRUE)
	c(y=m,ymin=m,ymax=m)
}
