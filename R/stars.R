#' stars
#' 
#' Outputs stars to denote p-values in graphics outputs. 
#'
#' @param t A htest (or similar) object 
#' @return A character string containing stars representing the significance of
#' the p-values
#' @author Alexander Hopkins
stars <- function(t){
	p <- t$p.value
	stopifnot(class(p) == 'numeric')
	if(p<=0.001){
		'***'
	}else if(p<=0.01){
		'**'
	}else if(p<=0.05){
		'*'
	} else if(p>0.05){
		''
	}
}
