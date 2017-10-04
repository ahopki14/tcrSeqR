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
