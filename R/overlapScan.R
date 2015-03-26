olScan <- function(x,y){
	ds <- data.frame(x=x,y=y)
	ds <- ds[order(-ds$x), ]
	ol <- vector()
	for(a in 10:nrow(ds)){
	  ol <- c(ol,overlap(ds$x[1:a],ds$y[1:a]))
	}
	ol
}
