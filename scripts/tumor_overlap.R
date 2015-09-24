out <- '~/Desktop/out-check.txt'
for(a in levels(as.factor(dict$patient))){
	w <- which(dict$patient==a & (dict$type=='PRE' | dict$type=='PDAC'))
	tds <- ds[ ,w]
	count <- nrow(tds[tds[ ,1]>0 & tds[ ,2]>0, ])
	total <- nrow(tds[tds[ ,2]>0 ,])
	pct <- round(100*(count/total),2)
	cat(paste0(as.character(a),
		   ': ', count, ' of ',total,' tumor clones...',pct,' percent\n' ),
		   file=out,append=TRUE)
}

