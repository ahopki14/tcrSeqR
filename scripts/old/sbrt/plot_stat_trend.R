
if(!is.na(stats[1,1])){
  stats <- rbind(rep(NA,ncol(stats)),rep(NA,ncol(stats)),stats) # pad the stats with NA so dict works
}
patients <- levels(dict$patient)
for(b in seq(ncol(stats))){
	pdf(paste0(path,'/',names(stats)[b],'/',names(stats)[b],'-tracking.pdf'),width=4.5,height=8)
	plot(
		c(0.5,2,3.25),
		c(min(stats[-c(1,2),b]),mean(stats[-c(1,2),b]),max(stats[-c(1,2),b])),
		type='n',ylab=names(stats)[b],xlab='',xaxt='n'
	    )
	axis(1,at=c(1,2,3),labels=c('Pre','Post SBRT','Post FOLF'))
	
	for(a in seq_along(patients)){
		w <- c(which(dict$patient==patients[a] & (dict$type=='PRETREAT')),
			which(dict$patient==patients[a] & (dict$type=='POSTSBRT')),
			which(dict$patient==patients[a] & (dict$type=='POSTFOLF'))
		)
	arm <- dict$arm[w][1]
	if(arm == 1){col <- 'red'} else {col <- 'black'}
	if(length(w)==3){x <- seq(3)} else {x <- c(1,2)}
	y <- stats[w,b]/stats[w,b][1]
	lines(x,stats[w,b],col=col)
	text(1,stats[w[1],b],as.character(patients[a]),cex=0.4,col='grey',pos=2)
	}
dev.off()
}

for(b in seq(ncol(stats))){
        pdf(paste0(path,'/',names(stats)[b],'/',names(stats)[b],'-tracking-norm.pdf'),width=4.5,height=8)
        plot(
                c(0.5,2,3.25),c(0,0,0),
                ylim=c(-3,3),
                type='n',
                ylab=paste0('log2(Normalized ',names(stats)[b],')'),
                xlab='',
                xaxt='n'
            )
        axis(1,at=c(1,2,3),labels=c('Pre','Post SBRT','Post FOLF'))

        for(a in seq_along(patients)){
                w <- c(which(dict$patient==patients[a] & (dict$type=='PRETREAT')),
                        which(dict$patient==patients[a] & (dict$type=='POSTSBRT')),
                        which(dict$patient==patients[a] & (dict$type=='POSTFOLF'))
                )
        arm <- dict$arm[w][1]
        if(arm == 1){col <- 'red'} else {col <- 'black'}
        if(length(w)==3){x <- seq(3)} else {x <- c(1,2)}
        y <- log(stats[w,b]/stats[w,b][1],2)
        #y <- 100*(stats[w,b]/stats[w,b][1])
#       lines(x,stats[w,b],col=col)
        lines(x,y,col=col)
#       text(1,stats[w[1],b],as.character(patients[a]),cex=0.4,col='grey',pos=2)
        text(3,tail(y,n=1),as.character(patients[a]),cex=0.4,col='grey',pos=4)
        }
dev.off()
}

