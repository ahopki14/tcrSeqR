library(scales)
load('ds_agg.Rda')
n <- 30
tds <- ds_agg[,c(2,3,5)]
tds <- tds[tds[,2]>0 | tds[,3]>0,]
rownames(tds) <- seq_len(nrow(tds))
tds$aa <- as.factor(as.character(tds$aa))
o <- overlap(tds[,2],tds[,3]) 
Pre <- rank(-tds[,2],ties.method='min') # this is not quite right...
Post <- rank(-tds[,3],ties.method='min')
#plot(Pre,Post,log='xy',col=alpha('blue',0.2),pch=19,cex=0.75)
pre_b <- which(Pre %in% 1:n)
post_b <- which(Post %in% 1:n)
#plot(Pre[b],Post[b])
desc <- c(rep("Pre",n),rep("Post",n))
desc <- as.factor(desc)
desc <- factor(desc, levels=rev(levels(desc)))
dat <- c(Pre[pre_b],Post[pre_b])
stripchart(
	dat ~ desc,
	ylim=c(5*n,1),
	vertical=TRUE,
	at=c(1.2,1.7),
	pch=19,
	col=alpha('red',0.5),
	ylab='Rank',
	main=paste0('Overlap: ',round(o,2))
	)
dat <- c(Pre[post_b],Post[post_b])
stripchart(
	dat ~ desc,
	add=TRUE,
	vertical=TRUE,
	at=c(1.3,1.8),
	pch=19,
	col=alpha('blue',0.5)
	)
for(a in 1:n){
	segments(1.2,Pre[pre_b][a],1.7,Post[pre_b][a],col=alpha('red',0.5),lend=0)
	segments(1.3,Pre[post_b][a],1.8,Post[post_b][a],col=alpha('blue',0.5),lend=0)
}
