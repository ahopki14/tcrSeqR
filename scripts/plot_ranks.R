library(scales)
library(immunoSeqR)
load('~/Documents/emj/ImmunoseqResults/sampleExport.2014-07-31_10-10-24/ds_agg.Rda')
n <- 100
tds <- ds_agg[,c(2,4,6,3)]
name <- '7.006'
tds <- tds[tds[,2]>0 | tds[,3]>0,]
tds <- tds[tds[,4]>0,] # restrict to tumor
rownames(tds) <- seq_len(nrow(tds))
tds$aa <- as.factor(as.character(tds$aa))
o <- overlap(tds[,2],tds[,3]) 
Pre <- rank(-tds[,2],ties.method='min') 
Post <- rank(-tds[,3],ties.method='min')
#plot(Pre,Post,log='xy',col=alpha('blue',0.2),pch=19,cex=0.75)
pre_b <- which(Pre %in% 1:n)
post_b <- which(Post %in% 1:n)
desc <- c(rep("Pre",length(pre_b)),rep("Post",length(post_b)))
desc <- as.factor(desc)
desc <- factor(desc, levels=rev(levels(desc)))
dat <- c(Pre[pre_b],Post[pre_b])
pdf(paste0('~/Documents/emj/ImmunoseqResults/R/plots/Rankplots/',name,'.pdf'),width=8.5,height=8)
stripchart( #makes an empty chart to draw lines on
	c(0,0)~factor(c('Pre','Post'),levels=c('Pre','Post')),
	at=c(1.25,1.75),
	ylim=c(3*n,1),
	xlim=c(1.2,1.8),
	vertical=TRUE,
	col='white',
	ylab='Rank',
	frame.plot=FALSE,
	main=name
	)
for(a in 1:n){
	segments(1.2,Pre[pre_b][a],1.7,Post[pre_b][a],col=alpha('red',0.4),lend=0)
	segments(1.3,Pre[post_b][a],1.8,Post[post_b][a],col=alpha('blue',0.4),lend=0)
}
dev.off()

