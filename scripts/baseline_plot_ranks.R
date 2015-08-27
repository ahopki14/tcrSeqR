library(scales)
#library(immunoSeqR)
load('/home/ahopkins/Documents/emj/ImmunoseqResults/sampleExport.2015-05-21_09-55-12.ds_agg.Rda')
load('/home/ahopkins/Documents/emj/ImmunoseqResults/sampleExport.2015-05-21_09-55-12/dict.Rda')

n <- 50
#ds <- ds_agg
#dict <- dict

pdf(paste0('~/Documents/emj/ImmunoseqResults/R/baseline_plots/Rankplots/all-',n,'-all_tp.pdf'), width=10,height=40)
par(mfrow=c(7,3))
par(mar=c(3,3,2,3))
for(tp in 1:7){
for(a in 1:3){
	#resp <- as.character(dict$response[which(dict$patient==name)][1])
	#samples <- c(2, # include the aa data
	#	which(dict$patient==name & dict$type=='Pre'),
	#	which(dict$patient==name & dict$type=='Post'),
	#	which(dict$patient==name & dict$type=='PDAC')
	#)
	#names(ds)[samples]
	days <- sort(dict$day[which(dict$patient==a)])[tp:(tp+1)] # gets the earliest 2 time points
	samples <- c(2,
	  which(dict$patient==a & dict$day==days[1]),
	  which(dict$patient==a & dict$day==days[2])
	) 
	tds <- ds[,samples]
	tds <- tds[tds[,2]>0 | tds[,3]>0,]
	rownames(tds) <- seq_len(nrow(tds))
	tds$aa <- as.factor(as.character(tds$aa))
	#o <- overlap(tds[,2],tds[,3]) 
	Pre <- rank(-tds[,2],ties.method='min') 
	Post <- rank(-tds[,3],ties.method='min')
	#plot(Pre,Post,log='xy',col=alpha('blue',0.2),pch=19,cex=0.75)
	pre_b <- which(Pre %in% 1:n)
	post_b <- which(Post %in% 1:n)
	#desc <- c(rep("Pre",length(pre_b)),rep("Post",length(post_b)))
	#desc <- as.factor(desc)
	#desc <- factor(desc, levels=rev(levels(desc)))
	dat <- c(Pre[pre_b],Post[pre_b])

	stripchart( #makes an empty chart to draw lines on
		c(0,0)~factor(days,levels=days),
		at=c(1.25,1.75),
		ylim=c((3*n),1),
		xlim=c(1.2,1.8),
		vertical=TRUE,
		col='white',
		ylab='Rank',
		frame.plot=FALSE,
		main=paste0('Donor ',a)#,'(',resp,')')#, '\n Present in Tumor')
		)
	for(i in 1:n){
		segments(1.2,Pre[pre_b][i],1.7,Post[pre_b][i],col=alpha('red',0.4),lend=0,lwd=0.1)
		segments(1.3,Pre[post_b][i],1.8,Post[post_b][i],col=alpha('blue',0.4),lend=0,lwd=0.1)
	}
}
}
dev.off()


