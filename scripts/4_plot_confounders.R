# use calculate_summary script to generate stats
# load stats and dictionary

#set path for file output
path <- '/home/ahopkins/Documents/emj/ImmunoseqResults/sampleExport.2014-07-31_10-10-24/rerun/plots/'
dir.create(path)

pdf(paste0(path,'number-vs-richness.pdf'),height=6,width=6, title='Number vs Richness')
plot(dict[-c(1,2),'cells'],stats$Richness,
  xlab='Input Cell Number',
  ylab='Richness')
dev.off()

pdf(paste0(path, 'number-vs-clonality.pdf'),height=6,width=6, title='Number vs Clonality')
plot(dict[-c(1,2),'cells'],stats$Clonality,
  xlab='Input Cell Number',
  ylab='Clonality')
dev.off()

pdf(paste0(path, 'age-vs-richness.pdf'),height=6,width=6, title='Age vs Richness')
plot(dict[-c(1,2),'age'],stats$Richness,
  xlab='Age at Diagnosis',
  ylab='Richness')
dev.off()

pdf(paste0(path, 'age-vs-clonality.pdf'),height=6,width=6, title='Age vs Clonality')
plot(dict[-c(1,2),'age'],stats$Clonality,
  xlab='Age at Diagnosis',
  ylab='Clonality')
dev.off()

pdf(paste0(path, 'age-vs-overlap.pdf'),height=6,width=6, title='Age vs Overlap')
plot(dict[-c(1,2),'age'],stats$overlap,
  xlab='Age at Diagnosis',
  ylab='Overlap')
dev.off()

pdf(paste0(path , 'number-vs-overlap.pdf'),height=6,width=6, title='Number vs Overlap')
plot(dict[-c(1,2),'cells'],stats$overlap,
  xlab='Input Cell Number',
  ylab='Overlap')
dev.off()

pdf(paste0(path, 'clonality-vs-overlap.pdf'),height=6,width=6, title='Clonality vs Overlap')
plot(stats$Clonality,stats$overlap,
  xlab='Clonality',
  ylab='Overlap')
dev.off()

pdf(paste0(path, 'richness-vs-overlap.pdf'),height=6,width=6, title='Richness vs Overlap')
plot(stats$Richness,stats$overlap,
  xlab='Richness',
  ylab='Overlap')
dev.off()

##################################Confounders######################################

t <- t.test(dict$age[dict$response=='NR' & dict$type=='PRE'],dict$age[dict$response=='R' & dict$type=='PRE'])

pdf(paste0(path, 'age-by-response.pdf'),height=8,width=4.5, title='Age by Response')
stripchart(dict$age[-c(1,2)] ~ factor(dict$response[-c(1,2)]),
	pch=19,
	main=paste('p=',round(t$p.value,4)),
	vertical=TRUE,
	at=c(1.25,1.75),
	xlim=c(1,2),
	ylab='Age',
	par(mai=c(1,1,1,1))
	)
dev.off()

t <- t.test(dict$cells[dict$response=='NR' & (dict$type=='PRE')],dict$cells[dict$response=='R' & (dict$type=='PRE')])
pdf(paste0(path, 'number-by-response.pdf'),height=8,width=4.5, title='Number by Respone')
stripchart(dict$cells[-c(1,2)] ~ factor(dict$response[-c(1,2)]),
	pch=19,
	main=paste('p=',round(t$p.value,4)),
	vertical=TRUE,
	at=c(1.25,1.75),
	xlim=c(1,2),
	ylab='Input Cell Number',
	par(mai=c(1,1,1,1))
	)
dev.off()
































