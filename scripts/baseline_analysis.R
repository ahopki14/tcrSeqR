load('/home/ahopkins/Documents/emj/ImmunoseqResults/sampleExport.2015-05-21_09-55-12/ds_agg.Rda')
load('/home/ahopkins/Documents/emj/ImmunoseqResults/sampleExport.2015-05-21_09-55-12/dict.Rda')
cl <- sapply(ds_agg[,-c(1,2)],clonality)
cl <- as.vector(cl)

#clonality plot
day <- dict$day[-c(1,2)]
patient <- as.factor(dict$patient[-c(1,2)])
pdf('/home/ahopkins/Documents/emj/ImmunoseqResults/R/baseline_plots/clonality.pdf',height=6,width=10)
trellis.par.set(strip.background=list(col="lightgrey"))
xyplot(cl ~ day | patient, 
      layout = c(3,1),
      xlab='Date',
      ylab='Clonality',
      scales=list(alternating=FALSE))
dev.off()


# 8x8 overlap heatmap
o <- matrix(data=rep(NA,64),nrow=8)
for(a in 1:3){
  o <- matrix(data=rep(NA,64),nrow=8)
  w <- which(dict$patient==a)
  ds <- ds_agg[,w]
  for(x in seq(8)){
    for(y in seq(8)[-x]){
      o[x,y] <- overlap(ds[ ,x],ds[ ,y])
      o[y,x] <- o[x,y] 
    }
  }  
pdf(paste0('/home/ahopkins/Documents/emj/ImmunoseqResults/R/baseline_plots/overlap_patient_',a,'.pdf'),height=6,width=10)
print(
levelplot(o,
	regions=TRUE,
	#at=seq(0.25,0.65,length=100),
	col.regions = gray(100:0/100),
	xlab="",
	ylab="",
	main=paste0('Overlap (Donor ',a,')')
	)
)		
dev.off()
}

		
		
