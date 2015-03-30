library(immunoSeqR)
load('~/Documents/emj/ImmunoseqResults/sampleExport.2014-07-31_10-10-24/ds_agg_reorder.Rda')
load('~/Documents/emj/ImmunoseqResults/sampleExport.2014-07-31_10-10-24/dict_reorder.Rda')
ds <- ds_agg_reorder
dict <- dict_reorder

scans <- list()
for(a in seq(length(levels(dict$patient)))){
	name <- levels(dict$patient)[a]
	resp <- as.character(dict$response[which(dict$patient==name)][1])
	samples <- c(which(dict$patient==name & dict$type=='Pre'),
		           which(dict$patient==name & dict$type=='Post')
	            )
	byN=100
  tmp <- olScan(ds[,samples[1]],ds[,samples[2]],byN=byN)
  scans[[name]] <- tmp
  #save(tmp,file=paste0('~/Documents/emj/ImmunoseqResults/olscan/',name,'_by_',byN,'.Rda'))
}
save(scans,file=paste0('~/Documents/emj/ImmunoseqResults/olscan/all_by_',byN,'.Rda'))
#make plotting stuff
biggest <- max(unlist(lapply(scans,length)))

pdf(file='~/Documents/emj/ImmunoseqResults/olscan/olscan.pdf',width=6,height=6)
plot(NA, 
  ylim=c(0,1),
  xlim=c(1,biggest),
  xlab='# of clones',ylab='Overlap',
  main='Overlap of top X clones')


for(a in seq(length(scans))){
resp <- as.character(dict$response[which(dict$patient==names(scans)[a])][1])
if(resp=='R'){color<-'blue'}else{color='red'}
points(seq(length(scans[[a]])),scans[[a]],col=color,pch='.')
}  
legend('topright',legend=c('R','NR'),col=c('blue','red'),pch=19)
dev.off()
