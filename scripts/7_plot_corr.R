# load the data and the dictionary
patients <- levels(as.factor(dict$patient))
#path <- '/home/ahopkins/Documents/emj/ImmunoseqResults/adjuvant_study/plots/'

for(a in seq_along(patients)){
  # make data subsets
  w <- which(dict$patient==patients[a])
  response <- as.character(dict$response[which(dict$patient==patients[a] & dict$type=='PRE')])
  tds <- ds[ ,w]
  tdict <- dict[w, ]
  # get pre and post
  pre <- tds[ ,which(tdict$type=='PRE')]
  post <- tds[ ,which(tdict$type=='POST')]
  x <- (pre!=0 | post!=0)
  pre <- pre[x]
  post <- post[x]
  # calculate overlap
  ol <- overlap(pre,post)
  # plot
  dir.create(paste0(path,'Corr/'))
  pdf(paste0(path,'Corr/',patients[a],'_corr.pdf'),
      width=8.5,height=8,
      title=paste('Correlation Plot for',patients[a]))
  par(oma=c(1,4,1,1))
  r <- range(pre,post)[2]*1.05	
	plot(post,pre,pch=".",asp=1,
	  		main=paste0(patients[a],'(',response,')'),
				xlim=c(1,r),ylim=c(1,r),
				xlab='Post',
				ylab='Pre',
				log='xy',
				cex.lab=1.5,
				cex.main=2)
	grid()
	text(r/100,r,paste0("o=",round(ol,4)),cex=1.5)
	dev.off()
}
