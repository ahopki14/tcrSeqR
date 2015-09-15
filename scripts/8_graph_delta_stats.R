# load dictionary and stats (see calculate_summary)

path <- '/home/ahopkins/Documents/emj/ImmunoseqResults/sampleExport.2014-07-31_10-10-24/rerun/plots/'

if(!is.na(stats[1,1])){
  stats <- rbind(rep(NA,ncol(stats)),rep(NA,ncol(stats)),stats) # pad the stats with NA so dict works
}
patients <- unique(dict$patient[-c(1,2)])
out <- matrix(NA,nrow=length(patients),ncol=ncol(stats)+2)
out <- as.data.frame(out)
names(out) <- c('patient','response',names(stats))

for(a in seq_along(patients)){
  w <- which(dict$patient==patients[a] & dict$type!='PDAC')
  tstats <- stats[w, ]
  tdict <- dict[w, ]
  for(b in seq_along(tstats)){
    # calculate p-values
     pct_change <- 100*((tstats[tdict$type=='POST' ,b] - tstats[tdict$type=='PRE',b]) /
                        tstats[tdict$type=='PRE',b])
     out$patient[a] <- patients[a]
     out$response[a] <- tdict$response[tdict$patient==patients[a] & tdict$type=='PRE']
     dir.create(paste0(path, colnames(tstats)[b],'/'), recursive=TRUE)
     pdf(paste0(path, colnames(tstats)[b],'/',types[a],'.pdf'),width=4.5,height=8)
       par(oma=c(1,2,1,1))
       stripchart(tstats[ ,b] ~ tdict$response,
                  at=c(1.25,1.75),
                  pch=19,
                  vertical=TRUE,
                  xlim=c(1,2),
                  ylab=colnames(tstats)[b],
                  main=paste0(types[a])
                  )
       text(1.5,max(c(r,nr)),paste0("p=",as.character(round(p,4))))
       dev.off()
  }
}

