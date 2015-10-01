# load dictionary and stats (see calculate_summary)

#path <- '/home/ahopkins/Documents/emj/ImmunoseqResults/adjuvant_study/plots/'

types <- levels(dict$type)
if(!is.na(stats[1,1])){
  stats <- rbind(rep(NA,ncol(stats)),rep(NA,ncol(stats)),stats) # pad the stats with NA so dict works
}
for(a in seq_along(types)){
  w <- which(dict$type==types[a])
  tstats <- stats[w, ]
  tdict <- dict[w, ]
  for(b in seq_along(tstats)){
    # calculate p-values
    r <- tstats[tdict$response=='R', b]
    nr <- tstats[tdict$response=='NR', b]
    p <- t.test(r,nr,var.equal=TRUE)$p.value
    dir.create(paste0(path, colnames(tstats)[b],'/'), recursive=TRUE)
    pdf(paste0(path, colnames(tstats)[b],'/',types[a],'.pdf'),
                   width=4.5,height=8,
                   title=paste0(colnames(tstats)[b],' ',types[a]))
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

## overlap
dir.create(paste0(path, 'Pre-Post/'), recursive=TRUE)
resp <- dict$response[which(dict$type=='PRE')]
r <- ol[resp=='R']
nr <- ol[resp=='NR']
p <- t.test(r,nr,var.equal=TRUE)$p.value
pdf(paste0(path,'Pre-Post/overlap.pdf'), width=4.5,height=8,title="Pre-Post Overlap")
stripchart(ol ~ resp,
          at=c(1.25,1.75),
          pch=19,
          vertical=TRUE,
          xlim=c(1,2),
          ylab='Pre-Post Overlap',
          )
text(1.5,max(c(r,nr)),paste0("p=",as.character(round(p,4))))
dev.off()

#Morisita
resp <- dict$response[which(dict$type=='PRE')]
r <- m[resp=='R']
nr <- m[resp=='NR']
p <- t.test(r,nr,var.equal=TRUE)$p.value
pdf(paste0(path,'Pre-Post/morisita_overlap.pdf'), width=4.5,height=8,title="Pre-Post Morisita Overlap")
stripchart(m ~ resp,
          at=c(1.25,1.75),
          pch=19,
          vertical=TRUE,
          xlim=c(1,2),
          ylab='Pre-Post Morisita Overlap',
          )
text(1.5,max(c(r,nr)),paste0("p=",as.character(round(p,4))))
dev.off()


