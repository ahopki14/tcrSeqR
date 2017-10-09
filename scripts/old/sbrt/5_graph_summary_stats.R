# load dictionary and stats (see calculate_summary)

#path <- '/home/ahopkins/Documents/emj/ImmunoseqResults/adjuvant_study/plots/'

#per sample stats
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
    r <- tstats[tdict$arm==1, b]
    nr <- tstats[tdict$arm==2, b]
    p <- t.test(r,nr,var.equal=TRUE)$p.value
    dir.create(paste0(path, colnames(tstats)[b],'/'), recursive=TRUE)
    pdf(paste0(path, colnames(tstats)[b],'/',types[a],'.pdf'),
                   width=4.5,height=8,
                   title=paste0(colnames(tstats)[b],' ',types[a]))
    par(oma=c(1,2,1,1))	
    stripchart(tstats[ ,b] ~ tdict$arm,
      at=c(1.25,1.75),
			pch=19,
			vertical=TRUE,
			xlim=c(1,2),
			ylab=colnames(tstats)[b],
			main=paste0(types[a])
		)
		text(1.5,max(c(r,nr)),
		     paste0("p=",as.character(round(p,4))))
		dev.off()
  } 
}


#comparative stats
patients <- levels(as.factor(dict$patient))
#make empty lists
ol_pre_post <- rep(NA,length(levels(dict$patient)))
names(ol_pre_post) <- as.character(levels(dict$patient))
ol_pre_postfolf <- ol_pre_post
ol_postsbrt_postfolf <- ol_pre_post
m_pre_post <- ol_pre_post
m_pre_postfolf <- ol_pre_post
m_postsbrt_postfolf <- ol_pre_post
arm <- ol_pre_post

#loop through patients and extract everything from the matricies
for(a in seq_along(patients)){
	pre <- which(dict$patient==patients[a] &
		     dict$type=='PRETREAT')-2
	post <- which(dict$patient==patients[a] & 
		      dict$type=='POSTSBRT')-2
	postfolf <- which(dict$patient==patients[a] & 
		       dict$type=='POSTFOLF')-2
	if(length(postfolf)<1){postfolf <- NA}
	ol_pre_post[a] <- olm[pre,post]
	m_pre_post[a] <- mm[pre,post]
	ol_pre_postfolf[a] <- olm[pre,postfolf]
        m_pre_postfolf[a] <- mm[pre,postfolf]
	ol_postsbrt_postfolf[a] <- olm[post,postfolf]
        m_postsbrt_postfolf[a] <- mm[post,postfolf]
	arm[a] <- as.character(
			dict$arm[which(dict$patient==patients[a] & 
			dict$type=='PRETREAT')]
				)
}
arm <- as.factor(arm)


## overlap
dir.create(paste0(path, 'Overlaps/'), recursive=TRUE)
r <- ol_pre_post[arm==1]
nr <- ol_pre_post[arm==2]
p <- t.test(r,nr,var.equal=TRUE)$p.value
pdf(paste0(path,'Overlaps/overlap-pre-postsbrt.pdf'), width=4.5,height=8,title="Pre-Post Overlap")
stripchart(ol_pre_post ~ arm,
          at=c(1.25,1.75),
          pch=19,
          vertical=TRUE,
          xlim=c(1,2),
          ylab='Pre-Post-SBRT Overlap',
          )
text(1.5,max(c(r,nr)),paste0("p=",as.character(round(p,4))))
dev.off()
r <- ol_pre_postfolf[arm==1]
nr <- ol_pre_postfolf[arm==2]
p <- t.test(r,nr,var.equal=TRUE)$p.value
pdf(paste0(path,'Overlaps/overlap-pre-postfolf.pdf'), width=4.5,height=8,title="Pre-Post Overlap")
stripchart(ol_pre_postfolf ~ arm,
          at=c(1.25,1.75),
          pch=19,
          vertical=TRUE,
          xlim=c(1,2),
          ylab='Pre-Post-FOLF Overlap',
          )
text(1.5,max(c(r,nr)),paste0("p=",as.character(round(p,4))))
dev.off()

r <- ol_postsbrt_postfolf[arm==1]
nr <- ol_postsbrt_postfolf[arm==2]
p <- t.test(r,nr,var.equal=TRUE)$p.value
pdf(paste0(path,'Overlaps/overlap-post-sbrt-post-folf.pdf'), width=4.5,height=8,title="Pre-Post Overlap")
stripchart(ol_postsbrt_postfolf ~ arm,
          at=c(1.25,1.75),
          pch=19,
          vertical=TRUE,
          xlim=c(1,2),
          ylab='Post-SBRT-Post-FOLF Overlap',
          )
text(1.5,max(c(r,nr)),paste0("p=",as.character(round(p,4))))
dev.off()

# Morisita
r <- m_pre_post[arm==1]
nr <- m_pre_post[arm==2]
p <- t.test(r,nr,var.equal=TRUE)$p.value
pdf(paste0(path,'Overlaps/morisita-pre-post-sbrt.pdf'), width=4.5,height=8,title="Pre-Post Overlap")
stripchart(m_pre_post ~ arm,
          at=c(1.25,1.75),
          pch=19,
          vertical=TRUE,
          xlim=c(1,2),
          ylab='Pre-Post-SBRT Morisita Overlap',
          )
text(1.5,max(c(r,nr)),paste0("p=",as.character(round(p,4))))
dev.off()

r <- m_pre_postfolf[arm==1]
nr <- m_pre_postfolf[arm==2]
p <- t.test(r,nr,var.equal=TRUE)$p.value
pdf(paste0(path,'Overlaps/morisita-pre-post-FOLF.pdf'), width=4.5,height=8,title="Pre-Post Overlap")
stripchart(m_pre_postfolf ~ arm,
          at=c(1.25,1.75),
          pch=19,
          vertical=TRUE,
          xlim=c(1,2),
          ylab='Pre-Post-FOLF Morisita Overlap',
          )
text(1.5,max(c(r,nr)),paste0("p=",as.character(round(p,4))))
dev.off()

r <- m_postsbrt_postfolf[arm==1]
nr <- m_postsbrt_postfolf[arm==2]
p <- t.test(r,nr,var.equal=TRUE)$p.value
pdf(paste0(path,'Overlaps/morisita-post-sbrt-post-folf.pdf'), width=4.5,height=8,title="Pre-Post Overlap")
stripchart(m_postsbrt_postfolf ~ arm,
          at=c(1.25,1.75),
          pch=19,
          vertical=TRUE,
          xlim=c(1,2),
          ylab='Post-SBRT-Post-FOLF Morisita Overlap',
          )
text(1.5,max(c(r,nr)),paste0("p=",as.character(round(p,4))))
dev.off()

