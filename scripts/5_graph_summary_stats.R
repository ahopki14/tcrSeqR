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


#comparative stats
patients <- levels(as.factor(dict$patient))
#make empty lists
ol_pre_post <- rep(NA,length(levels(dict$patient)))
names(ol_pre_post) <- as.character(levels(dict$patient))
ol_pre_tumor <- ol_pre_post
ol_post_tumor <- ol_pre_post
m_pre_post <- ol_pre_post
m_pre_tumor <- ol_pre_post
m_post_tumor <- ol_pre_post
resp <- ol_pre_post

#loop through patients and extract everything from the matricies
for(a in seq_along(patients)){
	pre <- which(dict$patient==patients[a] & dict$type=='PRE')-2
	post <- which(dict$patient==patients[a] & dict$type=='POST')-2
	tumor <- which(dict$patient==patients[a] & dict$type=='PDAC')-2
	ol_pre_post[a] <- olm[pre,post]
	m_pre_post[a] <- mm[pre,post]
	ol_pre_tumor[a] <- olm[pre,tumor]
        m_pre_tumor[a] <- mm[pre,tumor]
	ol_post_tumor[a] <- olm[post,tumor]
        m_post_tumor[a] <- mm[post,tumor]
	resp[a] <- as.character(dict$response[which(dict$patient==patients[a] & dict$type=='PRE')])
}
resp <- as.factor(resp)


## overlap
dir.create(paste0(path, 'Overlaps/'), recursive=TRUE)
r <- ol_pre_post[resp=='R']
nr <- ol_pre_post[resp=='NR']
p <- t.test(r,nr,var.equal=TRUE)$p.value
pdf(paste0(path,'Overlaps/overlap-pre-post.pdf'), width=4.5,height=8,
		title="Pre-Post Overlap")
stripchart(ol_pre_post ~ resp,
          at=c(1.25,1.75),
          pch=19,
          vertical=TRUE,
          xlim=c(1,2),
          ylab='Pre-Post Overlap',
          )
text(1.5,max(c(r,nr)),paste0("p=",as.character(round(p,4))))
dev.off()

r <- ol_pre_tumor[resp=='R']
nr <- ol_pre_tumor[resp=='NR']
p <- t.test(r,nr,var.equal=TRUE)$p.value
pdf(paste0(path,'Overlaps/overlap-pre-tumor.pdf'), width=4.5,height=8,title="Pre-Post Overlap")
stripchart(ol_pre_tumor ~ resp,
          at=c(1.25,1.75),
          pch=19,
          vertical=TRUE,
          xlim=c(1,2),
          ylab='Pre-Tumor Overlap',
          )
text(1.5,max(c(r,nr)),paste0("p=",as.character(round(p,4))))
dev.off()

r <- ol_post_tumor[resp=='R']
nr <- ol_post_tumor[resp=='NR']
p <- t.test(r,nr,var.equal=TRUE)$p.value
pdf(paste0(path,'Overlaps/overlap-post-tumor.pdf'), width=4.5,height=8,title="Pre-Post Overlap")
stripchart(ol_post_tumor ~ resp,
          at=c(1.25,1.75),
          pch=19,
          vertical=TRUE,
          xlim=c(1,2),
          ylab='Post-Tumor Overlap',
          )
text(1.5,max(c(r,nr)),paste0("p=",as.character(round(p,4))))
dev.off()

# Morisita
r <- m_pre_post[resp=='R']
nr <- m_pre_post[resp=='NR']
p <- t.test(r,nr,var.equal=TRUE)$p.value
pdf(paste0(path,'Overlaps/morisita-pre-post.pdf'), width=4.5,height=8,title="Pre-Post Overlap")
stripchart(m_pre_post ~ resp,
          at=c(1.25,1.75),
          pch=19,
          vertical=TRUE,
          xlim=c(1,2),
          ylab='Pre-Post Morisita Overlap',
          )
text(1.5,max(c(r,nr)),paste0("p=",as.character(round(p,4))))
dev.off()

r <- m_pre_tumor[resp=='R']
nr <- m_pre_tumor[resp=='NR']
p <- t.test(r,nr,var.equal=TRUE)$p.value
pdf(paste0(path,'Overlaps/morisita-pre-tumor.pdf'), width=4.5,height=8,title="Pre-Post Overlap")
stripchart(m_pre_tumor ~ resp,
          at=c(1.25,1.75),
          pch=19,
          vertical=TRUE,
          xlim=c(1,2),
          ylab='Pre-Tumor Morisita Overlap',
          )
text(1.5,max(c(r,nr)),paste0("p=",as.character(round(p,4))))
dev.off()

r <- m_post_tumor[resp=='R']
nr <- m_post_tumor[resp=='NR']
p <- t.test(r,nr,var.equal=TRUE)$p.value
pdf(paste0(path,'Overlaps/morisita-post-tumor.pdf'), width=4.5,height=8,title="Pre-Post Overlap")
stripchart(m_post_tumor ~ resp,
          at=c(1.25,1.75),
          pch=19,
          vertical=TRUE,
          xlim=c(1,2),
          ylab='Post-Tumor Morisita Overlap',
          )
text(1.5,max(c(r,nr)),paste0("p=",as.character(round(p,4))))
dev.off()

