library(scales)
#load data and dictionary

# set how many clones will be plotted
n <- 50
# Should only clones present in the tumor be plotted?
restrict_to_tumor <- FALSE
# Where to put the output
path <= '/home/ahopkins/Documents/emj/ImmunoseqResults/adjuvant_study/plots/'

resp <- dict$response[which(dict$type=='PRE')]
numr <- as.integer(summary(resp)['R'])
numnr <- as.integer(summary(resp)['NR'])
num <- length(resp)

# create plot layout matrix
if(numnr>=numr){
l <- matrix(c(seq(numnr),seq(numr)+numnr,rep(0,numnr-numr)),numnr,2)
} else {
l <- matrix(c(seq(numnr),rep(0,numr-numnr),seq(numr)+numnr),numr,2)
}

# plot
if(restrict_to_tumor){clones <- 'in_tumor'} else{clones <- 'all'}
dir.create(paste0(path,'Rankplots/'))
pdf(paste0(path,'Rankplots/',clones,'-',n,'.pdf'),width=15,height=7*max(numr,numnr), title='immunoSeqR Rank Plot')
layout(l,respect=TRUE)
for(b in c('NR','R')){
  tds <- ds[ ,which(dict$resp==b)]
  tdict <- dict[which(dict$resp==b), ] 
  
  for(a in seq(length(unique(tdict$patient)))){
  	name <- unique(tdict$patient)[a]
  	samples <- c(
  		which(tdict$patient==name & tdict$type=='PRE'),
  		which(tdict$patient==name & tdict$type=='POST'),
  		which(tdict$patient==name & tdict$type=='PDAC')
  	)
  	ttds <- tds[,samples]
  	ttds <- ttds[ttds[,1]>0 | ttds[,2]>0, ]
  	if(restrict_to_tumor){
  	  ttds <- ttds[ttds[,3]>0, ]  	
  	}
  	o <- overlap(ttds[,1],ttds[,2]) 
  	pre <- rank(-ttds[,1],ties.method='min') 
  	post <- rank(-ttds[,2],ties.method='min')
  	pre_b <- which(pre %in% 1:n)
  	post_b <- which(post %in% 1:n)
  	desc <- c(rep("Pre",length(pre_b)),rep("Post",length(post_b)))
  	desc <- as.factor(desc)
  	desc <- factor(desc, levels=rev(levels(desc)))
  	dat <- c(pre[pre_b],post[pre_b])
  
  	stripchart( #makes an empty chart to draw lines on
  		c(0,0)~factor(c('Pre','Post'),levels=c('Pre','Post')),
  		at=c(1.25,1.75),
  		ylim=c((3*n),1),
  		xlim=c(1.2,1.8),
  		vertical=TRUE,
  		col='white',
  		ylab='Rank',
  		frame.plot=FALSE,
  		main=paste0(name,'(',b,')')#, '\n Present in Tumor')
  		)
  	for(a in 1:n){
  		segments(1.2,pre[pre_b][a],1.7,post[pre_b][a],col=alpha('red',0.4),lend=0,lwd=0.1)
  		segments(1.3,pre[post_b][a],1.8,post[post_b][a],col=alpha('blue',0.4),lend=0,lwd=0.1)
	  }
  }
}
dev.off()


