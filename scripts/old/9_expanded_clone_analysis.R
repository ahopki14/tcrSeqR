#load data and dictionary
# Where to put the output
#path <- '/home/ahopkins/Documents/emj/ImmunoseqResults/adjuvant_study/plots/'
tr <- 2 # restrict reads>tr in tumor
exp_threshold <- 0.05 # p value cutoff for expanded clones
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
out <- matrix(NA, ncol=length(levels(as.factor((dict$patient)))),nrow=2)
out <- as.data.frame(out)
names(out) <- levels(as.factor(dict$patient))

dir.create(paste0(path,'Expanded_Clones/'))
pdf(paste0(path,'Expanded_Clones/expanded_clones_fisher.pdf'),width=15,height=7*max(numr,numnr), title='immunoSeqR: Expanded Clones')
layout(l,respect=TRUE)
for(b in c('NR','R')){
  tds <- ds[ ,which(dict$resp==b)]
  tdict <- dict[which(dict$resp==b), ]
  tdict$patient <- factor(tdict$patient)

    for(a in seq(length(levels(tdict$patient)))){
        name <- levels(tdict$patient)[a]
        samples <- c(
                     which(tdict$patient==name & tdict$type=='PRE'),
                     which(tdict$patient==name & tdict$type=='POST'),
                     which(tdict$patient==name & tdict$type=='PDAC')
                    )
        ttds <- tds[,samples] # restrict to patient
        ttds <- ttds[ttds[,1]>5 | ttds[,2]>5, ] # restrict to clones in Pre OR Post
#        ttds[ttds==0] <- 1 # add a pseudocount to avoid div by 0
#        lfc <- log2(ttds[ ,2]/ttds[ ,1]) # Post/Pre (bigger=expanded)
        p_vals <- exp_clone(ttds[,1],ttds[,2])
        p_adj <- p.adjust(p_vals,method='BH')
        exp_loc <- which(p_adj <= exp_threshold)
        pct_tumor_expanded <- round(100*(
					length(which(ttds[exp_loc,3] >= tr))/
					length(which(ttds[ ,3] >= tr))
					),2
				   )
        pct_expanded_tumor <- round(100*(
					length(which(ttds[exp_loc,3] >= tr))/
					length(exp_loc)
					),2
				   )
	out[levels(tdict$patient)[a]] <- c(pct_tumor_expanded,pct_expanded_tumor) 
        w <- which(ttds[ ,3] >= tr)
        if(length(w)>0){
        plot(p_adj[w],ttds[w,3],
             pch=19,cex=2,
             xlab='P Value',
             ylab='TCR Count in Tumor',
             main=paste0(name,' (',b,')','\n',pct_tumor_expanded,"% / ",pct_expanded_tumor,"%")
            ) 
	} else {plot.new()}
 }
}
dev.off()
# text(tmp[1],tmp[4],'text',pos=3,offset=1)

