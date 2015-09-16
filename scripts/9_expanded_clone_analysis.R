#load data and dictionary
# Where to put the output
path <- '/home/ahopkins/Documents/emj/ImmunoseqResults/sampleExport.2014-07-31_10-10-24/rerun/'
tr <- 5 # restrict reads>tr in tumor
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
dir.create(paste0(path,'Expanded_Clones/'))
pdf(paste0(path,'Expanded_Clones/expanded_clones.pdf'),width=15,height=7*max(numr,numnr), title='immunoSeqR: Expanded Clones')
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
        ttds <- ttds[ttds[,1]>0 | ttds[,2]>0, ] # restrict to clones in Pre OR Post
        ttds[ttds==0] <- 1 # add a pseudocount to avoid div by 0
        lfc <- log2(ttds[ ,2]/ttds[ ,1]) # Post/Pre (bigger=expanded)
        #pdf
        w <- which(ttds[ ,3] > tr)
        plot(lfc[w],ttds[w,3],
             pch=20,
             xlab='log2(Fold Change Pre to Post)',
             ylab='TCR Count in Tumor',
             main=paste0(name,' (',b,')')
            ) 
 }
}
dev.off()


