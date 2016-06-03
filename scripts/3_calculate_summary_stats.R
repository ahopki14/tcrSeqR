#sample level stats
# load data and dictionary
dir.create(path)
cl <- as.data.frame(sapply(ds[,-c(1,2)],clonality))
names(cl)[1] <- 'Clonality'
t <- as.data.frame(sapply(ds[,-c(1,2)],sum))
names(t)[1] <- 'Sum of Total Unique Sequences'
r <- as.data.frame(sapply(ds[,-c(1,2)],function(x){length(x[x>0])}))
names(r)[1] <- 'Richness'

# overlap
#dict$patient <- as.factor(dict$patient)
#ol <- numeric(length(levels(dict$patient)))
#m <- numeric(length(levels(dict$patient)))

#for(a in seq(length(levels(dict$patient)))){
#  tds <- ds[ ,which(dict$patient==levels(dict$patient)[a])]
#  tdict <- dict[which(dict$patient==levels(dict$patient)[a]),]
#  ol[a] <- overlap(tds[ ,which(tdict$type=='PRE')],tds[ 
#,which(tdict$type=='POST')])
#  m[a] <- morisita(tds[ ,which(tdict$type=='PRE')],tds[ 
#,which(tdict$type=='POST')])
#  names(ol)[a] <- levels(dict$patient)[a]
#  names(m)[a] <- levels(dict$patient)[a]
#}

stats <- cbind(r,t,cl)
stats$fn <- rownames(stats)

#save(stats,ol,m,file='stats.Rda')


# This should probably be re-done with apply()...
nsamp <- length(names(ds))-2
olm <- data.frame(matrix(NA,nrow=nsamp,ncol=nsamp))
mm <- data.frame(matrix(NA,nrow=nsamp,ncol=nsamp))
for(a in 1:(nsamp)){
	for(b in 1:(nsamp)){
		if(a!=b){
			olm[a,b] <- overlap(ds[ ,a],ds[ ,b])
			mm[a,b] <- morisita(ds[ ,a],ds[ ,b])
		}
		colnames(olm)[b] <- names(ds)[b]
		rownames(olm)[a] <- names(ds)[a] 
		colnames(mm)[b] <- names(ds)[b]
                rownames(mm)[a] <- names(ds)[a]
	}
}
olm <- as.matrix(olm)
mm <- as.matrix(mm)

save(stats,olm,mm,file=paste0(path,'stats.Rda'))

