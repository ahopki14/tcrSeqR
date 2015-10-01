#sample level stats
# load data and dictionary

cl <- as.data.frame(sapply(ds[,-c(1,2)],clonality))
names(cl)[1] <- 'Clonality'
t <- as.data.frame(sapply(ds[,-c(1,2)],sum))
names(t)[1] <- 'Sum of Total Unique Sequences'
r <- as.data.frame(sapply(ds[,-c(1,2)],function(x){length(x[x>0])}))
names(r)[1] <- 'Richness'

# overlap
dict$patient <- as.factor(dict$patient)
ol <- numeric(length(levels(dict$patient)))
m <- numeric(length(levels(dict$patient)))

for(a in seq(length(levels(dict$patient)))){
  tds <- ds[ ,which(dict$patient==levels(dict$patient)[a])]
  tdict <- dict[which(dict$patient==levels(dict$patient)[a]),]
  ol[a] <- overlap(tds[ ,which(tdict$type=='PRE')],tds[ ,which(tdict$type=='POST')])
  m[a] <- morisita(tds[ ,which(tdict$type=='PRE')],tds[ ,which(tdict$type=='POST')])
  names(ol)[a] <- levels(dict$patient)[a]
}

stats <- cbind(r,t,cl)
save(stats,ol,file='stats.Rda')

