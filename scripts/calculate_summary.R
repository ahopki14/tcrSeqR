#sample level stats
cl <- as.data.frame(sapply(ds[,-c(1,2)],clonality))
names(cl)[1] <- 'Clonality'
t <- as.data.frame(sapply(ds[,-c(1,2)],sum))
names(t)[1] <- 'Sum of Total Unique Sequences'
r <- as.data.frame(sapply(ds[,-c(1,2)],function(x){length(x[x>0])}))
names(r)[1] <- 'Richness'

# overlap
ol <- numeric(length(levels(dict$patient)))
count <- 0
for(a in seq(length(levels(dict$patient)))){
  count <- count + 1
  tds <- ds[ ,which(dict$patient==levels(dict$patient)[a])]
  tdict <- dict[which(dict$patient==levels(dict$patient)[a]),]
  ol[count] <- overlap(tds[ ,which(tdict$type=='Pre')],tds[ ,which(tdict$type=='Post')])
  names(ol)[count] <- levels(dict$patient)[a]
}

stats <- cbind(r,t,cl)

pdf(paste0('/home/ahopkins/Documents/emj/ImmunoseqResults/R/plots/number-vs-richness.pdf'),height=6,width=6)
plot(dict[-c(1,2),'cells'],stats$Richness,
  xlab='Input Cell Number',
  ylab='Richness')
dev.off()

pdf(paste0('/home/ahopkins/Documents/emj/ImmunoseqResults/R/plots/number-vs-clonality.pdf'),height=6,width=6)
plot(dict[-c(1,2),'cells'],stats$Clonality,
  xlab='Input Cell Number',
  ylab='Clonality')
dev.off()

pdf(paste0('/home/ahopkins/Documents/emj/ImmunoseqResults/R/plots/age-vs-richness.pdf'),height=6,width=6)
plot(dict[-c(1,2),'age'],stats$Richness,
  xlab='Age at Diagnosis',
  ylab='Richness')
dev.off()

pdf(paste0('/home/ahopkins/Documents/emj/ImmunoseqResults/R/plots/age-vs-clonality.pdf'),height=6,width=6)
plot(dict[-c(1,2),'age'],stats$Clonality,
  xlab='Age at Diagnosis',
  ylab='Clonality')
dev.off()

pdf(paste0('/home/ahopkins/Documents/emj/ImmunoseqResults/R/plots/age-by-response.pdf'),height=8,width=4.5)
stripchart(dict$age[-c(1,2)] ~ factor(dict$response[-c(1,2)]),
	pch=19,
	vertical=TRUE,
	at=c(1.25,1.75),
	xlim=c(1,2),
	ylab='Age')
dev.off()
# 24% difference, p=0.13226
t.test(dict$age[dict$response=='NR' & dict$type=='Pre'],dict$age[dict$response=='R' & dict$type=='Pre']































