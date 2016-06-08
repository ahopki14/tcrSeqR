# change to dir
setwd('~/Documents/emj/ImmunoseqResults/data/neoadjuvant/')

# load data and dictionary
ds <- readRDS('ds_agg.Rds')
dict <- readRDS('dict.Rds')

# Calculate per sample stats
if(any(grepl('nt',names(ds)))){stop('These calculations should be performed on the aggregated data')}
w <- grep('aa|syn',names(ds))
cl <- as.data.frame(sapply(ds[,-w],clonality))
names(cl)[1] <- 'Clonality'
t <- as.data.frame(sapply(ds[,-w],sum))
names(t)[1] <- 'Total Sequences'
r <- as.data.frame(sapply(ds[,-w],function(x){length(x[x>0])}))
names(r)[1] <- 'Richness'

stats <- cbind(cl,r,t)
stats$fn <- rownames(stats)

##save
saveRDS(stats,file='stats.Rds')

plot_ds <- merge(dict,stats)

#calculate change in metrics
plot_ds <- delta_stats(plot_ds,'PRE','POST','type',names(stats)[1:3],'patient',merge=TRUE)


# Calculate the overlap and morisitas
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

##save
save(olm,mm,file=paste0(path,'olm.Rda'))

# Calculate Expanded Clones
out <- iseqr_exp_cl(ds,dict,s1='PRE',s2='POST',category='type',by='patient',inc.all=FALSE)
plot_ds <- merge(plot_ds,exp_clones,by=c('patient','type'),all=TRUE)
saveRDS(plot_ds,'plot_ds.Rds')



##save
saveRDS(out,file='exp_clones.Rds')



