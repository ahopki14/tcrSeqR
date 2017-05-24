
# load data and dictionary

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
#saveRDS(stats,file='stats.Rds')

plot_ds <- merge(dict,stats)

plot_ds <- iseqr_order(plot_ds, ds)

#calculate change in metrics
plot_ds <- delta_stats(plot_ds, comps, 'Richness')
plot_ds <- delta_stats(plot_ds, comps, 'Clonality')

plot_ds <- iseqr_order(plot_ds, ds)


# Calculate overlaps

plot_ds <- iseqr_morisita(plot_ds, comps, ds, merge=T)
plot_ds <- iseqr_order(plot_ds, ds)
