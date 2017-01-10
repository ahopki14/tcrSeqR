# change to dir
setwd('~/Documents/emj/ImmunoseqResults/data/baseline/')

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

#for adjuvant and neoadjuvant
#plot_ds <- delta_stats(plot_ds,'PRE','POST','type',names(stats)[1:3],'patient',merge=TRUE)

#for sbrt
tmp2 <- delta_stats(plot_ds,'PRETREAT','POSTFOLF','type',names(stats)[1:3],'patient',merge=TRUE)
plot_ds <- delta_stats(plot_ds,'PRETREAT','POSTSBRT','type',names(stats)[1:3],'patient',merge=TRUE)

w <- which(plot_ds$type=='POSTFOLF')
m_ind_1 <- grep('Log2',names(plot_ds))
m_ind_2 <- grep('Log2',names(tmp2))

plot_ds[w,m_ind_1] <- tmp2[w,m_ind_2]

#for ipi
plot_ds <- merge(dict,stats)
tmp2 <- delta_stats(plot_ds,'PRE','POST3','type',names(stats)[1:3],'patient',merge=TRUE)
plot_ds <- delta_stats(plot_ds,'PRE','POST1','type',names(stats)[1:3],'patient',merge=TRUE)

w <- which(plot_ds$type=='POST3')
m_ind_1 <- grep('Log2',names(plot_ds))
m_ind_2 <- grep('Log2',names(tmp2))

plot_ds[w,m_ind_1] <- tmp2[w,m_ind_2]

#for pd1
tmp2 <- delta_stats(plot_ds,'PRE','POST3','type',names(stats)[2:4],'patient',merge=TRUE)
tmp3 <- delta_stats(plot_ds,'PDACPRE','PDACPOST','type',names(stats)[2:4],'patient',merge=TRUE)
plot_ds <- delta_stats(plot_ds,'PRE','POST2','type',names(stats)[2:4],'patient',merge=TRUE)

w <- which(plot_ds$type=='POST3')
m_ind_1 <- grep('Log2',names(plot_ds))
m_ind_2 <- grep('Log2',names(tmp2))

plot_ds[w,m_ind_1] <- tmp2[w,m_ind_2]

w <- which(plot_ds$type=='PDACPOST')
m_ind_1 <- grep('Log2',names(plot_ds))
m_ind_2 <- grep('Log2',names(tmp3))

plot_ds[w,m_ind_1] <- tmp3[w,m_ind_2]





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


# for adjuvant and neoadjuvant
#plot_ds <- merge(plot_ds,exp_clones,by=c('patient','type'),all=TRUE)
out <- iseqr_exp_cl(ds,dict,s1='PRETREAT',s2='POSTSBRT',category='type',by='patient',inc.all=FALSE)
exp_clones <- as.data.frame(out$num_exp)
exp_clones$patient <- rownames(exp_clones)

plot_ds[,'Number of Expanded Clones vs PRE'] <- rep(NA,nrow(plot_ds))
col <- grep('Number of Expanded Clones',names(plot_ds))
for(a in seq(nrow(exp_clones))){
  w <- which(plot_ds$patient==exp_clones$patient[a] & plot_ds$type==exp_clones$type[a])
  plot_ds[w,col] <- exp_clones[a,1]
}

# and again for POSTFOLF
out <- iseqr_exp_cl(ds,dict,s1='PRETREAT',s2='POSTFOLF',category='type',by='patient',inc.all=FALSE)
exp_clones <- out

col <- grep('Number of Expanded Clones',names(plot_ds))
for(a in seq(nrow(exp_clones))){
  w <- which(plot_ds$patient==exp_clones$patient[a] & plot_ds$type==exp_clones$type[a])
  plot_ds[w,col] <- exp_clones[a,1]
}


# for ipi
out1 <- iseqr_exp_cl(ds,dict,s1='PRE',s2='POST1',category='type',by='patient',inc.all=FALSE)
exp_clones <- out1

plot_ds[,'Number of Expanded Clones vs PRE'] <- rep(NA,nrow(plot_ds))
col <- grep('Number of Expanded Clones',names(plot_ds))
for(a in seq(nrow(exp_clones))){
  w <- which(as.character(plot_ds$patient)==as.character(exp_clones$patient[a]) &
      as.character(plot_ds$type)==exp_clones$type[a])
  plot_ds[w,col] <- exp_clones[a,1]
}

# and again for POST3
out2 <- iseqr_exp_cl(ds,dict,s1='PRE',s2='POST3',category='type',by='patient',inc.all=FALSE)
exp_clones <- out2

col <- grep('Number of Expanded Clones',names(plot_ds))
for(a in seq(nrow(exp_clones))){
  w <- which(as.character(plot_ds$patient)==as.character(exp_clones$patient[a]) &
      as.character(plot_ds$type)==exp_clones$type[a])
  plot_ds[w,col] <- exp_clones[a,1]
}



saveRDS(plot_ds,'plot_ds.Rds')



##save
saveRDS(out,file='exp_clones.Rds')



