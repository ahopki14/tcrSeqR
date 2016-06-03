pub <- find_clones(ds=ds)
write.table(t(dict[,c('patient','type','HLA.A')]),col.names=FALSE,
	    file=paste0(path,'viral_clones.csv'),sep=',')
write.table(pub,append=TRUE,file=paste0(path,'viral_clones.csv'),sep=',',col.names=FALSE)
