pub <- find_clones(ds=ds)
write.table(t(dict),col.names=FALSE,file=paste0(path,'viral_clones.csv'),sep=',')
write.table(pub,append=TRUE,file=paste0(path,'viral_clones.csv'),sep=',')
