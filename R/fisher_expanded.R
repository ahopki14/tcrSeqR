#mat <- matrix(ncol=2,sample(1:1000),100)


w <- which(dict$patient=='7.011' & dict$type!='PDAC')
mat <- as.matrix(ds[,w]) 
s <- apply(mat,MARGIN=2,FUN=sum)
exp_clone <- function(x,s){
	tab <- rbind(x,s-x)
	fisher.test(tab)$p.value
}
tm <- proc.time()
p_vals <- apply(mat,MARGIN=1,FUN=exp_clone,s=s)
el<- proc.time()-tm
p_adj <- p.adjust(p_vals,method='BH')
