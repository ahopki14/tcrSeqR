#mat <- matrix(ncol=2,sample(1:1000),100)
out <- numeric(length=nrow(mat))
for(a in seq(nrow(mat))){ 
	tab <- rbind(mat[a,],apply(mat[-a,],MARGIN=2,FUN=sum))
	out[a] <- fisher.test(tab)$p.value
}

s <- apply(mat,MARGIN=2,FUN=sum)
exp_clone <- function(x,s){
	tab <- rbind(x,s-x)
	fisher.test(tab)$p.value
}
tm <- proc.time()
out <- apply(mat,MARGIN=1,FUN=exp_clone,s=s)
el<- proc.time()-tm

