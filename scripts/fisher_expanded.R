#mat <- matrix(ncol=2,sample(1:1000),100)

patients <- levels(dict$patient)
out <- data.frame(p_adj=numeric(),ind=numeric(),patient=character())
num_exp <- numeric()
p <- matrix(nrow=nrow(ds),ncol=length(patients))
colnames(p) <- patients
for(a in seq_along(patients)){
	w <- which(dict$patient==patients[a] & dict$type!='PDAC')
	mat <- ds[,w] 
    rownames(mat) <- seq(nrow(mat))
	mat$p_adj <- rep(1,nrow(mat))
	mat <- as.matrix(mat)
	ind <- which((mat[,1]+mat[,2])>5) # at least 5 in both
	ind <- as.numeric(ind)
	s <- apply(mat[ind,],MARGIN=2,FUN=sum)
	tm <- proc.time()
	p_vals <- exp_clone(mat[ind,2],mat[ind,1])
	el<- proc.time()-tm
	p_adj <- p.adjust(p_vals,method='BH')
	mat[ind,3] <- p_adj
	num_exp[a] <- length(which(mat[,3]<0.05))
	names(num_exp)[a] <- patients[a]
	p[,a] <- mat[,3]	
}
resp <- sapply(names(num_exp),FUN=iseqr_lookup,dict=dict)
df <- data.frame(patient=names(num_exp),response=resp,num_exp=num_exp)
g <- iseqr_plot_factor(df,'num_exp','response',NA)
g + ylab('Number of Expanded Clones')
