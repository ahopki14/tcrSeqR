#mat <- matrix(ncol=2,sample(1:1000),100)

patients <- levels(dict$patient)
out <- data.frame(p_adj=numeric(),ind=numeric(),patient=character())
for(a in seq_along(patients)){
	w <- which(dict$patient==patients[a] & dict$type!='PDAC')
	mat <- ds[,w] 
	mat <- mat[mat[,1]>5 & mat[,2]>5,]
	mat <- as.matrix(mat)
	ps_mat <- mat
	ps_mat[ps_mat==0] <- 1
	fc <- (log(ps_mat[,2]/ps_mat[,1],2))
	mat <- mat[fc>1,] 
	s <- apply(mat,MARGIN=2,FUN=sum)
	tm <- proc.time()
	p_vals <- exp_clone(mat[,1],mat[,2])
	el<- proc.time()-tm
	p_adj <- p.adjust(p_vals,method='BH')
	exp_cl <- length(p_adj[p_adj<0.01])
	exp_cl_pct <- round(100*exp_cl/length(p_adj),2)
	pdf(file=paste0(path,'Expanded_Clones/',patients[a],'-hist.pdf'),
		width=8,height=8) 
	hist(p_adj[p_adj<0.05],breaks=500,
		main=paste0(patients[a],'\n',exp_cl,' clones (',exp_cl_pct,'%)'),
		xlab='P values')
	dev.off()
	p_adj_ord <- p_adj[order(p_adj)]
	tout <- data.frame(p_adj=p_adj_ord, ind=seq(length(p_adj)), 
			   patient=rep(patients[a],length(p_adj))) 
	out <- rbind(out,tout)
}


x <- seq(0,0.05,length=1000)

tmp <- function(x,y){
	length(y[y<x])
}
for(a in seq_along(patients)){
#	r <-  dict$response[dict$patient==patients[a] & 
#			    dict$type=='PRE'][3]
	col <- 'red'
#	if(r=='R'){col='black'}
	y <- out$p_adj[out$patient==patients[a]]
	num <- sapply(x,FUN=tmp,y=y)
pdf(file=paste0(path,'Expanded_Clones/',
	patients[a],'-exp_plot.pdf'))
	plot(x,num,pch=20,col=col,main=patients[a],
		xlab='FDR Cutoff',
		ylab='Number of Clones')
dev.off()
}
