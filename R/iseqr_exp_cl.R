iseqr_exp_cl <-
function(ds,dict,s1='PRE',s2='POST',category='type',by='patient',inc.all=FALSE){
    patients <- levels(dict[,by])
	if(length(patients)==0){stop("Is dict$patients a factor?")}
    out <- data.frame(p_adj=numeric(),ind=numeric(),patient=character())
    big_out <- list()
    num_exp <- numeric()
    p <- matrix(nrow=nrow(ds),ncol=length(patients))
    colnames(p) <- patients
    for(a in seq_along(patients)){
        w <- c(which(dict[,by]==patients[a] & (dict[,category]==s1)),
                which(dict[,by]==patients[a] & dict[,category]==s2))
        if(length(w)>2){stop(paste('Multiple matches for:',patients[a],s1,s2))}
        if(length(w)==2){
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
          p[,a] <- mat[,3]
          if(inc.all){
           big_out[[a]] <- mat
           names(big_out)[a] <- patients[a]
         }
       }else{
          mat <- NA
          big_out[[a]] <- NA
          num_exp[a] <- NA
       }
       names(num_exp)[a] <- patients[a]

    }
    out <- list(num_exp=num_exp,p=p,total=big_out)
    if(inc.all){out}else{
      tmp <- as.data.frame(out$num_exp)
      names(tmp) <- paste0('Number of Expanded Clones vs ', s1)
      tmp[,by] <- rownames(tmp)
      tmp[,category] <- rep(s2,nrow(tmp))
      tmp
    }
}
