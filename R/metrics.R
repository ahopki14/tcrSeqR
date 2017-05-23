#Signal to Noise
ss2n <- function(x,y){
    mean_x <- apply(x,1,mean)
    sd_x <- apply(x,1,sd)
    mean_y <- apply(y,1,mean)
    sd_y <- apply(y,1,sd)
    s2n <- (mean_x - mean_y) / (sd_x + sd_y)
    s2n[is.na(s2n)==TRUE] <- 0 # enforce sd > 0
    s2n[is.infinite(s2n)==TRUE] <- 0
#
  s2n
}



# Overlap
# A basic overlap function
overlap <- function(x,y){
    shared_x <- x[x>0 & y >0]
    shared_y <- y[x>0 & y >0]
    shared_sum <- sum(c(shared_x,shared_y))
    total <- sum(c(x,y))
    ol <- shared_sum/(total) # Adaptive has a + 1 in the denom...
#
    ol
}

# Overlap Scan
# Calculates the overlap of two samples in bigger and bigger chuncks
olScan <- function(x,y,byN=10){
    tds <- data.frame(x=x,y=y)
    tds <- tds[tds[,1]!=0 | tds[,2]!=0,]
    tds <- tds[order(-tds$x), ]
    byV <- seq(10,nrow(tds),byN)
    ol <- mapply(byV,FUN=function(n){overlap(tds$x[1:n],tds$y[1:n])})
    ol
}

#Synonymity
# x is a factor (coerced to char vector) or char vector
# syn is a named list of synonymities
synonymity <- function(x){
    x <- as.character(x)
    ind <- split(seq_len(length(x)),x)
    syn <- sapply(ind,length)
#
    syn
}

# clonality
clonality <- function(x,entropy=FALSE){
    x <- x[x != 0]
    x <- x/sum(x)
    cl <- -sum(x*log(x))/log(length(x)) # Entropy
    if(!entropy){
    cl <- 1-cl
    }
    cl
}



# Simsons Index
simpson <- function(x){
    x <- x/sum(x)
    x <- x^2
    l <- sum(x)
    l
}

# Morisita's overlap index
morisita <- function(x,y){
    if(length(x) != length(y)){stop("Vectors must be same length")}
    prod <- x*y
    m <- (2*(sum(prod))) / ((simpson(x) + simpson(y))*sum(x)*sum(y))
    m
}

# Wrapper for Morisita
iseqr_morisita <- function(dict,comps, ds){
#check if ds and dict are in correct order
stopifnot(dict$fn == names(ds)[1:(length(names(ds))-2)])
#
patients <- unique(dict$patient)
out <- dict[,c('patient','type')]
out$morisita <- rep(NA, nrow(out))
for(a in seq(length(patients))){
	for(b in comps){
		if(length(which(dict$patient==patients[a] & (dict$type==b[1]|dict$type==b[2]))) == 2){
			out[out$patient==patients[a] & out$type == b[2],'morisita'] <- 
			morisita(ds[,which(dict$patient==patients[a] & dict$type==b[1])],
				 ds[,which(dict$patient==patients[a] & dict$type==b[2])])
		}
	}

}
out
}
## Functions for Fisher Test

fisher <- function(x,s){ # x is the row of mat, s is the column sum of mat
    tab <- rbind(x,s-x)
    fisher.test(tab,alternative='less')$p.value
}


exp_clone <- function(x,y){ # x and y are vectors of counts in the samples compared
    mat <- as.matrix(cbind(x,y))
    s <- apply(mat,MARGIN=2,FUN=sum)
    p_vals <- apply(mat,MARGIN=1,FUN=fisher,s=s)
    p_vals
}


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

