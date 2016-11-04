
#pulls out clones from a text file
find_clones <- function(path='~/Documents/emj/ImmunoseqResults/published_clones.txt',
    ds=ds){
    clones <- read.table(path)
    w <- numeric(length(clones))
    for(a in seq(nrow(clones))){
        tmp <- grep(clones[a,1],ds$aa)
        if(length(tmp)>0){w[a]<-tmp}else{w[a] <- NA}
    }
    #w <- apply(clones,MARGIN=2,grep,ds$aa)
    #cl <- clones$V1[!is.na(w)]
    w <- w[!is.na(w)]
    pub <- ds[w,]
    pub$aa <- factor(pub$aa)
    pub
}


# an R terminal bell (require a shell script or alias which makes a noise and exits)
bleep <- function(){
  system('bleep &')
  system('notify-send -t 3000 "immunoSeqR" "Done"')
}

# after subsetting a dataframe, this refoactors all of the factor fields
refactor <- function(df){
	classes <- sapply(df,class)
	for(a in seq(classes)){
		if(classes[a]=='factor'){
			warn <- character()
			ord <- levels(df[,a])
			levs <- unique(as.character(df[,a]))
			loc <- rep(NA,length(levs))
			if(!any(is.na(levs))){
				for(b in seq_along(levs)){
					loc[b] <- which(levs[b]==ord)
				}
					loc <- order(loc)
					df[,a] <-factor(as.character(df[,a]),levels=levs[loc])
			}else{warn <- c(warn,names(classes)[a])}
		}
	}
if(length(warn)>0){
	warning(paste0('The following columns were omitted due to NAs:',paste0(warn,collapse=',')))
}
df
}

#old, not needed anymore?
iseqr_check <- function(dict,ds,stats=NA,v=FALSE){
	if(ncol(ds)==nrow(dict)){ #for padded dictionary
		if(v){warning('Dictionary is padded (likely using an old version)')}
		check_dict <- all(names(ds)==dict$fn,na.rm=TRUE)
		if(!all(is.na(stats))){check_stats <- all(dict$fn[-c(1,2)]==stats$fn)}
		else{check_stats <- TRUE}
	}else if(ncol(ds)==nrow(dict)+2){ # for un-padded dictionary
		if(v){print('Dictionary has correct size...')}
		w <- grep('aa|nt|syn',names(ds))
    check_dict <- all(names(ds)[-w]==dict$fn)
		if(!all(is.na(stats))){
        check_stats <- all(dict$fn==stats$fn)}
        else{check_stats <- NA}	
	}
out <- list(check_dict,check_stats)
names(out) <- c('ds','stats')
out
}



iseqr_clean_ds <- function(ds){
	w <- grep('syn|aa|nt',names(ds))
	ds <- ds[rowSums(ds[,-w])>0,]
	ds
}

mean_sd <- function(x){
	m <- mean(x,na.rm=TRUE)
	s <- sd(x,na.rm=TRUE)
	c(y=m,ymin=m-s,ymax=m+s)
}

mean_95 <- function(x){
	m <- mean(x,na.rm=TRUE)
	n <- length(x)
	i <- qt(0.975,df=n-1)*(sd(x)/sqrt(n))
	c(y=m,ymin=m-i,ymax=m+i)
}

iseqr_summarize <- function(plot_ds,split_on,metric){
    ds <- split(plot_ds,plot_ds[,split_on])
    out <- matrix(nrow=length(ds),ncol=length(ds))
    rownames(out) <- names(ds)
    colnames(out) <- names(ds)
    for(a in names(ds)){
        for(b in names(ds)[names(ds)!=a]){
            out[a,b] <- t.test(ds[[a]][metric],ds[[b]][metric])$p.value
        }
    }
out
}

iseqr_lookup <- function(i,dict,i_col='patient',o_col='response'){
    w <- which(dict[,i_col]==i)
    o <- dict[w,o_col]
    o <- as.character(o)
    if(length(unique(o))!=1){stop(paste('Multiple matches?',paste(o,collapse=',')))}else{
        unique(o)
    }
}


shared <- function(x,y){
	if(length(x)!=length(y)){stop('Vectors must be same length')}
	w <- which(x>0)
	out <- length(which(y[w]>0))
	out
}


tert <- function(x){
	q <- quantile(x,c(1/3,2/3))
	w1 <- which(x<=q[[1]])
	w2 <- which(x<=q[[2]] & x>q[[1]])
	w3 <- which(x>q[[2]])
	out <- character(length=length(x))
	out[w1] <- 1
	out[w2] <- 2
	out[w3] <- 3
	out <- factor(out,levels=c('1','2','3'))
	out
}

is.clean <- function(ds){
  has.stop <- grep('\\*',ds$aa)
  has.no.trans <- which(ds$aa=='')
  if(length(has.stop)==0 & length(has.no.trans)==0){out <- TRUE}
  else
    out <- FALSE
    warning(paste0(length(has.stop),' Stop Codons and ',length(has.no.trans),' Untranslated'))
  out
}

delta_stats <- function(plot_ds,s1,s2,col,metrics,split_on,merge=FALSE){
  orig <- plot_ds
  col.n <- which(names(plot_ds)==col)
  w <- sort(c(grep(s1,plot_ds[,col]),grep(s2,plot_ds[,col])))
  plot_ds <- plot_ds[w,]
  items <- levels(plot_ds[,split_on])
  out <- data.frame(matrix(nrow=length(items),ncol=length(metrics)))
  names(out) <- paste('Log2 Fold Change in',metrics)
  metric_ind <- grep(paste(metrics,collapse='|'),names(plot_ds))
    for(a in seq_along(items)){
      tds <- plot_ds[plot_ds[,split_on]==items[a],c(col.n,metric_ind)]
      d <- which(tds[,1]==s1)
      n <- which(tds[,1]==s2)
      # compute the log fold change of (s2/s1)
      if(length(d)+length(n)==2){
      out[a,] <- as.numeric(log(tds[n,-1]/tds[d,-1],2))
      }
      rownames(out)[a] <- items[a]
    }
  out[,split_on] <- rownames(out)
  if(!merge){
  out
  } else{
    out[,col] <- rep(s2,nrow(out))
    orig$sort_order <- seq(nrow(orig))
    w <- which(is.na(out[,1]))
    if(length(w)>0){
      orig  <- merge(orig,out[-w,],by=c(split_on,col),all=TRUE)
    }else{orig  <- merge(orig,out,by=c(split_on,col),all=TRUE)}
    orig <- orig[orig$sort_order,]
    w <- which(names(orig)=='sort_order')
    orig <- orig[,-w] 
  }
}

