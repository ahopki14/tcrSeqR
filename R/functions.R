# merge the exported ImmunoSeq data into a data frame
iseqr_merge <- function(all_files,use_expected_genomes=FALSE){
    start <- proc.time() # start the clock
	name <- gsub(".tsv","",all_files)
	name <- gsub("_","",name)
	name <- gsub("-","",name)
	name <- paste0('sample',name)
	cnames <- names(read.delim(all_files[1],sep='\t',nrow=1))
	if(!all(cnames[1:2] == c('nucleotide','aminoAcid'))){
        stop("Unexpected columns in file.")}
    if(use_expected_genomes){data_col <- grep('estimated',cnames,ignore.case=TRUE)}else{
        data_col <- grep('count',cnames)}
    if(length(data_col)!=1){stop('Unable to find data column.')}
# Verify that names are OK
    if(length(unique(name)) != length(name)){stop("Names do not appear to be unique")}
    nt_list <- vector()
    aa_list <- vector()
    for(a in seq_along(all_files)){
        tcnames <- names(read.delim(all_files[a],sep='\t',nrow=1))
        if(!all(tcnames==cnames)){stop(paste0('Error in ',all_files[a],'. Check the headers'))}
        ds <- read.delim(
            all_files[a],
            sep='\t',
            colClasses=c('character','character',rep('NULL',length(cnames)-2))
        ) 
        ds <- data.frame(ds)
        nt_list <- c(nt_list,as.vector(ds$nucleotide))
        aa_list <- c(aa_list,as.vector(ds$aminoAcid))
        cat(c("done loading ", all_files[a],"\n"))
    }
# clock
    t <- round((proc.time()-start)[['elapsed']]/60,2) 
    cat(paste0('Loaded ',length(all_files),' samples in ',t,' minutes\n'))
#
    list <- data.frame(nt=nt_list,aa=aa_list)
    list <- unique(list)
    ds <- list
    col_classes <- c('character','character',rep('NULL',length(cnames)-2))
    col_classes[data_col] <- 'integer'
    for(a in seq_along(all_files)){
        tmp_ds <- read.delim(
            all_files[a],
            sep='\t',
            colClasses=col_classes
            )
        tmp_ds <- data.frame(nt=tmp_ds$nucleotide,count=tmp_ds[,3])
        colnames(tmp_ds)[2] <- paste(name[a])
        ds <- merge(ds,tmp_ds,by='nt',all=TRUE)
        cat(c("done with", all_files[a],"\n"))
    }
    ds[is.na(ds)==TRUE] <- 0
# Print sanity check
    t <- round((proc.time()-start)[['elapsed']]/60,2) 
    cat("------------------\n")
    cat('Loaded the following data, please check against the raw files:\n')
    for(a in seq_along(colnames(ds))){
        cat(c(colnames(ds)[a],":", length(ds[,a][ds[,a]!=0]),"\n"))
    }
    cat(paste0('Completed in ',t,' minutes\n'))
    cat("------------------\n")
	#
    ds
}




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


# Aggregate
# ds is a data frame made with iseqr_merge (needs 'aa' and 'nt', others will be lost)
# ds_out is a data frame of the aggregated data
iseqr_aggregate <- function(ds, inc_nt=TRUE){
    start <- proc.time()
# restrict to productive and sanitize factor
    ds <- ds[ds$aa!='',] 
    ds$aa <- as.factor(as.character(ds$aa)) 
# make indicies and synonymity
    ind <- split(seq_len(nrow(ds)),ds$aa)
    syn <- sapply(ind,length)
# define data columns
    dc <- names(ds)[names(ds)!='nt' & names(ds)!='aa']
# do math
	gc()
    sum_rows <- function(x){sapply(ds[x,dc],FUN=sum)}
    a <- sapply(ind,FUN=sum_rows)
    b <- data.frame(t(a))
# put together ds_out
    rn <- data.frame(aa=rownames(b))
    rownames(b) <- seq_len(nrow(b))
    syn <- data.frame(syn)
    ds_out <- cbind(syn,rn,b)
    rownames(ds_out) <- seq_len(nrow(ds_out))
# combine synonymous nucleotide sequences into a single char (comma sep)
    if(inc_nt){
        paste_names <- function(x){toString(ds[x,'nt'])}
        nt <- sapply(ind,FUN=paste_names)
        nt <- data.frame(nt)
        ds_out <- cbind(nt, ds_out)
    }
# Print sanity check
    t <- round((proc.time()-start)[['elapsed']]/60,2) 
    cat(paste0('Collapsed ',length(syn[syn>1]),' TCRs\n', 'Left ',length(syn[syn==1]),' TCRs\n'))
    cat(paste0('Completed in ',t,' minutes\n'))
#
    ds_out
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

iseqr_load <- function(x){
    if(x=='nadj'){print('Loading Neoadjuvant Study')
    load('~/Documents/emj/ImmunoseqResults/neoadjuvant_study/dict_new.Rda',envir=.GlobalEnv)
    load('~/Documents/emj/ImmunoseqResults/neoadjuvant_study/ds_nadj_agg.Rda',envir=.GlobalEnv)
    load('~/Documents/emj/ImmunoseqResults/neoadjuvant_study/stats.Rda',envir=.GlobalEnv)
    assign('path','~/Documents/emj/ImmunoseqResults/neoadjuvant_study/plots/',envir=.GlobalEnv)
    }
    if(x=='adj'){print('Loading Adjuvant Study')
    load('~/Documents/emj/ImmunoseqResults/adjuvant_study/dict.Rda',envir=.GlobalEnv)
    load('~/Documents/emj/ImmunoseqResults/adjuvant_study/ds_adj_agg.Rda',envir=.GlobalEnv)
    load('~/Documents/emj/ImmunoseqResults/adjuvant_study/stats.Rda',envir=.GlobalEnv)
    assign('path','~/Documents/emj/ImmunoseqResults/adjuvant_study/plots/',envir=.GlobalEnv)
    }
    if(x=='pilot'){print('Loading Pilot Study')
    load('~/Documents/emj/ImmunoseqResults/sampleExport.2014-07-31_10-10-24/rerun/dict.Rda',envir=.GlobalEnv)
    load('~/Documents/emj/ImmunoseqResults/sampleExport.2014-07-31_10-10-24/rerun/ds_agg.Rda',envir=.GlobalEnv)
    load('~/Documents/emj/ImmunoseqResults/sampleExport.2014-07-31_10-10-24/rerun/stats.Rda',envir=.GlobalEnv)
    assign('path','~/Documents/emj/ImmunoseqResults/sampleExport.2014-07-31_10-10-24/rerun/plots/',envir=.GlobalEnv)
    }
    if(x=='sbrt'){print('Loading SBRT Study')
    load('~/Documents/emj/ImmunoseqResults/sbrt_study/dict.Rda',envir=.GlobalEnv)
    load('~/Documents/emj/ImmunoseqResults/sbrt_study/ds_sbrt_agg.Rda',envir=.GlobalEnv)
    load('~/Documents/emj/ImmunoseqResults/sbrt_study/stats.Rda',envir=.GlobalEnv)
    assign('path','~/Documents/emj/ImmunoseqResults/sbrt_study/plots/',envir=.GlobalEnv)
    }
    if(x=='baseline'){print('Loading Baseline Study')
    load('/home/ahopkins/Documents/emj/ImmunoseqResults/sampleExport.2015-05-21_09-55-12/dict.Rda',envir=.GlobalEnv)
    load('/home/ahopkins/Documents/emj/ImmunoseqResults/sampleExport.2015-05-21_09-55-12/ds_agg.Rda',envir=.GlobalEnv)
    load('/home/ahopkins/Documents/emj/ImmunoseqResults/sampleExport.2015-05-21_09-55-12/stats.Rda',envir=.GlobalEnv)
    assign('path','/home/ahopkins/Documents/emj/ImmunoseqResults/sampleExport.2015-05-21_09-55-12/plots/',envir=.GlobalEnv)
    }

#################################################### Expencted Genomes Data ########################################
    if(x=='nadj_expgen'){print('Loading Neoadjuvant Study Expected Genomes')
    load('~/Documents/emj/ImmunoseqResults/neoadjuvant_study/exp_gen/dict.Rda',envir=.GlobalEnv)
    load('~/Documents/emj/ImmunoseqResults/neoadjuvant_study/exp_gen/ds_nadj_agg_expgen.Rda',envir=.GlobalEnv)
    load('~/Documents/emj/ImmunoseqResults/neoadjuvant_study/exp_gen/stats.Rda',envir=.GlobalEnv)
    assign('path','~/Documents/emj/ImmunoseqResults/neoadjuvant_study/exp_gen/plots/',envir=.GlobalEnv)
    }
    if(x=='adj_expgen'){print('Loading Adjuvant Study Expected Genomes')
    load('~/Documents/emj/ImmunoseqResults/adjuvant_study/exp_gen/dict.Rda',envir=.GlobalEnv)
    load('~/Documents/emj/ImmunoseqResults/adjuvant_study/exp_gen/ds_adj_agg_expgen.Rda',envir=.GlobalEnv)
    load('~/Documents/emj/ImmunoseqResults/adjuvant_study/exp_gen/stats.Rda',envir=.GlobalEnv)
    assign('path','~/Documents/emj/ImmunoseqResults/adjuvant_study/exp_gen/plots/',envir=.GlobalEnv)
    }
    if(x=='pilot_expgen'){print('Loading Pilot Study Expected Genomes')
    load('~/Documents/emj/ImmunoseqResults/sampleExport.2014-07-31_10-10-24/exp_gen/dict.Rda',envir=.GlobalEnv)
    load('~/Documents/emj/ImmunoseqResults/sampleExport.2014-07-31_10-10-24/exp_gen/ds_pilot_agg_expgen.Rda',envir=.GlobalEnv)
    load('~/Documents/emj/ImmunoseqResults/sampleExport.2014-07-31_10-10-24/exp_gen/stats.Rda',envir=.GlobalEnv)
    assign('path','~/Documents/emj/ImmunoseqResults/sampleExport.2014-07-31_10-10-24/exp_gen/plots/',envir=.GlobalEnv)
    }
    if(x=='sbrt_expgen'){print('Loading SBRT Study Expected Genomes')
    load('~/Documents/emj/ImmunoseqResults/sbrt_study/exp_gen/dict.Rda',envir=.GlobalEnv)
    load('~/Documents/emj/ImmunoseqResults/sbrt_study/exp_gen/ds_sbrt_agg_exp_gen.Rda',envir=.GlobalEnv)
    load('~/Documents/emj/ImmunoseqResults/sbrt_study/exp_gen/stats.Rda',envir=.GlobalEnv)
    assign('path','~/Documents/emj/ImmunoseqResults/sbrt_study/exp_gen/plots/',envir=.GlobalEnv)
    }
    if(x=='baseline_expgen'){print('Loading Baseline Study Expected Genomes')
    load('/home/ahopkins/Documents/emj/ImmunoseqResults/sampleExport.2015-05-21_09-55-12/exp_gen/dict.Rda',envir=.GlobalEnv)
    load('/home/ahopkins/Documents/emj/ImmunoseqResults/sampleExport.2015-05-21_09-55-12/exp_gen/ds_baseline_agg_expgen.Rda',
        envir=.GlobalEnv)
    load('/home/ahopkins/Documents/emj/ImmunoseqResults/sampleExport.2015-05-21_09-55-12/exp_gen/stats.Rda',envir=.GlobalEnv)
    assign('path','/home/ahopkins/Documents/emj/ImmunoseqResults/sampleExport.2015-05-21_09-55-12/exp_gen/plots/',envir=.GlobalEnv)
    }
}

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
bleep <- function(){system('bleep &')}

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

iseqr_check <- function(ds,dict,stats=NA,v=FALSE){
	if(ncol(ds)==nrow(dict)){ #for padded dictionary
		if(v){print('Dictionary IS padded')}
		check_dict <- all(names(ds)==dict$fn,na.rm=TRUE)
		if(!all(is.na(stats))){check_stats <- all(dict$fn[-c(1,2)]==stats$fn)}
		else{check_stats <- TRUE}
	}else if(ncol(ds)==nrow(dict)+2){ # for un-padded dictionary
		if(v){print('Dictionary is NOT padded')}
		check_dict <- all(names(ds)[names(ds)!='aa' & names(ds)!='syn']==dict$fn)
		if(!all(is.na(stats))){check_stats <- all(dict$fn==stats$fn)}
        else{check_stats <- TRUE}	
	}
out <- list(check_dict,check_stats)
names(out) <- c('dict','stats')
out
}


iseqr_plot_metrics <- function(plot_ds,metric,x_val,type,sm=TRUE){
    # get the r^2
	l <- lm(plot_ds[which(plot_ds$type==type),metric] ~
            plot_ds[which(plot_ds$type==type),x_val])
    l <- summary(l)
    l <- l$r.squared
	# make the plot
	g <- ggplot(plot_ds[plot_ds$type==type,],
        aes_q(x=as.name(x_val),y=as.name(metric))) +
    geom_point() +
    xlab(bquote(R^2 == .(round(l,3)))) +
    ggtitle(as.character(type)) +
    theme_bw() + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + 
    geom_text(aes(label=patient),colour='grey',hjust=-0.5,size=1)
	if(sm){
		g <- g + geom_smooth(method=lm,alpha=0.1)
		}
	g
}

iseqr_plot_factor <- function(plot_ds,metric,x_val,type){
	if(!is.na(type)){
		plot_ds <- plot_ds[plot_ds$type==type,]
		title <- type
	}else{title <- ''}
    g <- ggplot(plot_ds,aes_q(x=as.name(x_val),y=as.name(metric))) +
    xlab('') +
    ggtitle(as.character(title)) +
    theme_bw() + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
    geom_text(aes(label=patient),colour='grey',hjust=-0.5,size=1) +
    stat_summary(fun.data=mean_95,geom='errorbar',width=0.05,colour='red')+
    geom_point()
	g
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

iseqr_summarize <- function(ds,split_on,metric){
    ds <- split(ds,ds[,split_on])
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

fisher <- function(x,s){ # x is the row of mat, s is the column sum of mat
    tab <- rbind(x,s-x)
    fisher.test(tab,alternative='greater')$p.value
}

exp_clone <- function(x,y){ # x and y are vectors of counts in the samples compared
    mat <- as.matrix(cbind(x,y))
    s <- apply(mat,MARGIN=2,FUN=sum)
    p_vals <- apply(mat,MARGIN=1,FUN=fisher,s=s)
    p_vals
}



iseqr_exp_cl <- function(ds,dict,s1='PRE',s2='POST',category='type',by='patient'){
    patients <- levels(dict[,by])
    out <- data.frame(p_adj=numeric(),ind=numeric(),patient=character())
	big_out <- list()
    num_exp <- numeric()
    p <- matrix(nrow=nrow(ds),ncol=length(patients))
    colnames(p) <- patients
    for(a in seq_along(patients)){
        w <- c(which(dict[,by]==patients[a] & (dict[,category]==s1)),
                which(dict[,by]==patients[a] & dict[,category]==s2))
        if(length(w)!=2){stop(paste('Multiple matches for:',patients[a],s1,s2))}
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
		big_out[[a]] <- mat
		names(big_out)[a] <- patients[a]
    }
    out <- list(num_exp=num_exp,p=p,total=big_out)
    out
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
