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
#         load('~/Documents/emj/ImmunoseqResults/neoadjuvant_study/exp_gen/stats.Rda',envir=.GlobalEnv)
         assign('path','~/Documents/emj/ImmunoseqResults/neoadjuvant_study/exp_gen/plots/',envir=.GlobalEnv)
        }
        if(x=='adj_expgen'){print('Loading Adjuvant Study Expected Genomes')
         load('~/Documents/emj/ImmunoseqResults/adjuvant_study/exp_gen/dict.Rda',envir=.GlobalEnv)
         load('~/Documents/emj/ImmunoseqResults/adjuvant_study/exp_gen/ds_adj_agg_expgen.Rda',envir=.GlobalEnv)
#         load('~/Documents/emj/ImmunoseqResults/adjuvant_study/exp_gen/stats.Rda',envir=.GlobalEnv)
         assign('path','~/Documents/emj/ImmunoseqResults/adjuvant_study/exp_gen/plots/',envir=.GlobalEnv)
        }
        if(x=='pilot_expgen'){print('Loading Pilot Study Expected Genomes')
         load('~/Documents/emj/ImmunoseqResults/sampleExport.2014-07-31_10-10-24/exp_gen/dict.Rda',envir=.GlobalEnv)
         load('~/Documents/emj/ImmunoseqResults/sampleExport.2014-07-31_10-10-24/exp_gen/ds_pilot_agg_expgen.Rda',envir=.GlobalEnv)
#         load('~/Documents/emj/ImmunoseqResults/sampleExport.2014-07-31_10-10-24/exp_gen/stats.Rda',envir=.GlobalEnv)
         assign('path','~/Documents/emj/ImmunoseqResults/sampleExport.2014-07-31_10-10-24/exp_gen/plots/',envir=.GlobalEnv)
        }
        if(x=='sbrt_expgen'){print('Loading SBRT Study Expected Genomes')
         load('~/Documents/emj/ImmunoseqResults/sbrt_study/exp_gen/dict.Rda',envir=.GlobalEnv)
         load('~/Documents/emj/ImmunoseqResults/sbrt_study/exp_gen/ds_sbrt_agg_exp_gen.Rda',envir=.GlobalEnv)
#         load('~/Documents/emj/ImmunoseqResults/sbrt_study/exp_gen/stats.Rda',envir=.GlobalEnv)
         assign('path','~/Documents/emj/ImmunoseqResults/sbrt_study/exp_gen/plots/',envir=.GlobalEnv)
        }
        if(x=='baseline_expgen'){print('Loading Baseline Study Expected Genomes')
         load('/home/ahopkins/Documents/emj/ImmunoseqResults/sampleExport.2015-05-21_09-55-12/exp_gen/dict.Rda',envir=.GlobalEnv)
         load('/home/ahopkins/Documents/emj/ImmunoseqResults/sampleExport.2015-05-21_09-55-12/exp_gen/ds_agg_exp_gen.Rda',envir=.GlobalEnv)
#         load('/home/ahopkins/Documents/emj/ImmunoseqResults/sampleExport.2015-05-21_09-55-12/exp_gen/stats.Rda',envir=.GlobalEnv)
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
bleep <- function(){system('bleep')}

