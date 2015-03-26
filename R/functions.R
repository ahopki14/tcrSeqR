# merge the exported ImmunoSeq data into a data frame
iseqr_merge <- function(all_files){
	start <- proc.time() # start the clock
	name <- gsub(".tsv","",all_files)
	name <- gsub("_","",name)
	name <- gsub("-","",name)
	name <- paste0('sample',name)
	cnames <- names(read.delim(all_files[1],sep='\t',nrow=1))
	if(!all(cnames[1:3] == c('nucleotide','aminoAcid','count'))){stop("Unexpected columns in file.")}
# Verify that names are OK
	if(length(unique(name)) != length(name)){stop("Names do not appear to be unique")}
	nt_list <- vector()
	aa_list <- vector()
	for(a in seq_along(all_files)){
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
	for(a in seq_along(all_files)){
		tmp_ds <- read.delim(
			all_files[a],
			sep='\t',
			colClasses=c('character','character','integer',rep('NULL',length(cnames)-3))
			)
		tmp_ds <- data.frame(nt=tmp_ds$nucleotide,count=tmp_ds$count)
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
olScan <- function(x,y){
	ds <- data.frame(x=x,y=y)
	ds <- ds[order(-ds$x), ]
	ol <- vector()
	for(a in 10:nrow(ds)){
	  ol <- c(ol,overlap(ds$x[1:a],ds$y[1:a]))
	}
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
