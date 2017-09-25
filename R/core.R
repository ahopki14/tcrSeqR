# merge the exported ImmunoSeq data into a data frame
# This loads Adaptive tsv files into a single data frame
iseqr_merge <- 
function(all_files,data='estimatedNumberGenomes',nucleotide='nucleotide',aminoAcid='aminoAcid'){
    start <- proc.time() # start the clock
#clean file names
  name <- gsub(".tsv","",all_files)
  name <- gsub("_","",name)
  name <- gsub("-","",name)
  name <- paste0('sample',name)
  cnames <- names(read.delim(all_files[1],sep='\t',nrow=1))
# picking which column in the tsv to use
  if(!all(cnames[1:2] == c(nucleotide,aminoAcid))){
        stop("Unexpected columns in file.")}
    data_col <- grep(data,cnames)
    if(length(data_col)!=1){stop('Unable to find data column.')}
# Verify that names are OK
    if(length(unique(name)) != length(name)){stop("Names do not appear to be unique")}
    nt_list <- vector()
    aa_list <- vector()
# read in the nucleotide and amino acid names only
    for(a in seq_along(all_files)){
        tcnames <- names(read.delim(all_files[a],sep='\t',nrow=1))
        if(!all(tcnames==cnames)){stop(paste0('Error in ',all_files[a],'. Check the headers'))}
        ds <- read.delim(
            all_files[a],
            sep='\t',
            colClasses=c('character','character',rep('NULL',length(cnames)-2))
        )
        ds <- data.frame(ds)
        nt_list <- c(nt_list,as.vector(ds[,nucleotide]))
        aa_list <- c(aa_list,as.vector(ds[,aminoAcid]))
        cat(c("done loading ", all_files[a],"\n"))
    }
# clock
    t <- round((proc.time()-start)[['elapsed']]/60,2)
    cat(paste0('Loaded ',length(all_files),' samples in ',t,' minutes\n'))
#
    list <- data.frame(nt=nt_list,aa=aa_list)
    names(list) <- c(nucleotide,aminoAcid)
    list <- unique(list)
    ds <- list
    col_classes <- c('character','character',rep('NULL',length(cnames)-2))
    col_classes[data_col] <- 'integer'
# now read in the counts, matching to the list made above
    for(a in seq_along(all_files)){
        tmp_ds <- read.delim(
            all_files[a],
            sep='\t',
            colClasses=col_classes
            )
        tmp_ds <- as.data.frame(tmp_ds[,c(nucleotide,data)])
        colnames(tmp_ds)[2] <- paste(name[a])
        ds <- merge(ds,tmp_ds,by=nucleotide,all=TRUE)
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
    names(ds)[which(names(ds)==nucleotide)] <- 'nt'
    names(ds)[which(names(ds)==aminoAcid)] <- 'aa'
    ds
}


# Aggregate
# ds is a data frame made with iseqr_merge (needs 'aa' and 'nt', others will be lost)
# ds_out is a data frame of the aggregated data
iseqr_aggregate <- function(ds){
	stopifnot(class(ds)=='tcr')
	stopifnot(assayNames(ds)=='tcr_nt')
	start <- proc.time()
	# restrict to productive and sanitize factor
	w_remove <- which(rownames(ds)=='')
	w_keep <- which(rownames(ds)!='')
	ds <- ds[w_keep,]
	#ds$aa <- as.factor(as.character(ds$aa))
	# make indicies and synonymity
	ind <- split(seq_len(nrow(ds)),rownames(ds))
	syn <- sapply(ind,length)
	#split out columns to be fixed
	f <- syn>1 # id clones which need to be collapsed
	to_fix <- as.vector(unlist(ind[f]))
	ok <- ds[-to_fix,]
	stopifnot(length(rownames(ok))==length(unique(rownames(ok))))
	to_fix <- ds[to_fix,]
	# do math
	gc()
	ind_to_fix <- split(seq_len(nrow(to_fix)),rownames(to_fix))
	sum_rows <- function(x){colSums(assay(to_fix)[unlist(ind_to_fix[x]),])}
	fixed_mat <- sapply(seq(length(ind_to_fix)),FUN=sum_rows)
	fixed_mat <- t(fixed_mat)
	rownames(fixed_mat) <- names(ind_to_fix)
	#assemble a SE of the fixed sequences
	fixed <- SummarizedExperiment(assays=list(tcr_nt=fixed_mat))
	rowData(fixed) <- DataFrame(nt='NA',aa=rownames(fixed_mat))
	ds_out <- rbind(ok, fixed)
	names(ds_out@assays$data) <- 'tcr'
	#Now deal with the dictionary
	dict <- iseqr_order(dict, ds,reorder=T)
	# Print sanity check
	t <- round((proc.time()-start)[['elapsed']]/60,2)
	cat(paste0('Collapsed ',sum(syn[syn>1]),' TCRs to ',length(syn[syn>1]),'\n',
		   'Left ',length(syn[syn==1]),' TCRs\n'))
	cat(paste0('Removed ',length(w_remove),' Non-productive TCRs\n'))
	stopifnot(nrow(ds_out) + sum(syn[syn>1]) - length(syn[syn>1]) == nrow(ds))
	cat(paste0('Completed in ',t,' minutes\n'))
	#
	ds_out
}

iseqr_make_tcr <- function(ds, dict){
	require(SummarizedExperiment)
	#reoder the dictionary to match the ds
	dict <- iseqr_order(dict, ds, reorder=T)
	#get metadata and data column locations
	w_m <- grep('aa|nt|syn', names(ds))
	w_d <- grep('aa|nt|syn', names(ds), invert=T)
	# make the SE object with the data
	if(!any(grepl('syn',names(ds)))){
		#if nt data
		ds_se <- SummarizedExperiment(assays=list(tcr_nt=as.matrix(ds[,w_d], rownames=F)))
				}else{
		#if agg data
		ds_se <- SummarizedExperiment(assays=list(tcr=as.matrix(ds[,w_d], rownames=F)))
	}
	# add some metadata
	rownames(ds_se) <- ds$aa
	rowData(ds_se) <- DataFrame(ds[,w_m])
	stopifnot(all(colnames(ds_se)==dict$fn))
	colData(ds_se) <- DataFrame(dict)
	colnames(ds_se) <- colnames(ds[,w_d])
	#return the SE object
	tcr(ds_se)
}
