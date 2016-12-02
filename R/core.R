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
        tmp_ds <- data.frame(nt=tmp_ds[,nucleotide],count=tmp_ds[,3])
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
    ds
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

