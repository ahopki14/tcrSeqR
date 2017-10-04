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

