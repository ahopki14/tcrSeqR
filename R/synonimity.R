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
