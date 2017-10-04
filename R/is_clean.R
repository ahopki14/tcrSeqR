is_clean <- function(ds){
  has.stop <- grep('\\*',ds$aa)
  has.no.trans <- which(ds$aa=='')
  if(length(has.stop)==0 & length(has.no.trans)==0){out <- TRUE}
  else
    out <- FALSE
    warning(paste0(length(has.stop),' Stop Codons and ',length(has.no.trans),' Untranslated'))
  out
}
