all_files <- list.files(pattern='.tsv')
for(a in all_files){
	ds <- read.csv(a,sep='\t')
	cols <- grep('nucleotide|aminoAcid|count|estimatedNumberGenomes',names(ds))
	ds <- ds[,cols]
	write.table(ds,file=paste0('mod/',a),sep='\t',quote=FALSE,row.names=FALSE)
}
