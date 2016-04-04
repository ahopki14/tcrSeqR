require(immunoSeqR)

ds <- readRDS('~/Documents/emj/ImmunoseqResults/mega/ds_Adjuvant.Rds')

ds$aa <- as.character(ds$aa)
w <- grep("\\*",ds$aa)
if(length(w)>0){ds <- ds[-w,]} #remove clones with stop codons
paste('Removed ',length(w),' clones with stop codons')
w <- which(ds$aa =='')
if(length(w)>0){ds <- ds[-w,]} #remove clones with no translation
paste('Removed ',length(w),' clones with no translation')
paste0('Leaving ',nrow(ds),' clones \n')

ds$nt <- as.character(ds$nt)
len <- sapply(ds$nt,nchar)
w <- which(len<87)
md <- grep('aa|nt',names(ds))
#colSums(ds[w,-md]) # should be zero for all type!=PDAC

ind <- split(seq_len(nrow(ds)),ds$aa)
syn <- sapply(ind,length)

w_syn <- which(syn>1)
to_fix <- ind[w_syn]

gc()

to_toss <- vector()
idk <- list()
for(a in seq_along(to_fix)){
#	tds <- ds[to_fix[[a]],]
	n <- length(to_fix[[a]])
	tlen <- as.vector(sapply(ds$nt[to_fix[[a]]],nchar))
	if(n==1){stop('Only one value found')}
	if(n==2){
		if(tlen[1]!=tlen[2]){
			keep <- which.max(tlen)
			toss <- which.min(tlen)
			k.ind <- to_fix[[a]][keep]
			t.ind <- to_fix[[a]][toss]
			if(grepl(ds$nt[t.ind],ds$nt[k.ind])){
				ds[k.ind,-md] <- colSums(ds[c(k.ind,t.ind),-md])
				to_toss <- c(to_toss,t.ind)
			}
		}
	}
	if(n>2){
		w_toss <- which(tlen<max(tlen))
		for(toss in w_toss){
			keep <- grep(ds$nt[to_fix[[a]][toss]],ds$nt[to_fix[[a]][-toss]])
			if(length(keep)>1){idk[[a]] <- to_fix[[a]]}
			if(length(keep)==1){
				k.ind <- to_fix[[a]][keep]	
				t.ind <- to_fix[[a]][toss]
				ds[k.ind,-md] <- colSums(ds[c(k.ind,t.ind),-md])
				to_toss <- c(to_toss,t.ind)
			}
		}	
	}
}
ds <- ds[-to_toss,]
paste('Ignored ',length(idk),' Due to multiple matches')
head(idk)
paste('Collapsed',length(to_toss),'clones')

saveRDS(ds,'~/Desktop/ds_Adjuvant_fixed.Rds')


