require(immunoSeqR)

dict <- readRDS('~/Documents/emj/ImmunoseqResults/mega/dict_Neoadjuvant.Rds')
ds_nt <- readRDS('~/Documents/emj/ImmunoseqResults/mega/ds_Neoadjuvant.Rds')

ds <- ds_nt
rm(ds_nt)
ds$aa <- as.character(ds$aa)
w <- grep("\\*",ds$aa)
if(length(w)>0){ds <- ds[-w,]} #remove clones with stop codons
w <- which(ds$aa =='')
if(length(w)>0){ds <- ds[-w,]} #remove clones with no translation


ds$nt <- as.character(ds$nt)
len <- sapply(ds$nt,nchar)
w <- which(len<87)
md <- grep('aa|nt',names(ds))
colSums(ds[w,-md]) # should be zero for all type!=PDAC

ind <- split(seq_len(nrow(ds)),ds$aa)
syn <- sapply(ind,length)

w_syn <- which(syn>1)
to_fix <- ind[w_syn]


to_toss <- vector()
for(a in seq_along(to_fix)){
	tds <- ds[to_fix[[a]],]
	tlen <- as.vector(sapply(tds$nt,nchar))
	if(nrow(tds)==1){stop('Only one value found')}
	if(nrow(tds)==2){
		if(tlen[1]!=tlen[2]){
			keep <- which.max(tlen)
			toss <- which.min(tlen)
			k.ind <- to_fix[[a]][keep]
			t.ind <- to_fix[[a]][toss]
			ds[k.ind,-md] <- colSums(tds[,-md])
			to_toss <- c(to_toss,t.ind)
		}
	}
	if(nrow(tds)>2){
		w <- which(tlen<max(tlen))
		for(b in w){
			#find match and collapse
			
		}	
	}
}
ds <- ds[-to_toss,]


