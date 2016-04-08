require(immunoSeqR)

ds <- readRDS('~/Documents/emj/ImmunoseqResults/mega/ds_Neoadjuvant.Rds')

ds$aa <- as.character(ds$aa)
w <- grep("\\*",ds$aa)
if(length(w)>0){ds <- ds[-w,]} #remove clones with stop codons
paste('Removed ',length(w),' clones with stop codons')
w <- which(ds$aa =='')
if(length(w)>0){ds <- ds[-w,]} #remove clones with no translation
paste('Removed ',length(w),' clones with no translation')
paste0('Leaving ',nrow(ds),' clones')

ds$nt <- as.character(ds$nt)
len <- sapply(ds$nt,nchar)
w <- which(len<87)
md <- grep('aa|nt',names(ds))
#colSums(ds[w,-md]) # should be zero for all type!=PDAC

ind <- split(seq_len(nrow(ds)),ds$aa)
syn <- sapply(ind,length)

w_syn <- which(syn>1)
to_fix <- ind[w_syn]

get_len <- function(x){as.vector(sapply(ds$nt[x],nchar))}
lens <- sapply(to_fix,get_len)
len_equiv <- function(x){min(lens[[x]])==max(lens[[x]])}
w <- sapply(seq(length(lens)),len_equiv)

to_fix <- to_fix[!w] # only groups where different lengths are present
lens <- lens[!w]
print(paste('Found ',length(lens), 'clones to collapse'))

gc()

to_toss <- vector()
idk <- list()
idk.c <- 1
nomatch <- list()
for(a in seq_along(to_fix)){
#	tds <- ds[to_fix[[a]],]
#	n <- length(to_fix[[a]])
	tlen <- lens[[a]]
	w.toss <- which(tlen<max(tlen))
	for(b in w.toss){
		w.match <- grep(ds$nt[to_fix[[a]][b]],ds$nt[to_fix[[a]]])
		w.match <- w.match[which(w.match!=b)] #so it doesn't match itself
		if(length(w.match)==0){
			nomatch <- c(nomatch,to_fix[[a]][b])
		}else{
			counts <- as.numeric(rowSums(ds[to_fix[[a]][w.match],-md]))
			w.big <- which(counts==max(counts))
			if(length(w.big)==1){
				k.ind <- to_fix[[a]][w.match[w.big]]
				t.ind <- to_fix[[a]][b]
				ds[k.ind,-md] <- colSums(ds[c(k.ind,t.ind),-md])
				to_toss <- c(to_toss,t.ind)
			}else{
					idk[idk.c] <- to_fix[[a]][b]
					idk.c <- idk.c + 1
				}
		}
	}	
}
ds <- ds[-to_toss,]
paste('Ignored ',length(idk),' Due to multiple matches')
paste('Collapsed',length(to_toss),'clones (',round(100*length(to_toss)/length(to_fix),2),'%)')
paste('Failed to collapse', length(nomatch), 'clones')

saveRDS(ds,'~/Documents/emj/ImmunoseqResults/mega/fixed/ds_Neodjuvant_fixed.Rds')


