#pulls out clones from a text file
find_clones <- function(path='~/Documents/emj/ImmunoseqResults/published_clones.txt',
    ds=ds){
    clones <- read.table(path)
    w <- numeric(length(clones))
    for(a in seq(nrow(clones))){
        tmp <- grep(clones[a,1],ds$aa)
        if(length(tmp)>0){w[a]<-tmp}else{w[a] <- NA}
    }
    #w <- apply(clones,MARGIN=2,grep,ds$aa)
    #cl <- clones$V1[!is.na(w)]
    w <- w[!is.na(w)]
    pub <- ds[w,]
    pub$aa <- factor(pub$aa)
    pub
}
