#set path to adaptive tsv files
path <- '~/Documents/emj/ImmunoseqResults/immunoSeqR/data/ex_tsv/'
setwd(path)
all_files <- list.files(pattern=".tsv")

# construct the dataset with iseqr_merge
ds <- iseqr_merge(all_files)



# Remove Empty Sequences
length(which(ds$aa==''))
ds <- ds[ds$aa!='',]

# Remove Sequences with stop codons
length(grep('\\*',ds$aa))
ds <- ds[grep('\\*',ds$aa,invert=TRUE),]

# aggregate the data
# this collapses synonymous nucleotide sequences
ds_agg <- iseqr_aggregate(ds,inc_nt=FALSE)

#moving the aa and nt columns to the end makes life easier
ds <- ds[,c(3:8,1,2)]
ds_agg <- ds_agg[,c(3:8,1,2)]

#Construct a dictionary
dict <- data.frame(fn =names(ds)[1:6],
                   patient=as.factor(c(rep("1",3),rep("2",3))),
                   type=factor(rep(c("Tumor","Pre","Post"),2),levels=c("Pre","Post","Tumor")),
                   response=factor(c(rep("R",3),rep("NR",3)))
                  )

#going forward, use on the aggregated dataset
ds <- ds_agg

# rows of dictionary correspond to cols of ds
names(ds)
names(ds)[-c(7,8)] == dict$fn

# so you can do this:
w <- which(dict$type=='Pre') # restrict to Pre-vac sample types
head(ds[,w])

# or this:
w <- which(dict$patient=='1') # restrict to patient 1
head(ds[,w])

# Clonality of first sample
clonality(ds[,1])

# quickly calculate for all samples
sapply(ds[,-c(7,8)],clonality)


# What about richness?
length(which(ds[,1]>0))


# overlap
overlap(ds[,1],ds[,2])
morisita(ds[,1],ds[,2])

