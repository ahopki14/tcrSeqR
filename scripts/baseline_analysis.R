load('/home/ahopkins/Documents/emj/ImmunoseqResults/sampleExport.2015-05-21_09-55-12/ds_agg.Rda')
load('/home/ahopkins/Documents/emj/ImmunoseqResults/sampleExport.2015-05-21_09-55-12/dict.Rda')
cl <- sapply(ds_agg[,-c(1,2)],clonality)
cl <- as.vector(cl)

#clonality plot
day <- dict$day[-c(1,2)]
patient <- as.factor(dict$patient[-c(1,2)])
pdf('/home/ahopkins/Documents/emj/ImmunoseqResults/R/baseline_plots/clonality.pdf',height=6,width=10)
trellis.par.set(strip.background=list(col="lightgrey"))
xyplot(cl ~ day | patient, 
      layout = c(3,1),
      xlab='Date',
      ylab='Clonality',
      scales=list(alternating=FALSE))
dev.off()

# Clonalities
cl_stdevs <- rep(0,3)
cl_avg <- rep(0,3)
for(i in 1:3){
  cl_stdevs[i] <- sd(cl[patient==i])
  cl_avg[i] <- mean(cl[patient==i])
  }
mean(cl_stdevs) # 0.008528761
mean(cl_stdevs/cl_avg) #  0.1091


# Richness plot
day <- dict$day[-c(1,2)]
patient <- as.factor(dict$patient[-c(1,2)])
r <- sapply(ds_agg,function(x){length(x[x!=0])})
r <- as.vector(r)
r <- r[-c(1,2)]
pdf('/home/ahopkins/Documents/emj/ImmunoseqResults/R/baseline_plots/richness.pdf',height=6,width=10)
trellis.par.set(strip.background=list(col="lightgrey"))
 xyplot(r ~ day | patient,
      layout = c(3,1),
      xlab='Date',
      ylab='Richness',
      scales=list(alternating=FALSE))
dev.off()

#Richnesses 
r_stdevs <- rep(0,3)
r_avg <- rep(0,3)
for(i in 1:3){
  r_stdevs[i] <- sd(r[patient==i])
  r_avg[i] <- mean(r[patient==i])
  }
mean(r_stdevs) #  32154.03
mean(r_stdevs/r_avg) #  0.1267


# Total Reads
day <- dict$day[-c(1,2)]
patient <- as.factor(dict$patient[-c(1,2)])
tot <- sapply(ds_agg[,-c(1,2)],sum)
tot <- as.vector(tot)
pdf('/home/ahopkins/Documents/emj/ImmunoseqResults/R/baseline_plots/total_reads.pdf',height=6,width=10)
trellis.par.set(strip.background=list(col="lightgrey"))
 xyplot(tot ~ day | patient,
      layout = c(3,1),
      xlab='Date',
      ylab='Total Reads',
      scales=list(alternating=FALSE))
dev.off()

# Reads vs. Richness
pdf('/home/ahopkins/Documents/emj/ImmunoseqResults/R/baseline_plots/richness-vs-total_reads.pdf',height=10,width=10)
plot(tot,r,
  xlab='Total Reads',
  ylab='Richness',
)
dev.off()

# Reads vs. Clonality
pdf('/home/ahopkins/Documents/emj/ImmunoseqResults/R/baseline_plots/clonality-vs-total_reads.pdf',height=10,width=10)
plot(tot,cl,
  xlab='Total Reads',
  ylab='Clonality',
)
dev.off()

# Richness vs. Clonality
pdf('/home/ahopkins/Documents/emj/ImmunoseqResults/R/baseline_plots/richness-vs-clonality.pdf',height=10,width=10)
plot(r,cl,
  xlab='Richness',
  ylab='Clonality'
)
dev.off()


# 8x8 overlap heatmap
o <- matrix(data=rep(NA,64),nrow=8)
for(a in 2:3){
  o <- matrix(data=rep(NA,64),nrow=8)
  w <- which(dict$patient==a)
  ds <- ds_agg[,w]
  for(x in seq(8)){
    for(y in seq(8)[-x]){
      o[x,y] <- overlap(ds[ ,x],ds[ ,y])
      o[y,x] <- o[x,y] 
    }
  }
  save(o,file=paste0('patient_',a,'_overlap_matrix.Rda'))
 }

for(a in 1:3){
load(paste0('/home/ahopkins/Documents/emj/ImmunoseqResults/sampleExport.2015-05-21_09-55-12/','patient_',a,'_overlap_matrix.Rda'))

pdf(paste0('/home/ahopkins/Documents/emj/ImmunoseqResults/R/baseline_plots/overlap_patient_',a,'.pdf'),height=10,width=10)
print(
levelplot(o,
	regions=TRUE,
	#at=seq(0.25,0.65,length=100),
	col.regions = gray(100:0/100),
	xlab="",
	ylab="",
	main=paste0('Overlap (Donor ',a,')')
	)
)		
dev.off()
}

all_ol <- vector()
for(a in 1:3){
load(paste0('/home/ahopkins/Documents/emj/ImmunoseqResults/sampleExport.2015-05-21_09-55-12/','patient_',a,'_overlap_matrix.Rda'))
all_ol <- c(all_ol,o[1,])
}
all_ol[c(1,9,17)] <- 1
day <- dict$day[-c(1,2)]
patient <- as.factor(dict$patient[-c(1,2)])
pdf('/home/ahopkins/Documents/emj/ImmunoseqResults/R/baseline_plots/overlap.pdf',height=6,width=10)
trellis.par.set(strip.background=list(col="lightgrey"))
xyplot(all_ol ~ day | patient, 
      layout = c(3,1),
      xlab='Date',
      ylab='Overlap with First Blood Draw',
      scales=list(alternating=FALSE))
dev.off()
#Richnesses 
ol_stdevs <- rep(0,3)
ol_avg <- rep(0,3)
for(i in 1:3){
  ol_stdevs[i] <- sd(all_ol[patient==i],na.rm=TRUE)
  ol_avg[i] <- mean(all_ol[patient==i],na.rm=TRUE)
  }
mean(ol_stdevs) #  0.0205294
mean(ol_stdevs/ol_avg) #  0.05375939


		
		
