library(lattice)
patients <- levels(as.factor(patient))[1:6]
types <- levels(as.factor(type))[2:5]
responses <- c('NR','R','NR','NR','R','R')

# Loop through each sample type (pre, post, TIL, PDAC) and plot each vs all others.
# Do this for every patient. 
for(x in 1:4){
	for(y in seq(1:4)[-x]){
			for(a in seq_along(patients)){
				tmp_ds <- ds[patient==patients[a]]
				tmp_ds <- tmp_ds[tmp_ds[,x]!=0 | tmp_ds[,y]!=0,]
				p <- cor(tmp_ds[,x],tmp_ds[,y])
				ol <- overlap(tmp_ds[,x],tmp_ds[,y])
				r <- range(tmp_ds[,x],tmp_ds[,y])[2]*1.05
				dir.create(paste0('~/Documents/emj/ImmunoseqResults/R/plots/corr/',type[x+1],'-v-',type[y+1],'/'),recursive=TRUE)
				pdf(paste0('~/Documents/emj/ImmunoseqResults/R/plots/corr/',type[x+1],'-v-',type[y+1],'/',patients[a],'.pdf'),width=8.5,height=8)
					par(oma=c(1,4,1,1))	
					plot(tmp_ds[,x],tmp_ds[,y],pch=".",asp=1,
						main=paste0(patients[a],'(',responses[a],')'),
						xlim=c(1,r),ylim=c(1,r),
						xlab=type[x+1],#also a little hacky...
						ylab=type[y+1],
						log='xy',
						cex.lab=1.5,
						cex.main=2)
					grid()
					text(r/100,r,paste0("o=",round(ol,4)),cex=2)
				dev.off()
		}
	}
}
# Loop through each patient and make a heatmap of each sample's overlap with the others.
for(a in seq_along(patients)){
	tmp_ds <- ds[patient==patients[a]]
	tmp_ds <- tmp_ds[tmp_ds[,x]!=0 | tmp_ds[,y]!=0,]
	tmp_cor <- cor(tmp_ds)
	tmp_cor <- tmp_cor[c(1,3,2,4),c(1,3,2,4)]
	colnames(tmp_cor) <- c('PDAC','TIL','Pre','Post')
	rownames(tmp_cor) <- c('PDAC','TIL','Pre','Post')
	pdf(paste0('~/Documents/emj/ImmunoseqResults/R/plots/corr/cor_plot_',patients[a],'.pdf'),width=8,height=8)
print(
	levelplot(tmp_cor,
		regions=TRUE,
		col.regions = gray(100:0/100),
		xlab="",
		ylab="",
		main=paste0('Sample correlations for ',patients[a],'(',responses[a],')'))
)
	dev.off()
	
}
