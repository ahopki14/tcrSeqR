#calculate the change in any statistic over any comparison
delta_stats <- function(plot_ds,comps, metric, merge=T){
	patients <- unique(plot_ds$patient)
	name <- paste('Log2 Fold Change in',metric) 
	out <- data.frame(patient=character(), type=character(), name=numeric())
	for(a in seq(length(patients))){
		for(b in comps){
			fc <- filter(plot_ds, type==b[[2]] & patient==patients[a])[,metric] /
				filter(plot_ds, type==b[[1]] & patient==patients[a])[,metric]
			lfc <- log2(fc)
			if(length(lfc)>0){
				out <- rbind(out, data.frame(patient=patients[a],type=b[[2]],name=lfc))
		}
		}
	}
	names(out)[3] <- name
	if(merge){
	out <- merge(plot_ds, out, by=c('patient','type'), all=T)
	}
	out
}
