plot_ds$arm <- rep(NA,nrow(plot_ds))


for(a in levels(plot_ds$patient)){
	plot_ds$arm[which(plot_ds$patient==a)] <- rep(arm$arm[which(arm$patient==a)],
							length(plot_ds$arm[which(plot_ds$patient==a)]))
}
