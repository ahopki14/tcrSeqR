iseqr_plot_metrics <- function(plot_ds,metric,x_val,type){	
	g <- ggplot(plot_ds[plot_ds$type==type,],
		aes_q(x=as.name(x_val),y=as.name(metric))) +
    geom_point() +
	xlab('') +
	ggtitle(as.character(type)) +
    theme_bw() +
    geom_text(aes(label=patient),colour='grey',hjust=-0.5,size=2)
	g
}


k <- list(
		iseqr_plot_metrics(plot_ds,metric,x_val,'PRE'),
		iseqr_plot_metrics(plot_ds,metric,x_val,'POST'),
		iseqr_plot_metrics(plot_ds,metric,x_val,'PDAC')
		)


l <- vector('list')
for(a in levels(plot_ds$type)){
	l[[a]] <- iseqr_plot_metrics(plot_ds,metric,x_val,a)
}
do.call(grid.arrange,c(l,ncol=length(levels(plot_ds$type))))
