iseqr_plot_metrics <- function(plot_ds,x_vals){	
	metrics <- names(stats)[names(stats)!='fn']
    #plots <- vector('list')
    for(metric in metrics){
        for(a in seq_along(x_vals)){
            for(type in levels(plot_ds$type)){
            g <- ggplot(plot_ds[plot_ds$type==type,],
                    aes_q(x=as.name(x_vals[a]),y=as.name(metric))) +
                geom_point() +
                ggtitle(as.character(type)) +
                theme_bw() +
                geom_text(aes(label=patient),colour='grey',
                            hjust=-0.5,size=2)
	        #plots[[metric]] <- g
			print(g)
			}
	    }
	}
}

