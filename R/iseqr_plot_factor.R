iseqr_plot_factor <- function(plot_ds,metric,x_val,type=NA,hjust=-0.5,labels=TRUE,...){
  if(!is.na(type)){
    plot_ds <- plot_ds[plot_ds$type==type,]
    title <- type
  }else{title <- ''}
    g <- ggplot(plot_ds,aes_q(x=as.name(x_val),y=as.name(metric))) +
    xlab('') + ylab(gsub('\\.',' ', metric)) + 
    ggtitle(as.character(title)) +
    theme_bw() + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
    stat_summary(fun.data=mean_95,geom='errorbar',width=0.05,colour='red')+
    stat_summary(fun.data=mean_only,geom='errorbar',width=0.25,colour='red')+
    geom_point()
  if(labels){g + geom_text(aes(label=patient),colour='grey',hjust=hjust,size=1)}else{g}
}



setMethod("iseqr_plot_factor",signature="tcr", 
	  definition=function(plot_ds,metric,x_val,...){
	iseqr_plot_factor(as.data.frame(colData(plot_ds)),metric,x_val, ...)
	}
)

