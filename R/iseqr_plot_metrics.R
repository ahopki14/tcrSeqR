#' iseqr_plot_metrics
#' 
#' Generated ggplot objects from a data.frame or tcr object
#'
#' @param plot_ds A data.frame or tcr object to be plotted 
#' @param metric The metric to be plotted (column name)
#' @param by The column to split the data by. This must be a numeric.
#' @param type A value to restric the data by (must correspond to a value in the
#' "type" field of plot_ds)
#' @param sm Logical indicating if trend line should be added with geom_smooth
#' @param hjust To be passed to geom_text
#' @param labels Logical indicating if labels in the "patient" field of plot_ds
#' should be printed
#' @return A ggplot object graphing the data
#' @author Alexander Hopkins
#' @export
iseqr_plot_metrics <- function(plot_ds,metric,by,type=NA,sm=TRUE,hjust=-0.5,labels=TRUE,...){
  if(!is.na(type)){
    plot_ds <- plot_ds[plot_ds$type==type,]
    title <- type
  }else{title <- ''}
    # get the r^2
  l <- lm(plot_ds[,metric] ~ plot_ds[,by])
    l <- summary(l)$r.squared
  # make the plot
  g <- ggplot(plot_ds, aes_q(x=as.name(by),y=as.name(metric))) +
    geom_point() +
    xlab(bquote(R^2 == .(round(l,3)))) +
    ggtitle(title) +
    theme_bw() + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
  if(labels){g <- g + geom_text(aes(label=patient),colour='grey',hjust=hjust,size=1)}
  if(sm){
    g <- g + geom_smooth(method=lm,alpha=0.1)
    }
  g
}

setMethod("iseqr_plot_metrics",signature="tcr", 
	  definition=function(plot_ds,metric,by,...){
	iseqr_plot_metrics(as.data.frame(colData(plot_ds)),metric,by, ...)
	}
)
