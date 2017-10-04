#' iseqr_plot_factor
#' 
#' Generated ggplot objects from a data.frame or tcr object, splitting the data
#' appropriately and adding error bars
#'
#' @param plot_ds A data.frame or tcr object to be plotted 
#' @param metric The metric to be plotted (column name)
#' @param by The column to split the data by. This must be a factor
#' @param type A value to restric the data by (must correspond to a value in the
#' "type" field of plot_ds)
#' @param hjust To be passed to geom_text
#' @param labels Logical indicating if labels in the "patient" field of plot_ds
#' should be printed
#' @return A ggplot object graphing the data
#' @author Alexander Hopkins
#' @export
iseqr_plot_factor <- function(plot_ds,metric,by,type=NA,hjust=-0.5,labels=TRUE,...){
  if(!is.na(type)){
    plot_ds <- plot_ds[plot_ds$type==type,]
    title <- type
  }else{title <- ''}
    g <- ggplot(plot_ds,aes_q(x=as.name(by),y=as.name(metric))) +
    xlab('') + ylab(gsub('\\.',' ', metric)) + 
    ggtitle(as.character(title)) +
    theme_bw() + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
    stat_summary(fun.data=mean_95,geom='errorbar',width=0.05,colour='red')+
    stat_summary(fun.data=mean_only,geom='errorbar',width=0.25,colour='red')+
    geom_point()
  if(labels){g + geom_text(aes(label=patient),colour='grey',hjust=hjust,size=1)}else{g}
}



setMethod("iseqr_plot_factor",signature="tcr", 
	  definition=function(plot_ds,metric,by,type,...){
	iseqr_plot_factor(as.data.frame(colData(plot_ds)),metric,by,type=type, ...)
	}
)

