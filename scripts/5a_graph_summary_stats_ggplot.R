plot_ds <- merge(dict,stats)
dir.create(paste0(path,'ggplot'))
# set the range of variable that go on the x axis
x_vals <- c('response','Number of Lymphoid Aggregates','Lymphoid Aggregate Tertile')
metrics <- names(stats)[names(stats)!='fn']
for(metric in metrics){
    for(x_val in x_vals){
        for(type in levels(dict$type)){
            g <- ggplot(plot_ds[plot_ds$type==type,],aes_q(x=as.name(x_val),y=as.name(metric))) +
                geom_point() +
#                xlab('') +
                ggtitle(as.character(type)) +
                theme_bw() +
                geom_text(aes(label=patient),colour='grey',
                        hjust=-0.5,size=2)
            plot_width= length(unique(plot_ds[,x_val]))
            ggsave(g,file=paste0(path,'ggplot/',metric,'-by_',x_val,'-',type,'.pdf'),
                    width=plot_width+0.5,height=6,units='in')
        }
    }
}
