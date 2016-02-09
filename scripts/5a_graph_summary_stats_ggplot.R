plot_ds <- merge(dict,stats)
if(nrow(plot_ds)!=nrow(ds)){stop('Could not merge, check the column names')}

dir.create(paste0(path,'ggplot'))
# set the range of variable that go on the x axis
x_vals <- c('response','arm')
metrics <- names(stats)[names(stats)!='fn']
for(metric in metrics){
    for(x_val in x_vals){
        for(type in levels(plot_ds$type)){
		g <- ggplot(plot_ds[plot_ds$type==type,],aes_q(x=as.name(x_val),y=as.name(metric))) +
            geom_point() +
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

plot_ds <- plot_ds[plot_ds$type=='PDAC' | plot_ds$type=='FFPE',]
for(metric in metrics){
    type_plot <- ggplot(plot_ds,
                    aes_q(x=plot_ds$type,y=as.name(metric))) +
        geom_point() +
        xlab('') +
        theme_bw() +
        theme(legend.position='none') +
        geom_text(aes(label=patient),
           hjust=-0.5,size=2)
        plot_width=length(levels(plot_ds$type))
        ggsave(type_plot,file=paste0(path,'ggplot/',metric,'-by_type_tumors.pdf'),
                    width=plot_width+0.5,height=6,units='in')
}


common_fields <- names(plot_ds)[3:5]
tmp <- split(plot_ds,plot_ds$type)
df1 <- as.data.frame(tmp[1])
names(df1)[3:5] <- common_fields
df2 <- as.data.frame(tmp[2])
names(df2)[3:5] <- common_fields
df3 <- as.data.frame(tmp[3])
names(df3)[3:5] <-  common_fields

tmp2 <- merge(df1,df2,common_fields,all.x=TRUE)
tmp3 <- merge(tmp2,df3,all.x=TRUE)
comp_ds <- tmp3[,c(1,2,3,grep('Clonality|Richness|Sum',names(tmp3)))]
names(comp_ds)[1:5] <- names(plot_ds)[1:5]
