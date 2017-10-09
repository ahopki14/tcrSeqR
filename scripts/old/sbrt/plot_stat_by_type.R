#requires plot_ds from 5a


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


#plots

## Clonality
g1 <- ggplot(comp_ds,aes(x=PRETREAT.Clonality,y=POSTSBRT.Clonality)) +
            geom_point() +
            geom_text(aes(label=patient,colour=arm),
            hjust=-0.5,size=2)+
            geom_line(aes(x=c(0,0.45),y=c(0,0.45))) +
            xlim(0,0.45) + ylim(0,0.45) +
            xlab('Pretreatment Clonality') + ylab('Post SBRT Clonality')+
            theme_bw() +
            theme(legend.position='none') 
ggsave(g1,file=paste0(path,'/ggplot/Clonality_PRETREAT_vs_POSTSBRT.pdf'))  

g2 <- ggplot(comp_ds,aes(x=POSTSBRT.Clonality,y=POSTFOLF.Clonality)) +
            geom_point() +
            geom_text(aes(label=patient,colour=arm),
            hjust=-0.5,size=2)+
            geom_line(aes(x=c(0,0.45),y=c(0,0.45))) +
            xlim(0,0.45) + ylim(0,0.45) +
            xlab('Post SBRT Clonality') + ylab('Post FOLF Clonality') +
            theme_bw() +
            theme(legend.position='none')
ggsave(g2,file=paste0(path,'/ggplot/Clonality_POSTSBRT_vs_POSTFOLF.pdf'))

## Richness

g1 <- ggplot(comp_ds,aes(x=PRETREAT.Richness,y=POSTSBRT.Richness)) +
            geom_point() +
            geom_text(aes(label=patient,colour=arm),
            hjust=-0.5,size=2)+
            geom_line(aes(x=c(0,214000),y=c(0,214000))) +
            xlim(0,214000) + ylim(0,214000) +
            xlab('Pretreatment Richness') + ylab('Post SBRT Richness')+
            theme_bw() +
            theme(legend.position='none')
ggsave(g1,file=paste0(path,'/ggplot/Richness_PRETREAT_vs_POSTSBRT.pdf'))

g2 <- ggplot(comp_ds,aes(x=POSTSBRT.Richness,y=POSTFOLF.Richness)) +
            geom_point() +
            geom_text(aes(label=patient,colour=arm),
            hjust=-0.5,size=2)+
            geom_line(aes(x=c(0,214000),y=c(0,214000))) +
            xlim(0,214000) + ylim(0,214000) +
            xlab('Post SBRT Richness') + ylab('Post FOLF Richness') +
            theme_bw() +
            theme(legend.position='none')
ggsave(g2,file=paste0(path,'/ggplot/Richness_POSTSBRT_vs_POSTFOLF.pdf'))


