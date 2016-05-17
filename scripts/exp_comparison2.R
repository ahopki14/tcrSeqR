require(plyr)
path <- '~/Documents/emj/ImmunoseqResults/new/'

load('~/Desktop/nadj_plot_ds.Rda')
nadj <- plot_ds
load('~/Desktop/adj_plot_ds.Rda')
adj <- plot_ds
load('~/Desktop/sbrt_plot_ds.Rda')
sbrt <- plot_ds

plot_ds <- rbind.fill(nadj,adj,sbrt)





pdac <- iseqr_plot_factor(plot_ds,'Clonality','experiment','PDAC')
pre <- iseqr_plot_factor(plot_ds,'Clonality','experiment','PRE')
post <- iseqr_plot_factor(plot_ds,'Clonality','experiment','POST')
cl <- grid.arrange(pre,post,pdac,ncol=3)
ggsave(cl,file=paste0(path,'Clonality.pdf'),width=8,height=4)

pdac <- iseqr_plot_factor(plot_ds,'Richness','experiment','PDAC')
pre <- iseqr_plot_factor(plot_ds,'Richness','experiment','PRE')
post <- iseqr_plot_factor(plot_ds,'Richness','experiment','POST')
r <- grid.arrange(pre,post,pdac,ncol=3)
ggsave(r,file=paste0(path,'Richness.pdf'),width=8,height=4)

pdac <- iseqr_plot_factor(plot_ds,'Total Sequences','experiment','PDAC')
pre <- iseqr_plot_factor(plot_ds,'Total Sequences','experiment','PRE')
post <- iseqr_plot_factor(plot_ds,'Total Sequences','experiment','POST')
ts <- grid.arrange(pre,post,pdac,ncol=3)
ggsave(ts,file=paste0(path,'TotalSequences.pdf'),width=8,height=4)


# filter by responders only
w <- which(plot_ds$response == 'R')
plot_ds_r <- plot_ds[w,]

pdac <- iseqr_plot_factor(plot_ds_r,'Clonality','experiment','PDAC')
pre <- iseqr_plot_factor(plot_ds_r,'Clonality','experiment','PRE')
post <- iseqr_plot_factor(plot_ds_r,'Clonality','experiment','POST')
cl <- grid.arrange(pre,post,pdac,ncol=3)
ggsave(cl,file=paste0(path,'Clonality_r.pdf'),width=8,height=4)

pdac <- iseqr_plot_factor(plot_ds_r,'Richness','experiment','PDAC')
pre <- iseqr_plot_factor(plot_ds_r,'Richness','experiment','PRE')
post <- iseqr_plot_factor(plot_ds_r,'Richness','experiment','POST')
r <- grid.arrange(pre,post,pdac,ncol=3)
ggsave(r,file=paste0(path,'Richness_r.pdf'),width=8,height=4)

pdac <- iseqr_plot_factor(plot_ds_r,'Total Sequences','experiment','PDAC')
pre <- iseqr_plot_factor(plot_ds_r,'Total Sequences','experiment','PRE')
post <- iseqr_plot_factor(plot_ds_r,'Total Sequences','experiment','POST')
ts <- grid.arrange(pre,post,pdac,ncol=3)
ggsave(ts,file=paste0(path,'TotalSequences_r.pdf'),width=8,height=4)

