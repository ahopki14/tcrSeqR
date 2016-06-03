iseqr_load('adj')
dict$fn <-  gsub('\\.','',as.character(dict$fn))
stats$fn <- rownames(stats)
adj_stats <- merge(dict,stats,by='fn')
adj_stats$exp <- 'Adjuvant' 

iseqr_load('nadj')
stats$fn <- rownames(stats)
nadj_stats <- merge(dict,stats,by='fn')
nadj_stats$exp <- 'Neoadjuvant'
nadj_stats <- nadj_stats[,-c(9,10)]

iseqr_load('pilot')
stats$fn <- rownames(stats)
pilot_stats <- merge(dict,stats,by='fn')
pilot_stats$exp <- 'Pilot'


stats <- rbind(adj_stats,nadj_stats,pilot_stats)
save(stats,file='~/Desktop/stats_experiment_comparison.Rda')
library(ggplot2)

metrics <- names(stats)[c(9,10,11)]
for(a in metrics){
    g <- ggplot(stats[stats$type=='PDAC' & stats[,10]>125000,],
    aes_q(x=as.name('exp'),y=as.name(a))) +
        geom_point(aes(colour=response)) +
        xlab('') +
        theme_bw() + 
        geom_text(aes(label=patient),colour='grey',hjust=-0.5,size=2)
    ggsave(g,file=paste0('~/Desktop/',a,'-experiment_comparison_125k_restricted.pdf'))
}



