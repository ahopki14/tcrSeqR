
adj <- readRDS('~/Desktop/adj.Rds')
nadj <- readRDS('~/Desktop/nadj.Rds')
adj.dict <- readRDS('~/Documents/emj/ImmunoseqResults/mega/dict_Adjuvant.Rds')
nadj.dict <- readRDS('~/Documents/emj/ImmunoseqResults/mega/dict_Neoadjuvant.Rds')

resp_nadj <- sapply(names(nadj),FUN=iseqr_lookup,dict=nadj.dict)
resp_adj <- sapply(names(adj),FUN=iseqr_lookup,dict=adj.dict)

plot_ds <- data.frame(patients = c(names(adj),names(nadj)))
plot_ds$exp_cl <- c(unlist(adj),unlist(nadj))
plot_ds$resp <- c(unlist(resp_adj),unlist(resp_nadj))
plot_ds$experiment <- c(rep('Adjuvant',length(adj)),
						rep('Neoadjuvant',length(nadj)))

g <- ggplot(plot_ds, aes(x=experiment,y=exp_cl,color=resp)) + 
	geom_point() + 
	theme_bw() + 
	xlab('') + 
	ylab('Number of Expanded Clones')
ggsave(g,file='~/Desktop/exp_comparison.pdf')
