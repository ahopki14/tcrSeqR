resp <- sapply(names(out$num_exp),FUN=iseqr_lookup,dict=dict)
df <- data.frame(patient=names(out$num_exp),response=resp,num_exp=out$num_exp)
g <- iseqr_plot_factor(df,'num_exp','response',NA)
g + ylab('Number of Expanded Clones')
