


# For PD1
comps <- list(c('PRE','POST2'),c('PRE','POST3'),c('PDACPRE','PDACPOST'))
out <- iseqr_morisita(dict,comps)
plot_ds <-  merge(plot_ds,out,by=c('patient','type'))
 

#For Ipi

comps <- list(c('PRE','POST1'),c('PRE','POST3'))
out <- iseqr_morisita(dict,comps)
plot_ds <-  merge(plot_ds,out,by=c('patient','type'))



