#load data and dictionary
exp_threshold <- 0.05 # p value cutoff for expanded clones

dir.create(paste0(path,'Expanded_Clones/'))
    for(a in seq(length(levels(as.factor(dict$patient))))){
        name <- levels(as.factor(dict$patient))[a]
        samples <- c(
                     which(dict$patient==name & dict$day_seq==1),
                     which(dict$patient==name & dict$day_seq==2)
                    )
        tds <- ds[,samples] # restrict to patient
        tds <- tds[tds[,1]>5 | tds[,2]>5, ] # restrict to clones in Pre OR Post
#        ttds[ttds==0] <- 1 # add a pseudocount to avoid div by 0
#        lfc <- log2(ttds[ ,2]/ttds[ ,1]) # Post/Pre (bigger=expanded)
        p_vals <- exp_clone(tds[,1],tds[,2])
        p_adj <- p.adjust(p_vals,method='BH')
        exp_loc <- which(p_adj <= exp_threshold)
        pct <- round(100*(length(exp_loc)/length(p_adj)),2)
        pdf(paste0(path,'Expanded_Clones/expanded_clones_fisher-',name,'.pdf'),
            width=8,height=8, title='immunoSeqR: Expanded Clones')
        hist(p_adj[p_adj<0.05], 
            main=paste0('Distribution of Adjusted P Values \n',
            as.character(length(exp_loc)),
            '(',as.character(pct),')'
            )
        )
        dev.off()
}
# text(tmp[1],tmp[4],'text',pos=3,offset=1)

