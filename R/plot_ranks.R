load('ds_agg.Rda')
tds <- ds_new[,c(4,6)]
tds <- tds[tds[,1]>0 | tds[,2]>0,]
Pre <- rank(-tds[,1])
Post <- rank(-tds[,2])
plot(Pre,Post,log='xy',col=alpha('black',0.3),pch=19)
b <- which(Pre %in% 1:10)
