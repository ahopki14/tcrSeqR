fisher <- function(x,s){ # x is the row of mat, s is the column sum of mat
    tab <- rbind(x,s-x)
    fisher.test(tab,alternative='less')$p.value
}
