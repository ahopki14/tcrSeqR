# Simsons Index
simpson <- function(x){
    x <- x/sum(x)
    x <- x^2
    l <- sum(x)
    l
}
