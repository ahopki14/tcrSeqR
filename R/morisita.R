# Morisita's overlap index
morisita <- function(x,y){
    stopifnot(length(x) == length(y))
    prod <- x*y
    m <- (2*(sum(prod))) / ((simpson(x) + simpson(y))*sum(x)*sum(y))
    m
}
