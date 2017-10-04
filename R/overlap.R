# Overlap
# A basic overlap function
overlap <- function(x,y){
    shared_x <- x[x>0 & y >0]
    shared_y <- y[x>0 & y >0]
    shared_sum <- sum(c(shared_x,shared_y))
    total <- sum(c(x,y))
    ol <- shared_sum/(total) # Adaptive has a + 1 in the denom...
#
    ol
}
