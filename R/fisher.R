#' fisher
#' 
#' Performs a fisher test using a contingency table where the top row is the
#' abundance of one clone and the bottom row is the total counts of all other
#' clones in that sample
#'
#' @param x A 1x2 matrix containing the counts of one clone (typically before
#' and after a treatment)
#' @param s A 1x2 matrix containing the sum of all receptors in each sample
#' clonality (FALSE), given by 1 - Entropy.
#' @return A p value derived from a one sided fisher test (testing the
#' hypothesis that the clone expanded from the left sample to the right sample)
#' @author Alexander Hopkins
fisher <- function(x,s){ # x is the row of mat, s is the column sum of mat
    tab <- rbind(x,s-x)
    fisher.test(tab,alternative='less')$p.value
}
