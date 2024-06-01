##' Scorings
##' @title Scoring functions to choose splits
##' @description These functions define scores to evaluate candidate
##' splits along a single margin within partition.
##' @details Each of these functions accepts `boundss`, an ordered
##' numeric vector containing the candidate splits within a bin and
##' the bin bounds all in increasing order, and `nbelow` which gives
##' the count of points below each split. `n` is used to determine the
##' number of points above the split.
##' @param bounds numeric vector giving candidate split bounds in
##' increasing order
##' @param nbelow integer vector giving the number of points below
##' each candidate split
##' @param n the total number of points in the bin to be split
##' @return A vector of scores.
##' @examples
##' vals <- c(2, 5, 12, 16, 19)
##' chiScores(vals, 1:3, 3)
##' ## same for the miScores
##' miScores(vals, 1:3, 3)
##' ## random scoring produces different output every time
##' randScores(vals, 1:3, 3)
##' randScores(vals, 1:3, 3)
##' @author Chris Salahub
##' @describeIn scorings A chi-squared statistic score
chiScores <- function(bounds, nbelow, n) {
    total <- bounds[2:length(bounds)] - bounds[1]
    h1 <- total[1:(length(nbelow))] # length below
    h2 <- total[length(total)] - h1 # length above
    d <- n/total[length(total)] # density
    i <- nbelow # number below split i
    ni <- n - i # number above i
    scr <- (i - d*h1)^2/(h1*d) + (ni - d*h2)^2/(h2*d)
    scr[is.na(scr)] <- 0
    scr
}
##' @describeIn scorings A mutual information score
miScores <- function(bounds, nbelow, n) {
    total <- bounds[2:length(bounds)] - bounds[1]
    h1 <- total[1:(length(nbelow))] # length below
    h2 <- total[length(total)] - h1 # length above
    d <- n/total[length(total)] # density
    i <- nbelow # number below point i
    ni <- n - i # number above i
    below <- (i/n)*log(i/(d*h1))
    above <- (ni/n)*log(ni/(d*h2)) # split expectation
    below[is.na(below)] <- 0
    above[is.na(above)] <- 0 # areas with zero expectation
    scr <- below + above
    scr
}
##' @describeIn scorings A random score for random splitting
randScores <- function(bounds, nbelow, n) {
    scores <- runif(length(nbelow))
    scores
}
