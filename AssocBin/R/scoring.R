##' Scorings
##' @title Scoring functions to choose splits
##' @description These functions define scores to evaluate candidate
##' splits along a single margin within partition.
##' @details Each of these functions accepts `vals`, an ordered
##' numeric vector containing the candidate splits within a bin and
##' the bin bounds all in increasing order. To restrict splitting,
##' they also accept `expn` and `minExp`, which provide the expected
##' count within the split and minimum value of this count,
##' respectively. Any split which produces an expected value less
##' than `minExp` (assuming a uniform density within the bin) is given
##' a score of zero.
##' @param vals numeric vector candidate splits and bounds
##' @param expn the expected number of points in the bin
##' @param minExp the minimum number of points allowed in a bin
##' @return A vector of scores.
##' @examples
##' vals <- c(2, 5, 12, 16, 19)
##' ## restricting the minExp changes output
##' chiScores(vals, 4, minExp = 0)
##' chiScores(vals, 4, minExp = 2)
##' ## same for the miScores
##' miScores(vals, 4, minExp = 0)
##' miScores(vals, 4, minExp = 2)
##' ## random scoring produces different output every time
##' randScores(vals, 4, minExp = 0)
##' randScores(vals, 4, minExp = 0)
##' @author Chris Salahub
##' @describeIn scorings A chi-squared statistic score
chiScores <- function(vals, expn, minExp = 0) {
    diffs <- diff(vals)
    n <- length(vals) - 2
    total <- c(diffs[1]-1, cumsum(diffs))
    h1 <- total[1:(n+1)] # length below
    h2 <- total[n+2] - h1 # length above
    d <- n/total[n+2] # density
    i <- 0:n # number below split i
    ni <- n - i # number above i
    scr <- (i - d*h1)^2/(h1*d) + (ni - d*h2)^2/(h2*d)
    scr[is.na(scr)] <- 0
    scr[pmin(expn*h1/total[n+2],
             expn*h2/total[n+2]) < minExp] <- 0 # minimum size limit
    scr
}
##' @describeIn scorings A mutual information score
miScores <- function(vals, expn, minExp = 0) {
    diffs <- diff(vals)
    n <- length(vals) - 2
    total <- c(diffs[1] - 1, cumsum(diffs))
    h1 <- total[1:(n+1)] # length below
    h2 <- total[n+2] - h1 # length above
    d <- n/total[n+2] # density
    i <- 0:n # number below point i
    ni <- n - i # number above i
    below <- (i/n)*log(i/(d*h1))
    above <- (ni/n)*log(ni/(d*h2)) # split expectation
    below[1] <- 0
    above[n+1] <- 0 # handle known zeros
    scr <- below + above
    scr[pmin(expn*h1/total[n+2],
             expn*h2/total[n+2]) < minExp] <- 0 # minimum size limit
    scr
}
##' @describeIn scorings A random score for random splitting
randScores <- function(vals, expn, minExp = 0) {
    diffs <- diff(vals)
    n <- length(vals) - 2
    scores <- runif(length(diffs))
    ## compute the expected counts at each split
    total <- c(diffs[1]-1, cumsum(diffs))
    h1 <- total[1:(n+1)] # length below
    h2 <- total[n+2] - h1 # length above
    ## if the difference is one, splitting creates bin with area 0
    scores[1] <- min(diffs[1]-1, scores[1])
    ## difference of zero here does the same
    scores[length(scores)] <- min(diffs[length(diffs)],
                                  scores[length(scores)])
    scores[pmin(expn*h1/total[n+2],
                expn*h2/total[n+2]) < minExp] <- 0 # minimum size limit
    scores
}
