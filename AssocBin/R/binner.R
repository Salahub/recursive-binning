##' @title Wrapper for recursive binning
##' @description `binner` is an iterative implementation of a
##' recursive binary partitioning algorithm which accepts the
##' splitting and stopping functions that guide partitioning as
##' arguments.
##' @details `binner` creates a two-dimensional histogram of the
##' sample space of `x` and `y` by recursively splitting partitions of
##' the data using `splitter` until `stopper` indicates that all
##' partitions are not to be split. An optional argument `init` gives
##' the function applied to the first bin containing all points to
##' initialize the binning algorithm.
##' @param x numeric vector of the first variable to be binned
##' @param y numeric vector of the second variable to be binned
##' @param stopper function which accepts a list with elements
##' `x`, `y`, `bnds`, `expn`, and `n` and returns a logical indicating
##' whether a split should occur for the bin defined by that list
##' @param splitter function which accepts a list of lists with
##' elements `x`, `y`, `bnds`, `expn`, and `n` and returns a list
##' where each element is a list of two corresponding to a split of
##' the bin at that position in the original list
##' @param init function like `splitter` applied to the sole first
##' bin
##' @return A list of lists each with elements `x`, `y`, `bnds`,
##' `expn`, and `n`.
##' @examples
##' ## necessary set up
##' crits <- makeCriteria(depth >= 4, n < 10, expn <= 5)
##' stopFn <- function(bns) stopper(bns, crits)
##' spltFn <- function(bn) maxScoreSplit(bn, chiScores)
##' ## generate data
##' x <- sample(1:100)
##' y <- sample(1:100)
##' ## run binner
##' bins <- binner(x, y, stopper = stopFn, splitter = spltFn)
##' @author Chris Salahub
binner <- function(x, y, stopper, splitter, init = halfSplit) {
    ## initialize bin with all the data contained
    bin <- list(x = x, y = y, # x and y points
                 bnds = list(x = c(0, max(x, na.rm = TRUE)),
                             y = c(0, max(y, na.rm = TRUE))),
                 expn = length(x), # default expectation is n
                 n = length(x), depth = 0) # size, depth
    binList <- init(bin) # first split, otherwise score max fails
    stopStatus <- stopper(binList) # initialize logical vector

    while (any(!stopStatus)) { # check the stop criteria
        oldBins <- binList[stopStatus] # stopped bins
        oldStop <- stopStatus[stopStatus] # all TRUE
        newBins <- lapply(binList[!stopStatus], splitter) # split bins
        newBins <- unlist(newBins, recursive = FALSE) # simplify
        newStop <- stopper(newBins) # get stop values
        binList <- c(oldBins, newBins) # update list of bins
        stopStatus <- c(oldStop, newStop) # update stop status
    }

    binList # return the final list of bins
}
