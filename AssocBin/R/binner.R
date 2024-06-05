##' @title Many split recursive binning
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
##' @param init function like `splitter` applied to the first bin
##' @return A list of lists each with elements `x`, `y`, `bnds`,
##' `expn`, `n`, and `stopped`.
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
    bin <- makeBin(x = x, y = y, # x and y points
                   bnds = list(x = c(0, max(x, na.rm = TRUE)),
                               y = c(0, max(y, na.rm = TRUE))),
                   expn = length(x), # default expectation is n
                   n = length(x), depth = 0, # size, depth
                   stopped = FALSE) # premature stopping 
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

##' @title Single margin binning
##' @description `uniBinner` is an iterative implementation of a
##' recursive binary partitioning algorithm which accepts the
##' splitting and stopping functions that guide partitioning as
##' arguments and applies them to the margin `y` alone.
##' @details `binner` creates a one-dimensional histogram of `y` for
##' each categorical value of `x` by recursively splitting partitions
##' of the data using `splitter` until `stopper` indicates that all
##' partitions are not to be split.
##' @param x factor vector for the the first variable
##' @param y numeric vector of the second variable (to be split)
##' @param stopper function which accepts a list with elements
##' `x`, `y`, `bnds`, `expn`, and `n` and returns a logical indicating
##' whether a split should occur for the bin defined by that list
##' @param splitter function which accepts a list of lists with
##' elements `x`, `y`, `bnds`, `expn`, and `n` and returns a list
##' where each element is a list of two corresponding to a split of
##' the bin at that position in the original list
##' @return A list of lists each with elements `x`, `y`, `bnds`,
##' `expn`, `n`, and `stopped`.
##' @author Chris Salahub
uniBinner <- function(x, y, stopper, splitter) {
    xtab <- table(x)
    xbrks <- rbind(cumsum(c(0, xtab[-length(xtab)])),
                   cumsum(xtab))
    colnames(xbrks) <- names(xtab)
    ## initial bin splits data by category
    binList <- lapply(names(xtab),
                      function(lev) {
                          makeBin(x[x == lev], y[x == lev],
                                  bnds = list(x = xbrks[, lev],
                                              y = range(y) - c(1,0)),
                                  expn = sum(x == lev),
                                  n = sum(x == lev), depth = 1,
                                  stopped = FALSE)
                      })
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

##' @title Binning of categorical variable pairs
##' @description `catBinner` converts the cross-tabulation of two
##' categorical variables into bins which work with all of the
##' functionality on bins built into `AssocBin`.
##' @details As both variables are already categorical, `catBinner`
##' performs no splits and does not merge any categories by default.
##' @param x factor vector for the first categorical variable
##' @param y factor vector for the second categorical variable
##' @return A list of lists each with elements `x`, `y`, `bnds`,
##' `expn`, `n`, and `stopped`.
##' @author Chris Salahub
catBinner <- function(x, y) {
    n <- length(x) # computational set up
    xtab <- table(x)
    xbrks <- rbind(cumsum(c(0, xtab[-length(xtab)])),
                   cumsum(xtab))
    colnames(xbrks) <- names(xtab)
    xlevels <- names(xtab); nx <- length(xlevels)
    ytab <- table(y)
    ybrks <- rbind(cumsum(c(0, ytab[-length(ytab)])),
                   cumsum(ytab))
    colnames(ybrks) <- names(ytab)
    ylevels <- names(ytab); ny <- length(ylevels)
    binList <- vector(mode = "list", length = nx*ny)
    
    for (ii in seq_len(nx)) { # iterate through all combinations
        for (jj in seq_len(ny)) {
            inds <- (x == xlevels[ii]) & (y == ylevels[jj])
            binList[[(ii - 1)*ny + jj]] <-
                list(x = x[inds], y = y[inds],
                     bnds = list(x = xbrks[, ii],
                                 y = ybrks[, jj]),
                     expn = xtab[[ii]]*ytab[[jj]]/length(x),
                     n = sum(inds), depth = 1, stopped = TRUE)
        }
    }

    binList # return the final list of bins
}

##' @title Single split recursive binning
##' @description `singleBinner` is an iterative implementation of a
##' recursive binary partitioning algorithm which accepts the
##' splitting and stopping functions that guide partitioning as
##' arguments.
##' @details `singleBinner` creates a two-dimensional histogram of the
##' sample space of `x` and `y` by recursively splitting partitions of
##' the data using `splitter` until `stopper` indicates that all
##' partitions are not to be split. An optional argument `init` gives
##' the function applied to the first bin containing all points to
##' initialize the binning algorithm. Unlike `binner`, it does this by
##' splitting one bin at a time, and so accepts an argument to specify
##' exactly how many bins to produce.
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
##' @param maxK integer giving the number of bins where splitting
##' is stopped regardless of stop criteria
##' @return A list of lists each with elements `x`, `y`, `bnds`,
##' `expn`, `n`, and `stopped`.
##' @examples
##' ## necessary set up
##' crits <- makeCriteria(depth >= 4, n < 10, expn <= 5)
##' stopFn <- function(bns) stopper(bns, crits)
##' spltFn <- function(bn) rIntSplit(bn, minExp = 5)
##' ## generate data
##' x <- sample(1:100)
##' y <- sample(1:100)
##' ## run binner
##' bins <- singleBinner(x, y, stopper = stopFn, splitter = spltFn)
##' @author Chris Salahub
singleBinner <- function(x, y, stopper, splitter, init = halfSplit,
                         maxK = 5) {
    ## initialize bin with all the data contained
    bin <- makeBin(x = x, y = y, # x and y points
                   bnds = list(x = c(0, max(x, na.rm = TRUE)),
                               y = c(0, max(y, na.rm = TRUE))),
                   expn = length(x), # default expectation is n
                   n = length(x), depth = 0, # size, depth
                   stopped = FALSE) # premature stopping 
    N <- length(x)
    binList <- init(bin)
    K <- 2
    areas <- sapply(binList, function(el) el$expn)
    binList <- binList[order(areas, decreasing = TRUE)]
    areas <- areas[order(areas, decreasing = TRUE)]
    stopStatus <- stopper(binList) # initialize logical vector

    while (any(!stopStatus) & K < maxK) { # check the stop criteria
        continue <- !stopStatus
        chosen <- sample(seq_len(sum(continue)), size = 1,
                         prob = areas[continue])
        newBins <- splitter(binList[continue][[chosen]])
        areas <- areas[-chosen] # remove split bin
        binList <- binList[-chosen]
        stopStatus <- stopStatus[-chosen]
        ord <- order(c(newBins[[1]]$expn, newBins[[2]]$expn,
                       areas))
        areas <- c(newBins[[1]]$expn, newBins[[2]]$expn,
                   areas)[ord]
        binList <- c(newBins, binList)[ord]
        stopStatus <- c(stopper(newBins), stopStatus)[ord]
        K <- K + 1
    }

    binList # return the final list of bins
}
