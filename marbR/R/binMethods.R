## a simple constructor
## data has elements x and y with no NAs
makeBin <- function(data) {
    list(x = data$x, y = data$y,
         bnds = list(x = range(data$x),
                     y = range(data$y)), # boundaries
         area = prod(diff(range(data$x)), diff(range(data$y))),
         expn = length(data$x), # expected count is n
         n = length(data$x), depth = 0) # bin size, depth
}

##' bin splitter given the indices below and the split value:
##' if the indices are specified without a value, choose the
##' maximal value in these indices
splitBin <- function(bin, belowInds, mar, value) {
    below <- rep(FALSE, bin$n)
    below[belowInds] <- TRUE # handle 0 case
    lowbnds <- upbnds <- bin$bnds
    lowbnds[[mar]][2] <- value
    upbnds[[mar]][1] <- value # changing bounds
    uparea <- prod(sapply(upbnds, diff))
    lowarea <- bin$area - uparea # compute areas
    upper <- list(x = bin$x[!below], y = bin$y[!below],
                  bnds = upbnds, area = uparea,
                  expn = (bin$expn)*(uparea/bin$area),
                  n = bin$n - length(belowInds),
                  depth = bin$depth + 1)
    lower <- list(x = bin$x[below], y = bin$y[below],
                  bnds = lowbnds, area = lowarea,
                  expn = (bin$expn)*(lowarea/bin$area),
                  n = length(belowInds),
                  depth = bin$depth + 1) # construct
    list(lower, upper) # return
}

## something to plot/visualize bins
plotBinning <- function(bins, xlab = "x", ylab = "y",
                        main = "Bins", ...) {
    nbins <- length(binList) # number of bins
    xbnds <- sapply(binList, function(bn) bn$bnds$x)
    ybnds <- sapply(binList, function(bn) bn$bnds$y)
    plot(NA, type = "n", xlim = range(xbnds), ylim = range(ybnds),
         xlab = xlab, ylab = ylab, main = main)
    for (ii in 1:nbins) {
        rect(xbnds[1,ii], ybnds[1,ii], xbnds[2,ii], ybnds[2,ii])
        points(binList[[ii]]$x, binList[[ii]]$y, ...)
    }
}
