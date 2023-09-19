
## a wrapper to plot a list of bins using BinPlot
plotBinning <- function(bins, fill, add = FALSE, xlab = "x",
                        ylab = "y", ...) {
    if (missing(fill)) fill <- rep(NA, length(bins)) # custom fill option
    nbins <- length(bins)
    xbnds <- sapply(bins, function(bn) bn$bnds$x)
    ybnds <- sapply(bins, function(bn) bn$bnds$y)
    if (!add) {
        plot(NA, xlim = range(xbnds), ylim = range(ybnds), xlab = xlab,
             ylab = ylab, ...)
    }
    for (ii in seq_along(bins)) {
        rect(xbnds[1,ii], ybnds[1,ii], xbnds[2,ii], ybnds[2,ii],
             col = fill[ii])
        points(bins[[ii]]$x, bins[[ii]]$y, ...) # disable with pch = ""
    }
}

## colour bins by their depth in the algorithm
depthFill <- function(bins, colrng = c("floralwhite", "firebrick")) {
    depths <- sapply(bins, function(bn) bn$depth)
    colorRampPalette(colrng)(max(depths))[depths]
}
## colour bins by their residual value
## the colour cuts have not been determined statistically, this would
## be a nice extension
residualFill <- function(bins, resFun = binChi, maxRes,
                         colrng = c("steelblue", "floralwhite",
                                    "firebrick"),
                         breaks = NA, nbr = 50) {
    residuals <- resFun(bins)$residuals # get residuals
    if (missing(maxRes)) maxRes <- 1.01*max(abs(residuals))
    if (is.na(breaks)) {
        breaks <- seq(-maxRes, maxRes, length.out = nbr)
    }
    residCols <- cut(residuals, breaks) # distribute colors
    colorRampPalette(colrng)(nbr)[as.numeric(residCols)]
}
