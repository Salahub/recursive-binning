##' @title Plot a binning using shaded rectangles
##' @description Use a binning and vector of fill
##' colours to visualize the sample space of pairwise data.
##' @details `plotBinning` plots each bin within a list of bins with
##' custom shading to communicate large residuals, the depth of bins,
##' or highlight particular bins
##' @param bins list of lists each with a named elements `x`, `y`, and
##' `bnds`, the last of which is a list having named elements `x` and
##' `y`
##' @param fill vector of values which can be interpreted as colours
##' of the same length as `bins`
##' @param jitter vector of length two where the first entry is true
##' if the `x` points should be jittered and the second is true if the
##' `y` margin should be jittered (both false by default)
##' @param add logical, should the plot of bins be added to the
##' current plot area?
##' @param xlab string, the label to be placed on the x axis
##' @param ylab string, the label to be placed on the y axis
##' @param border argument to be passed to `rect` internally giving
##' the border colour
##' @param ... optional additional arguments to be passed to `plot`,
##' `points`
##' @return A list of lists each with elements `x`, `y`, `bnds`,
##' `expn`, and `n`.
##' @examples
##' bin <- list(x = 1:10, y = sample(1:10),
##'             bnds = list(x = c(0, 10), y = c(0, 10)),
##'             expn = 10, n = 10, depth = 0)
##' bin2 <- halfSplit(bin, "x")
##' bin3 <- unlist(lapply(bin2, maxScoreSplit, scorer = chiScores),
##'                recursive = FALSE)
##' plotBinning(bin3)
##' @author Chris Salahub
plotBinning <- function(bins, fill, jitter = c(F, F), add = FALSE,
                        xlab = "x", ylab = "y", border = "black",
                        ...) {
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
             col = fill[ii], border = border)
        if (all(jitter)) {
            xa <- diff(bins[[ii]]$bnds$x)/2
            ya <- diff(bins[[ii]]$bnds$y)/2
            points(jitter(bins[[ii]]$x, amount = xa),
                   jitter(bins[[ii]]$y, amount = ya),
                   ...)
        } else if (jitter[1]) {
            xa <- diff(bins[[ii]]$bnds$x)/2
            points(jitter(bins[[ii]]$x, amount = xa),
                   bins[[ii]]$y, ...)
        } else if (jitter[2]) {
            ya <- diff(bins[[ii]]$bnds$y)/2
            points(bins[[ii]]$x,
                   jitter(bins[[ii]]$y, amount = ya),
                   ...)
        } else {
            points(bins[[ii]]$x, bins[[ii]]$y, ...)
        }
    }
}

##' Shadings
##' @title Generate fills encoding bin features
##' @description These functions all accept a list of bins and return
##' a vector of colours of the same length that encode some feature of
##' the bins.
##' @details Two functions are provided by default: one which
##' generates a fill based on bin depth and the other based on a
##' residual function applied to each bin.
##' @param bins list of bins to be visualized
##' @param colrng hue range to be passed to `colorRampPalette` to
##' generate the final hue scale
##' @param resFun function which returns a result with a name element
##' `residuals` that is a numeric vector of the same length as `bins`
##' @param maxRes numeric maximum value of the residuals to maintain
##' the correct origin, taken to be the maximum observed residual if
##' not provided
##' @param breaks numeric vector of breakpoints to control hues
##' @param nbr number of breakpoints for automatic breakpoint
##' generation if `breaks` is not provided
##' @return A vector of colours the same length as `bins`.
##' @examples
##' bin <- makeBin(x = 1:10, y = sample(1:10))
##' bin2 <- halfSplit(bin, "x")
##' bin3 <- unlist(lapply(bin2, maxScoreSplit,
##'                       scorer = chiScores, minExp = 2),
##'                recursive = FALSE)
##' plotBinning(bin3, fill = depthFill(bin3)) # all the same depth
##' plotBinning(bin3, fill = residualFill(bin3)) # diff resids
##' @author Chris Salahub
##' @describeIn shadings Fill by depth
depthFill <- function(bins, colrng = c("floralwhite", "firebrick")) {
    depths <- sapply(bins, function(bn) bn$depth)
    colorRampPalette(colrng)(max(depths))[depths]
}
##' @describeIn shadings Fill by residual values
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
