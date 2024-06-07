##' @title Plot a binning using shaded rectangles
##' @description Use a binning and vector of fill
##' colours to visualize the sample space of pairwise data.
##' @details `plotBinning` plots each bin within a list of bins with
##' custom shading to communicate large residuals, the depth of bins,
##' or highlight particular bins. It automatically jitters points
##' within categorical levels to avoid overplotting.
##' @param bins list of lists each with a named elements `x`, `y`, and
##' `bnds`, the last of which is a list having named elements `x` and
##' `y`
##' @param fill vector of values which can be interpreted as colours
##' of the same length as `bins`
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
plotBinning <- function(bins, fill, add = FALSE, xlab = "x",
                        ylab = "y", border = "black", ...) {
    if (missing(fill)) fill <- rep(NA, length(bins)) # custom fill
    nbins <- length(bins)
    xbnds <- t(sapply(bins, function(bn) bn$bnds$x))
    ybnds <- t(sapply(bins, function(bn) bn$bnds$y))
    xfac <- is.factor(bins[[1]]$x)
    yfac <- is.factor(bins[[1]]$y) # check bin 1 for factor status
    if (!add) { # create new plot area
        if (xfac & yfac) { # depends on what is a factor
            plot(NA, xlim = range(xbnds), ylim = range(ybnds),
                 xlab = xlab, ylab = ylab, xaxt = "n", yaxt = "n",
                 ...)
            unqx <- unique(xbnds)
            unqy <- unique(ybnds)
            xlocs <- sort((unqx[,1] + unqx[,2])/2)
            ylocs <- sort((unqy[,1] + unqy[,2])/2)
            mtext(levels(bins[[1]]$x), at = xlocs, side = 1, line = 1)
            mtext(levels(bins[[1]]$y), at = ylocs, side = 2, line = 1)
        } else if (xfac) {
            plot(NA, xlim = range(xbnds), ylim = range(ybnds),
                 xlab = xlab, ylab = ylab, xaxt = "n", ...)
            unqx <- unique(xbnds)
            xlocs <- sort((unqx[,1] + unqx[,2])/2)
            mtext(levels(bins[[1]]$x), at = xlocs, side = 1, line = 1)
        } else if (yfac) {
            plot(NA, xlim = range(xbnds), ylim = range(ybnds),
                 xlab = xlab, ylab = ylab, yaxt = "n", ...)
            unqy <- unique(ybnds)
            ylocs <- sort((unqy[,1] + unqy[,2])/2)
            mtext(levels(bins[[1]]$y), at = ylocs, side = 2, line = 1)
        } else {
            plot(NA, xlim = range(xbnds), ylim = range(ybnds),
                 xlab = xlab, ylab = ylab, ...)
        }
    } # add bins and points
    for (ii in seq_along(bins)) {
        rect(xbnds[ii,1], ybnds[ii,1], xbnds[ii,2], ybnds[ii,2],
             col = fill[ii], border = border)
        if (xfac) {
            xa <- diff(bins[[ii]]$bnds$x)/2
            pltx <- jitter(rep((bins[[ii]]$bnds$x[1] +
                                bins[[ii]]$bnds$x[2])/2,
                               bins[[ii]]$n),
                           amount = xa)
        } else pltx <- bins[[ii]]$x
        if (yfac) {
            ya <- diff(bins[[ii]]$bnds$y)/2
            plty <- jitter(rep((bins[[ii]]$bnds$y[1] +
                                bins[[ii]]$bnds$y[2])/2,
                               bins[[ii]]$n),
                           amount = ya)
        } else plty <- bins[[ii]]$y
        points(pltx, plty, ...)
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
