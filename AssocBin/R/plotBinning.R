##' @title Plot a binning using shaded rectangles
##' @description Use a binning and vector of fill colours to
##' visualize the sample space of pairwise data.
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
##' @param factor number between 0 and 1, what factor should be
##' applied to jittering of categorical variables?
##' @param xlab string, the label to be placed on the x axis
##' @param ylab string, the label to be placed on the y axis
##' @param suppressLabs logical, should axis labels be suppressed or
##' displayed?
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
plotBinning <- function(bins, fill, add = FALSE, factor = 0.5,
                        xlab = "x", ylab = "y", suppressLabs = FALSE,
                        border = "black", ...) {
    if (missing(fill)) fill <- rep(NA, length(bins)) # custom fill
    nbins <- length(bins)
    xbnds <- t(sapply(bins, function(bn) bn$bnds$x))
    ybnds <- t(sapply(bins, function(bn) bn$bnds$y))
    xfac <- is.factor(bins[[1]]$x)
    yfac <- is.factor(bins[[1]]$y) # check bin 1 for factor status
    if (suppressLabs) {
        xaxt <- yaxt <- "n"
    } else xaxt <- yaxt <- "s"
    if (!add) { # create new plot area
        if (xfac & yfac) { # depends on what is a factor
            plot(NA, xlim = range(xbnds), ylim = range(ybnds),
                 xlab = xlab, ylab = ylab, xaxt = "n", yaxt = "n",
                 ...)
            unqx <- unique(xbnds)
            unqy <- unique(ybnds)
            xlocs <- sort((unqx[,1] + unqx[,2])/2)
            ylocs <- sort((unqy[,1] + unqy[,2])/2)
            if (!suppressLabs) {
                mtext(levels(bins[[1]]$x), at = xlocs, side = 1,
                      line = 1)
                mtext(levels(bins[[1]]$y), at = ylocs, side = 2,
                      line = 1)
            }
        } else if (xfac) {
            plot(NA, xlim = range(xbnds), ylim = range(ybnds),
                 xlab = xlab, ylab = ylab, xaxt = "n", yaxt = yaxt,
                 ...)
            unqx <- unique(xbnds)
            xlocs <- sort((unqx[,1] + unqx[,2])/2)
            if (!suppressLabs) {
                mtext(levels(bins[[1]]$x), at = xlocs, side = 1,
                      line = 1)
            }
        } else if (yfac) {
            plot(NA, xlim = range(xbnds), ylim = range(ybnds),
                 xlab = xlab, ylab = ylab, yaxt = "n", xaxt = xaxt,
                 ...)
            unqy <- unique(ybnds)
            ylocs <- sort((unqy[,1] + unqy[,2])/2)
            if (!suppressLabs) {
                mtext(levels(bins[[1]]$y), at = ylocs, side = 2,
                      line = 1)
            }
        } else {
            plot(NA, xlim = range(xbnds), ylim = range(ybnds),
                 xlab = xlab, ylab = ylab, xaxt = xaxt,
                 yaxt = yaxt, ...)
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
                           amount = xa,
                           factor = factor)
        } else pltx <- bins[[ii]]$x
        if (yfac) {
            ya <- diff(bins[[ii]]$bnds$y)/2
            plty <- jitter(rep((bins[[ii]]$bnds$y[1] +
                                bins[[ii]]$bnds$y[2])/2,
                               bins[[ii]]$n),
                           amount = ya,
                           factor = factor)
        } else plty <- bins[[ii]]$y
        points(pltx, plty, ...)
    }
}

##' Shadings
##' @title Generate fills encoding bin features
##' @description These functions all accept a list of bins and return
##' a vector of colours of the same length that encode some feature of
##' the bins. standardizedChiFill is a special case which adjusts the
##' residuals obtained by the binChi function by the variance of each
##' bin to obtain a better normal approximation and more even
##' gradient.
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
##' @param breaks numeric vector of breakpoints to control hues,
##' defaults to breakpoints that indicate Pearson residuals outside
##' the asymptotic 95 percent confidence interval around zero under
##' the null
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
depthFill <- function(bins, colrng = c("white", "firebrick")) {
    depths <- sapply(bins, function(bn) bn$depth)
    colorRampPalette(colrng)(max(depths))[depths]
}
##' @describeIn shadings Fill by residual values
residualFill <- function(bins, resFun = binChi, maxRes,
                         colrng = c("steelblue", "white",
                                    "firebrick"),
                         breaks = NA, nbr = NA) {
    residuals <- resFun(bins)$residuals # get residuals
    if (missing(maxRes)) maxRes <- 1.01*max(abs(residuals))
    if (is.na(breaks)) {
        if (is.na(nbr)) { # default: sig. residuals
            breaks <- sort(c(-maxRes, -1.96, 1.96, maxRes))
        } else {
            breaks <- seq(-maxRes, maxRes, length.out = nbr)
        }
    }
    residCols <- cut(residuals, breaks) # distribute colors
    colorRampPalette(colrng)(length(breaks)-1)[as.numeric(residCols)]
}
##' @describeIn shadings Fill by variance-adjusted chi residuals
standardizedChiFill <- function(bins, nbr = NA, breaks = NA,
                                colrng = c("steelblue", "white",
                                           "firebrick")) {
    wids <- sapply(bins, function(x) diff(x$bnds$x))
    hgts <- sapply(bins, function(x) diff(x$bnds$y))
    obs <- sapply(bins, function(x) x$n)
    N <- sum(obs)
    expn <- wids*hgts/N
    widAdj <- 1 - wids/N # avoid zero denom
    widAdj[widAdj < 2*.Machine$double.eps] <- 1/N
    hgtAdj <- 1 - hgts/N
    hgtAdj[hgtAdj < 2*.Machine$double.eps] <- 1/N
    denom <- N/(N-1)*widAdj*hgtAdj*expn
    stRes <- (obs - expn)^2/sqrt(denom)*sign(obs - expn)
    maxRes <- 1.01*max(abs(stRes))
    if (is.na(breaks)) {
        if (is.na(nbr)) { # default: 16 splits over -4, 4
            nbr <- 16
        }
        if (maxRes > 4) {
            breaks <- c(-maxRes,
                        seq(-4, 4, length.out = nbr),
                        maxRes)
        } else { # arbitrary larger point to keep them constant
            breaks <- c(-10, seq(-4, 4, length.out = nbr), 10)
        }
    } else {
        breaks <- sort(breaks)
    }
    residCols <- cut(stRes, breaks) # distribute colors
    colorRampPalette(colrng)(length(breaks)-1)[as.numeric(residCols)]
}
