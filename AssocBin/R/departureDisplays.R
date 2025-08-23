##' depDisplay
##' @title Generate a departure display
##' @description This is a generic function which generates a
##' departure display to show the dependence between pairs of
##' variables for several common data structures.
##' @details `depDisplay` is a wrapper of the `plotBinning` function
##' with defaults set to be informative for most investigations.
##' @param x a `data.frame`, `DepSearch` object, or a vector
##' @param y an optional vector, only used if `x` is a vector
##' @param ... additional arguments to pass to plot
##' @param pair the pair of variables to display when `x` is a
##' `data.frame` or an `DepSearch`. If `x` is a `data.frame`, pair can
##' be specified in three ways: as a string with format "<y>:<z>",
##' as a character vector of length two, or as a numeric vector of
##' length two specifying the pair of variables to bin. If `x` is
##' an `DepSearch`, pair must be either a number or a string of the
##' format "<y>:<z>" specifying which binned pair of `x` to display.
##' @param quants list of two named vectors `x` and `y` providing the
##' quantiles to display on the corresponding axis in the case it is
##' a continuous variable. Defaults to the five number summary.
##' @param border string providing the colour of bin borders to draw,
##' NA suppresses borders
##' @return Invisibly returns the binning obtained and generates a
##' departure display of the pairwise dependence.
##' @examples
##' x <- rnorm(100)
##' y <- factor(abs(round(x*2)))
##' depDisplay(x, y)
##'
##' ## on the iris data
##' data(iris)
##' firstPair <- depDisplay(iris, pair = c(1,2))
##' ## another way
##' firstPair2 <- depDisplay(iris, pair = c("Sepal.Length", "Sepal.Width"))
##' ## a final way
##' firstPair2 <- depDisplay(iris, pair = "Sepal.Length:Sepal.Width")
##' @author Chris Salahub
depDisplay <- function(x, y, ..., pair, quants) {
    UseMethod("depDisplay")
}
##' @describeIn depDisplay Default depDisplay method
depDisplay.default <- function(x, y, ..., quants, border) {
    if (missing(quants)) {
        quants = list(x = seq(0, 1, length.out=5),
                      y = seq(0, 1, length.out=5))
    }
    if (missing(border)) {
        border = NA
    }
    colrng <- c("steelblue", "white", "firebrick")
    crits <- "depth >= 6 | n < 1 | expn <= 10 | stopped"
    stopFn <- function(bns) stopper(bns, crits)
    xcat <- class(x) %in% c("factor", "character", "logical")
    ycat <- class(y) %in% c("factor", "character", "logical")
    if (xcat & ycat) {
        binned <- catBinner(x, y)
        plotBinning(binned, factor = 0.9, border = border,
                    fill = importanceFill(binned, colrng = colrng,
                                          nbr = NA),
                    showXax = TRUE, showYax = TRUE,
                    ...)
    } else if (xcat) {
        y_ <- rank(y, ties.method="random")
        binned <- uniBinner(x, y_, on = "y",
                            stopper = stopFn,
                            splitter = uniRIntSplit)
        plotBinning(binned, factor = 0.9, border = border,
                    fill = importanceFill(binned, colrng = colrng,
                                          nbr = NA),
                    showXax = TRUE, showYax = FALSE,
                    ...)
        axis(side = 2, at = length(y)*quants$y, labels = FALSE)
        mtext(side = 2, at = length(y)*quants$y, line = 1,
              text = quantile(y, quants$y))
    } else if (ycat) {
        x_ <- rank(x, ties.method="random")
        binned <- uniBinner(x_, y, on = "x",
                            stopper = stopFn,
                            splitter = uniRIntSplit)
        plotBinning(binned, factor = 0.9, border = border,
                    fill = importanceFill(binned, colrng = colrng,
                                          nbr = NA),
                    showXax = FALSE, showYax = TRUE,
                    ...)
        axis(side = 1, at = length(x)*quants$x, labels = FALSE)
        mtext(side = 1, at = length(x)*quants$x, line = 1,
              text = quantile(x, quants$x))
    } else {
        x_ <- rank(x, ties.method="random")
        y_ <- rank(y, ties.method="random")
        binned <- binner(x_, y_, stopper = stopFn,
                         splitter = rIntSplit)
        plotBinning(binned, factor = 0.9, border = border,
                    fill = importanceFill(binned, colrng = colrng,
                                          nbr = NA),
                    showXax = FALSE, showYax = FALSE,
                    ...)
        axis(side = 1, at = length(x)*quants$x, labels = FALSE)
        axis(side = 2, at = length(y)*quants$y, labels = FALSE)
        mtext(side = 1, at = length(x)*quants$x, line = 1,
              text = quantile(x, quants$x))
        mtext(side = 2, at = length(y)*quants$y, line = 1,
              text = quantile(y, quants$y))
    }
    invisible(binned)
}
##' @describeIn depDisplay data.frame method for depDisplay
depDisplay.data.frame <- function(x, ..., pair, quants, border) {
    if (missing(quants)) {
        quants = list(x = seq(0, 1, length.out=5),
                      y = seq(0, 1, length.out=5))
    }
    if (missing(border)) {
        border = NA
    }
    if (missing(pair)) {
        warning("'pair' not provided, binning first two variables from 'names(x)'.")
        pair <- c(1, 2)
    } else if (!is.character(pair)) {
        if (is.numeric(pair)) {
            if (length(pair) < 2) {
                stop("'pair' must specify two variables to bin.")
            } else if (length(pair) != 2) {
                warning("More than two variables provided, binning only the first two.")
                pair <- pair[1:2]
            }
        } else {
            stop("'pair' must be a character or numeric vector.")
        }
    } else if (length(pair) == 1) {
        pair <- strsplit(pair, split = ":", fixed = TRUE)[[1]]
        if (length(pair) < 2) {
            stop("'pair' must have the format '<y>:<z>' if it is a string.")
        }
    } else if (length(pair) < 2) {
        stop("'pair' must specify two variables to bin.")
    } else if (length(pair) != 2) {
        warning("More than two variables provided, binning only the first two.")
        pair <- pair[1:2]
    }
    y <- x[[pair[1]]]
    z <- x[[pair[2]]]
    colrng <- c("steelblue", "white", "firebrick")
    crits <- "depth >= 6 | n < 1 | expn <= 10 | stopped"
    stopFn <- function(bns) stopper(bns, crits)
    ycat <- class(y) %in% c("factor", "character", "logical")
    zcat <- class(z) %in% c("factor", "character", "logical")
    if (ycat & zcat) { # in binner notation, x = y, y = z
        binned <- catBinner(y, z)
        plotBinning(binned, factor = 0.9, border = border,
                    fill = importanceFill(binned, colrng = colrng,
                                          nbr = NA),
                    showXax = TRUE, showYax = TRUE,
                    ...)
    } else if (ycat) {
        z <- rank(z, ties.method="random")
        binned <- uniBinner(y, z_, on = "y",
                            stopper = stopFn,
                            splitter = uniRIntSplit)
        plotBinning(binned, factor = 0.9, border = border,
                    fill = importanceFill(binned, colrng = colrng,
                                          nbr = NA),
                    showXax = TRUE, showYax = FALSE,
                    ...)
        axis(side = 2, at = length(z)*quants$y, labels = FALSE)
        mtext(side = 2, at = length(z)*quants$y, line = 1,
              text = quantile(z, quants$y))

    } else if (zcat) {
        y_ <- rank(y, ties.method="random")
        binned <- uniBinner(y_, z, on = "x",
                            stopper = stopFn,
                            splitter = uniRIntSplit)
        plotBinning(binned, factor = 0.9, border = border,
                    fill = importanceFill(binned, colrng = colrng,
                                          nbr = NA),
                    showXax = FALSE, showYax = TRUE,
                    ...)
        axis(side = 1, at = length(y)*quants$x, labels = FALSE)
        mtext(side = 1, at = length(y)*quants$x, line = 1,
              text = quantile(y, quants$x))
    } else {
        y_ <- rank(y, ties.method="random")
        z_ <- rank(z, ties.method="random")
        binned <- binner(y_, z_,
                         stopper = stopFn,
                         splitter = rIntSplit)
        plotBinning(binned, factor = 0.9, border = border,
                    fill = importanceFill(binned, colrng = colrng,
                                          nbr = NA),
                    showXax = FALSE, showYax = FALSE,
                    ...)
        axis(side = 1, at = length(y)*quants$x, labels = FALSE)
        axis(side = 2, at = length(z)*quants$x, labels = FALSE)
        mtext(side = 1, at = length(y)*quants$x, line = 1,
              text = quantile(y, quants$x))
        mtext(side = 2, at = length(z)*quants$y, line = 1,
              text = quantile(z, quants$y))
    }
    invisible(binned)
}
##' @describeIn depDisplay DepSearch method for depDisplay
depDisplay.DepSearch <- function(x, ..., pair, quants, border) {
    if (missing(quants)) {
        quants = list(x = seq(0, 1, length.out=5),
                      y = seq(0, 1, length.out=5))
    }
    if (missing(border)) {
        border = NA
    }
    if (missing(pair)) {
        pair <- 1
    } else if (!is.character(pair)) {
        if (is.numeric(pair)) {
            if (length(pair) > 1) {
                warning("More than one pair provided, displaying only the first.")
                pair <- pair[1]
            }
        } else {
            stop("'pair' must be a character or numeric.")
        }
    } else {
        clear <- grepl(":", pair)
        if (length(clear) > 1) {
            warning("More than one pair provided, displaying only the first.")
            clear <- clear[1]
        }
        if (!clear) {
            stop("'pair' is unclear, its format should be '<y>:<z>'.")
        }
        pair <- which(x$pairs == pair)
    }

    colrng <- c("steelblue", "white", "firebrick")
    depData <- tryCatch(get(x$datName),
                        error = function(e) NA)
    if (is.na(depData)) {
        print("x data not found, displaying without marginal quantiles.")
        plotBinning(x$binnings[[pair]], factor = 0.9, border = border,
                    fill = importanceFill(x$binnings[[pair]],
                                          colrng = colrng,
                                          nbr = NA),
                    ...)
    } else {
        pairVars <- strsplit(x$pairs[pair], ":")[[1]]
        y <- depData[[pairVars[1]]]
        z <- depData[[pairVars[2]]]
        ycat <- class(y) %in% c("factor", "character", "logical")
        zcat <- class(z) %in% c("factor", "character", "logical")
        if (ycat & zcat) {
            plotBinning(x$binnings[[pair]], factor = 0.9,
                        border = border,
                        fill = importanceFill(binned, colrng = colrng,
                                              nbr = NA),
                        showXax = TRUE, showYax = TRUE,
                        ...)
        } else if (ycat) {
            plotBinning(x$binnings[[pair]], factor = 0.9,
                        border = border,
                        fill = importanceFill(binned, colrng = colrng,
                                              nbr = NA),
                        showXax = TRUE, showYax = FALSE,
                        ...)
            axis(side = 2, at = length(z)*quants$y, labels = FALSE)
            mtext(side = 2, at = length(z)*quants$y, line = 1,
                  text = quantile(z, quants$y))
        } else if (zcat) {
            plotBinning(x$binnings[[pair]], factor = 0.9,
                        border = border,
                        fill = importanceFill(binned, colrng = colrng,
                                              nbr = NA),
                        showXax = FALSE, showYax = TRUE,
                        ...)
            axis(side = 1, at = length(y)*quants$x, labels = FALSE)
            mtext(side = 1, at = length(y)*quants$x, line = 1,
                  text = quantile(y, quants$x))
        } else {
            plotBinning(x$binnings[[pair]], factor = 0.9,
                        border = border,
                        fill = importanceFill(binned, colrng = colrng,
                                              nbr = NA),
                        showXax = FALSE, showYax = FALSE,
                        ...)
            axis(side = 1, at = length(y)*quants$x, labels = FALSE)
            axis(side = 2, at = length(z)*quants$y, labels = FALSE)
            mtext(side = 1, at = length(y)*quants$x, line = 1,
                  text = quantile(y, quants$x))
            mtext(side = 2, at = length(z)*quants$y, line = 1,
                  text = quantile(z, quants$y))
        }
    }
}
