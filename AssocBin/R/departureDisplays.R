##' Departure display
##' @title Generate a departure display
##' @description This is a generic function which generates a
##' departure display to show the dependence between pairs of
##' variables for several common data structures.
##' @details `depDisplay` is a wrapper of the `plotBinning` function
##' with defaults set to be informative for most investigations.
##' @param x a `data.frame`, `inDep` object, or a vector
##' @param y an optional vector, only used if `x` is a vector
##' @param pair the pair of variables to display when `x` is a
##' `data.frame` or an `inDep`. If `x` is a `data.frame`, pair can
##' be specified in three ways: as a string with format "<y>:<z>",
##' as a character vector of length two, or as a numeric vector of
##' length two specifying the pair of variables to bin. If `x` is
##' an `inDep`, pair must be either a number or a string of the
##' format "<y>:<z>" specifying which binned pair of `x` to display.
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
depDisplay <- function(x, y, ..., pair) {
    UseMethod("depDisplay")
}

##' @describeIn Departuredisplay Default depDisplay method
depDisplay.default <- function(x, y, ...) {
    colrng <- c("steelblue", "white", "firebrick")
    crits <- "depth >= 6 | n < 1 | expn <= 10 | stopped"
    stopFn <- function(bns) stopper(bns, crits)
    xcat <- class(x) %in% c("factor", "character", "logical")
    ycat <- class(y) %in% c("factor", "character", "logical")
    if (xcat & ycat) {
        binned <- catBinner(x, y)
    } else if (xcat) {
        binned <- uniBinner(x, y, on = "y",
                            stopper = stopFn,
                            splitter = uniRIntSplit)
    } else if (ycat) {
        binned <- uniBinner(x, y, on = "x",
                            stopper = stopFn,
                            splitter = uniRIntSplit)
    } else {
        binned <- binner(x, y, rIntSplit)
    }
    plotBinning(binned, factor = 0.9, border = NA,
                fill = importanceFill(binned, colrng = colrng,
                                      nbr = NA),
                ...)
    invisible(binned)
}

##' @describeIn Departuredisplay data.frame method for depDisplay
depDisplay.data.frame <- function(x, ..., pair) {
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

    y <- x[pair[1]]
    z <- x[pair[2]]
    colrng <- c("steelblue", "white", "firebrick")
    crits <- "depth >= 6 | n < 1 | expn <= 10 | stopped"
    stopFn <- function(bns) stopper(bns, crits)
    ycat <- class(y) %in% c("factor", "character", "logical")
    zcat <- class(z) %in% c("factor", "character", "logical")
    if (ycat & zcat) {
        binned <- catBinner(y, z)
    } else if (ycat) { # in binner notation, x = y, y = z
        binned <- uniBinner(y, z, on = "y",
                            stopper = stopFn,
                            splitter = uniRIntSplit)
    } else if (zcat) {
        binned <- uniBinner(y, z, on = "x",
                            stopper = stopFn,
                            splitter = uniRIntSplit)
    } else {
        binned <- binner(y, z, rIntSplit)
    }
    plotBinning(binned, factor = 0.9, border = NA,
                fill = importanceFill(binned, colrng = colrng,
                                      nbr = NA),
                ...)
    invisible(binned)
}

##' @describeIn Departuredisplay inDep method for depDisplay
depDisplay.inDep <- function(x, ..., pair) {
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
    plotBinning(x$binnings[[pair]], factor = 0.9, border = NA,
                fill = importanceFill(x$binnings[[pair]], colrng = colrng,
                                      nbr = NA),
                ...)
}
