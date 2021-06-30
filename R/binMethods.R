## define the formal class "Bin" to support later processing,
## store the bounds, x coordinates inside, y coordinates inside,
## number of points inside, and the depth it was created
setClass("Bin", representation(x = "numeric", y = "numeric",
                               bounds = "list", size = "numeric",
                               depth = "numeric"))
## the constructor to create a new bin
Bin <- function(x, y, bounds, size, depth) {
    stopifnot(length(x) == length(y))
    if (missing(bounds)) bounds = list(x = range(x), y = range(y))
    new("Bin", x = x, y = y, bounds = bounds, size = length(x),
        depth = depth)
}
## a lighter constructor without checks
quickBin <- function(x, y, bounds, size, depth) {
    new("Bin", x = x, y = y, bounds = bounds, size = size,
        depth = depth)
}

## calculate the bin area using the bounds
setGeneric("BinArea",
           function(Bin) standardGeneric("BinArea"))
setMethod("BinArea", signature(Bin = "Bin"),
          function(Bin) prod(diff(Bin@bounds$x),
                             diff(Bin@bounds$y)))

## given a list of bins, calculate the observed and expected counts
getObsExp <- function(bins) {
    allBounds <- lapply(bins, attr, "bounds")
    xrng <- range(sapply(allBounds, function(bd) bd$x))
    yrng <- range(sapply(allBounds, function(bd) bd$y)) # to compute density
    ob <- sapply(bins, attr, which = "size")
    n <- sum(ob)
    ex <- n*(sapply(bins, BinArea)/(diff(xrng)*diff(yrng))) # n times area over density
    list(ob = ob, ex = ex, n = n)
}

## wrapper functions to compute different statistics over the bins
## a chi^2 measure
binChi <- function(bins, agg = sum) {
    feats <- getObsExp(bins) # observed and expected
    resids <- with(feats, (ob-ex)^2/ex)
    signs <- with(feats, sign(ob-ex)) # signs of points
    list(residuals = signs*sqrt(resids), stat = agg(resids))
}

## a mutual information measure
binMI <- function(bins, agg = sum) {
    feats <- getObsExp(bins)
    resids <- with(feats, log(ob/ex))
    resids[feats$ob == 0] <- 0 # remove NAs
    probs <- with(feats, ob/n)
    list(residuals = resids, stat = agg(resids*probs))
}

## simply sum the absolute values of the counts
binCount <- function(bins, agg = sum) {
    feats <- getObsExp(bins)
    resids <- with(feats, abs(ob - ex))
    signs <- with(feats, sign(ob-ex))
    list(residuals = signs*resids, stat = agg(resids))
}



## a function to plot a single bin
setGeneric("BinPlot", function(Bin, bg, ...) standardGeneric("BinPlot"))
setMethod("BinPlot", signature(Bin = "Bin"),
          function(Bin, bg, ...) {
              if (is.null(dev.list())) { # check if plot already active
                  plot(NA, xlim = range(Bin@x), ylim = range(Bin@y), ...)
              }
              polygon(x = rep(Bin@bounds$x, each = 2), # add polygons to plot
                      y = c(Bin@bounds$y, rev(Bin@bounds$y)),
                      border = "gray50", col = bg)
              points(Bin@x, Bin@y, ...) # add points to plot
              })
## a wrapper to plot a list of bins using BinPlot
plotBinning <- function(bins, fill, ...) {
    if (missing(fill)) fill <- rep(NA, length(bins)) # custom fill option
    allBounds <- lapply(bins, attr, "bounds")
    xs <- sapply(allBounds, function(bd) bd$x)
    ys <- sapply(allBounds, function(bd) bd$y)
    plot(NA, xlim = range(xs), ylim = range(ys), xlab = "x", ylab = "y",
         ...)
    for (ii in seq_along(bins)) BinPlot(bins[[ii]], bg = fill[ii],...)
    points(unlist(sapply(bins, attr, "x")), unlist(sapply(bins, attr, "y")),
           ...) # can be disabled with pch = ""
}

## colour bins by their depth in the algorithm
depthFill <- function(bins, colrng = c("floralwhite", "firebrick")) {
    depths <- sapply(bins, attr, "depth")
    colorRampPalette(colrng)(max(depths))[depths]
}
## colour bins by their residual value
## the colour cuts have not been determined statistically, this would
## be a nice extension
residualFill <- function(bins, resFun = binChi, maxRes,
                         colrng = c("steelblue", "floralwhite", "firebrick"),
                         ncol = 50) {
    residuals <- resFun(bins)$residuals # get residuals
    if (missing(maxRes)) maxRes <- 1.01*max(abs(residuals))
    coldivs <- seq(-maxRes, maxRes, length.out = ncol)
    residCols <- cut(residuals, coldivs) # add the option for custom colour cuts?
    colorRampPalette(colrng)(ncol)[as.numeric(residCols)]
}
