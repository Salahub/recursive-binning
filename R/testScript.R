## custom plotting function with narrow margins
narrowPlot <- function(xgrid, ygrid, main = "", xlab = "", ylab = "",
                       xticks = xgrid, yticks = ygrid,
                       mars = c(2.1, 2.1, 1.1, 1.1),
                       xlim = range(xgrid), ylim = range(ygrid),
                       addGrid = TRUE, ...) {
    par(mar = mars) # set narrow margins
    plot(NA, ylim = ylim, xlim = xlim, xaxt = 'n', xlab = "",
         yaxt = 'n', ylab = "", main = "", ...)
    ## add labels
    mtext(main, side = 3, line = 0, cex = 0.8) # main
    mtext(ylab, side = 2, line = 1, cex = 0.8) # ylab
    mtext(xlab, side = 1, line = 1, padj = 0, cex = 0.8) # xlab
    ## add grid lines
    if (addGrid) {
        abline(h = ygrid, v = xgrid, lty = 1,
               col = adjustcolor("gray", alpha.f = 0.4))
    }
    ## and ticks
    mtext(side = 1, at = xgrid, text = "|", line = 0, cex = 0.5,
          padj = -2)
    mtext(text = xticks, at = xgrid, side = 1, cex = 0.8)
    mtext(side = 2, at = ygrid, text = "|", line = 0, cex = 0.5,
          padj = 1)
    mtext(text = yticks, at = ygrid, side = 2, cex = 0.8)
}

## consider the general form of a binner for the data
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

## bin: list with x, y, bnds, expn, n, depth

## some scoring functions to be maximized
chiScores <- function(vals, expn, minExp = 0) {
    diffs <- diff(vals)
    n <- length(vals) - 2
    total <- c(diffs[1]-1, cumsum(diffs))
    h1 <- total[1:(n+1)] # length below
    h2 <- total[n+2] - h1 # length above
    d <- n/total[n+2] # density
    i <- 0:n # number below split i
    ni <- n - i # number above i
    scr <- (i - d*h1)^2/(h1*d) + (ni - d*h2)^2/(h2*d)
    scr[is.na(scr)] <- 0
    scr[pmin(expn*h1/total[n+2],
             expn*h2/total[n+2]) < minExp] <- 0 # minimum size limit
    scr
}
## mutual informations
miScores <- function(vals, expn, minExp = 0) {
    diffs <- diff(vals)
    n <- length(vals) - 2
    total <- c(diffs[1] - 1, cumsum(diffs))
    h1 <- total[1:(n+1)] # length below
    h2 <- total[n+2] - h1 # length above
    d <- n/total[n+2] # density
    i <- 0:n # number below point i
    ni <- n - i # number above i
    below <- (i/n)*log(i/(d*h1))
    above <- (ni/n)*log(ni/(d*h2)) # split expectation
    below[1] <- 0
    above[n+1] <- 0 # handle known zeros
    scr <- below + above
    scr[pmin(expn*h1/total[n+2],
             expn*h2/total[n+2]) < minExp] <- 0 # minimum size limit
    scr
}
## random scoring for random splits
randScores <- function(vals, expn, minExp = 0) {
    diffs <- diff(vals)
    n <- length(vals) - 2
    scores <- runif(length(diffs))
    ## compute the expected counts at each split
    total <- c(diffs[1]-1, cumsum(diffs))
    h1 <- total[1:(n+1)] # length below
    h2 <- total[n+2] - h1 # length above
    ## if the difference is one, splitting creates bin with area 0
    scores[1] <- min(diffs[1]-1, scores[1])
    ## difference of zero here does the same
    scores[length(scores)] <- min(diffs[length(diffs)],
                                  scores[length(scores)])
    scores[pmin(expn*h1/total[n+2],
                expn*h2/total[n+2]) < minExp] <- 0 # minimum size limit
    scores
}

## helper functions to split margins at bound based on indices
splitX <- function(bin, bd, above, below) {
    belowfac <- (bd - bin$bnds$x[1])/diff(bin$bnds$x)
    abovefac <- (bin$bnds$x[2] - bd)/diff(bin$bnds$x)
    list(list(x = bin$x[below], y = bin$y[below],
              bnds = list(x = c(bin$bnds$x[1], bd),
                          y = bin$bnds$y),
              expn = bin$expn*belowfac,
              n = bin$n-length(above), depth = bin$depth + 1),
         list(x = bin$x[above], y = bin$y[above],
              bnds = list(x = c(bd, bin$bnds$x[2]),
                          y = bin$bnds$y),
              expn = bin$expn*abovefac,
              n = length(above), depth = bin$depth + 1))
}
splitY <- function(bin, bd, above, below) {
    belowfac <- (bd - bin$bnds$y[1])/diff(bin$bnds$y)
    abovefac <- (bin$bnds$y[2] - bd)/diff(bin$bnds$y)
    list(list(x = bin$x[below], y = bin$y[below],
              bnds = list(x = bin$bnds$x,
                          y = c(bin$bnds$y[1], bd)),
              expn = bin$expn*belowfac,
              n = bin$n-length(above), depth = bin$depth + 1),
         list(x = bin$x[above], y = bin$y[above],
              bnds = list(x = bin$bnds$x,
                          y = c(bd, bin$bnds$y[2])),
              expn = bin$expn*abovefac,
              n = length(above), depth = bin$depth + 1))
}

## halve a bin
halfSplit <- function(bin, margin = "x") {
    if (margin == "x") {
        xsort <- order(bin$x)
        hind <- floor(bin$n/2) # middle index
        newbnd <- bin$x[xsort][hind] # middle value
        above <- xsort[(hind+1):(bin$n)] # points above
        below <- xsort[1:hind] # and below
        splitX(bin, bd = newbnd, above = above, below = below)
    } else if (margin == "y") {
        ysort <- order(bin$y)
        hind <- floor(bin$n/2) # middle index
        newbnd <- bin$y[ysort][hind] # middle value
        above <- ysort[(hind+1):(bin$n)] # points above
        below <- ysort[1:hind]
        splitY(bin, bd = newbnd, above = above, below = below)
    } else stop("Margin must be one of x or y")
}

## choose maximum values subject to constraints
sizeLimMax <- function(scores, lim = 10) {
    which.max(scores[(lim+1):(length(scores)-lim)]) + lim
}

## splitter maximizing a score function
maxScoreSplit <- function(bin, scorer = diff, pickMax = which.max,
                          ...) {
  xsort <- order(bin$x)
  ysort <- order(bin$y) # get marginal ordering
  xscore <- scorer(c(bin$bnds$x[1], bin$x[xsort], bin$bnds$x[2]),
                   expn = bin$expn, ...)
  yscore <- scorer(c(bin$bnds$y[1], bin$y[ysort], bin$bnds$y[2]),
                   expn = bin$expn, ...)
  xmax <- pickMax(xscore)
  ymax <- pickMax(yscore) # the score values
  xallEq <- all(abs(xscore - xscore[1]) < sqrt(.Machine$double.eps))
  yallEq <- all(abs(yscore - yscore[1]) < sqrt(.Machine$double.eps))
  if (xallEq & yallEq) { # no variation, halve a random margin
      hind <- floor(bin$n/2)
      u <- runif(1)
      if (u < 0.5) {
          newbnd <- bin$y[ysort][hind]
          above <- ysort[(hind+1):(bin$n)]
          below <- ysort[1:hind]
          splitY(bin, bd = newbnd, above = above, below = below)
      } else {
          newbnd <- bin$x[xsort][hind]
          above <- xsort[(hind+1):(bin$n)]
          below <- xsort[1:hind]
          splitX(bin, bd = newbnd, above = above, below = below)
      }
  } else if (xscore[xmax] >= yscore[ymax]) { # ties go to x
      xsplts <- bin$x[xsort]
      newbnd <- c(xsplts[1]-1, xsplts)[xmax] # new boundary
      below <- xsort[seq_len(xmax-1)] # get indices of points below
      above <- if (xmax == bin$n+1) integer(0) else xsort[xmax:bin$n]
      splitX(bin, bd = newbnd, above = above, below = below)
  } else { # do the same on y
      ysplts <- bin$y[ysort]
      newbnd <- c(ysplts[1]-1, ysplts)[ymax]
      below <- ysort[seq_len(ymax-1)]
      above <- if (ymax == bin$n+1) integer(0) else ysort[ymax:bin$n]
      splitY(bin, bd = newbnd, above = above, below = below)
  }
}

## a univariate version
uniMaxScoreSplit <- function(bin, scorer = diff,
                             pickMax = which.max,
                             ...) {
  xsort <- order(bin$x)
  xscore <- scorer(c(bin$bnds$x[1], bin$x[xsort], bin$bnds$x[2]),
                   expn = bin$expn, ...)
  xmax <- pickMax(xscore)
  xsplts <- bin$x[xsort]
  newbnd <- c(xsplts[1]-1, xsplts)[xmax]  # new bin boundary
  below <- xsort[seq_len(xmax-1)]
  above <- if (xmax == bin$n+1) integer(0) else xsort[xmax:bin$n]
  splitX(bin, bd = newbnd, above = above, below = below)
}

##' make criteria by capturing expressions
makeCriteria <- function(...) {
    cl <- match.call() # capturing inputs
    crits <- as.list(cl) # change to a list
    ## remove self reference, collapse into single OR
    paste(sapply(crits[-1], deparse), collapse = " | ")
}

##' take critera, apply to all bins in list
stopper <- function(binList, criteria) {
    sapply(binList,
           function(b) eval(parse(text = criteria), envir = b))
}

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
        points(bins[[ii]]$x, bins[[ii]]$y, ...) # diasble with pch = ""
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

## wrapper functions to compute different statistics over the bins
## a chi^2 measure
binChi <- function(bins, agg = sum) {
    obs <- sapply(bins, function(bn) bn$n)
    ex <- sapply(bins, function(bn) bn$expn)
    resids <- (obs - ex)^2/ex
    signs <- sign(obs - ex) # signs of residuals
    list(residuals = signs*sqrt(resids), stat = agg(resids))
}

## a mutual information measure
binMI <- function(bins, agg = sum) {
    obs <- sapply(bins, function(bin) bin$n)
    ex <- sapply(bins, function(bin) bin$expn)
    n <- sum(obs)
    resids <- log(obs/ex)
    resids[obs == 0] <- 0
    probs <- obs/n
    list(residuals = resids, stat = agg(resids*probs))
}

## absolute values of the counts
binAbsDif <- function(bins, agg = sum) {
    obs <- sapply(bins, function(bin) bin$n)
    ex <- sapply(bins, function(bin) bin$expn)
    resids <- abs(obs - ex)
    signs <- sign(obs - ex)
    list(residuals = signs*resids, stat = agg(resids))
}

plotShading <- function(bins, resFun, brks = c(-1, 1, by = 0.1),
                        pal = colorRampPalette(
                            c("steelblue", "white",
                              "firebrick"))(length(brks)),
                        xlab = "x", ylab = "y") {
    nbins <- length(bins)
    xbnds <- sapply(bins, function(bn) bn$bnds$x)
    ybnds <- sapply(bins, function(bn) bn$bnds$y)
    res <- cut(sapply(bins, resFun), brks) # residual scores binned
    plot(NA, type = "n", xlim = range(xbnds), ylim = range(ybnds),
         xlab = xlab, ylab = ylab, main = main)
    for (ii in 1:nbins) {
        rect(xbnds[1,ii], ybnds[1,ii], xbnds[2,ii], ybnds[2,ii],
             col = pal[unclass(res)])
    }
}

## TRY THESE OUT
sizeLim10 <- function(scores) sizeLimMax(scores, lim = 10)
criteria <- makeCriteria(expn <= 10, n <= 20, depth >= 10)
## monotonic associations
test <- binner(1:1000, 1:1000,
               stopper = function(bns) stopper(bns, criteria),
               splitter = function(bin) {
                   maxScoreSplit(bin, chiScores)
               })



## SIMPLE EXAMPLES ###################################################
## generate some random data
set.seed(16062021)
randx <- sample(1:1e3)
randy <- sample(1:1e3) # random ranks
## plots the progression of the algorithm for several different depths
dep <- c(1, 2, 10)
## R's lexical scoping allows for dynamic criteria construction: d is
## not defined, but will be within a later loop
criteria <- makeCriteria(expn <= 10, n <= 1, depth >= d)
## define a stopper based on this
stopFn <- function(bns) stopper(bns, criteria)
## and some splitters for different score functions
chiSplit <- function(bn) maxScoreSplit(bn, chiScores, minExp = 5)
miSplit <- function(bn) maxScoreSplit(bn, miScores, minExp = 5)
rndSplit <- function(bn) maxScoreSplit(bn, randScores, minExp = 5)

## random data plot (Fig 4.1(a))
png("randomData.png", width = 3, height = 3, units = "in", res = 480)
narrowPlot(xgrid = seq(0, 1000, by = 250),
           ygrid = seq(0, 1000, by = 250),
           addGrid = FALSE, xlab = "x", ylab = "y")
points(randx, randy, pch = 19, cex = 0.5,
       col = adjustcolor("gray50", 0.5))
dev.off()

## allocate storage
randBin.chi <- vector(mode = "list", length = length(dep))
## split random data using different rules and depths
for (ii in seq_along(dep)) {
    d <- dep[ii]
    randBin.chi[[ii]] <- binner(randx, randy, stopper = stopFn,
                                splitter = chiSplit)
}

## the early progression of splitting (Figs 4.2(a), 4.3(a))
png("randomSplit1.png", width = 3, height = 3, units = "in", res = 480)
narrowPlot(xgrid = seq(0, 1000, by = 250),
           ygrid = seq(0, 1000, by = 250),
           addGrid = FALSE, xlab = "x", ylab = "y")
plotBinning(randBin.chi[[1]], pch = 19, cex = 0.5, add = TRUE,
            col = adjustcolor("gray50", 0.5))
dev.off()
png("randomSplit2.png", width = 3, height = 3, units = "in", res = 480)
narrowPlot(xgrid = seq(0, 1000, by = 250),
           ygrid = seq(0, 1000, by = 250),
           addGrid = FALSE, xlab = "x", ylab = "y")
plotBinning(randBin.chi[[2]], pch = 19, cex = 0.5, add = TRUE,
            col = adjustcolor("gray50", 0.5))
dev.off()

## plot bins coloured by depth (Fig 4.4(a))
png("randomSplitDepth.png", width = 3, height = 3, units = "in",
    res = 480)
narrowPlot(xgrid = seq(0, 1000, by = 250),
           ygrid = seq(0, 1000, by = 250),
           addGrid = FALSE, xlab = "x", ylab = "y")
plotBinning(randBin.chi[[3]], pch = 19, cex = 0.5, add = TRUE,
            fill = depthFill(randBin.chi[[3]]),
            col = adjustcolor("black", 0.5))
dev.off()

## data with perfect rank agreement
linex <- 1:1e3
liney <- 1:1e3

## allocate storage
lineBin.chi <- vector(mode = "list", length = length(dep))
## split random data using different rules and depths
for (ii in seq_along(dep)) {
    d <- dep[ii]
    lineBin.chi[[ii]] <- binner(linex, liney, stopper = stopFn,
                                splitter = chiSplit)
}

## plot the line data (Fig 4.1(b))
png("lineData.png", width = 3, height = 3, units = "in", res = 480)
narrowPlot(xgrid = seq(0, 1000, by = 250),
           ygrid = seq(0, 1000, by = 250),
           addGrid = FALSE, xlab = "x", ylab = "y")
points(linex, liney, pch = 19, cex = 0.5,
       col = adjustcolor("gray50", 0.5))
dev.off()

## the early progression of splitting (Figs 4.2(b), 4.3(b))
png("lineSplit1.png", width = 3, height = 3, units = "in", res = 480)
narrowPlot(xgrid = seq(0, 1000, by = 250),
           ygrid = seq(0, 1000, by = 250),
           addGrid = FALSE, xlab = "x", ylab = "y")
plotBinning(lineBin.chi[[1]], pch = 19, cex = 0.5, add = TRUE,
            col = adjustcolor("gray50", 0.5))
dev.off()
png("lineSplit2.png", width = 3, height = 3, units = "in", res = 480)
narrowPlot(xgrid = seq(0, 1000, by = 250),
           ygrid = seq(0, 1000, by = 250),
           addGrid = FALSE, xlab = "x", ylab = "y")
plotBinning(lineBin.chi[[2]], pch = 19, cex = 0.5, add = TRUE,
            col = adjustcolor("gray50", 0.5))
dev.off()

## plot bins coloured by depth (Fig 4.4(b))
png("lineSplitDepth.png", width = 3, height = 3, units = "in",
    res = 480)
narrowPlot(xgrid = seq(0, 1000, by = 250),
           ygrid = seq(0, 1000, by = 250),
           addGrid = FALSE, xlab = "x", ylab = "y")
plotBinning(lineBin.chi[[3]], pch = 19, cex = 0.5, add = TRUE,
            fill = depthFill(lineBin.chi[[3]]),
            col = adjustcolor("black", 0.5))
dev.off()

## compute pearson residuals for both and take the maximum magnitude
## to ensure the shading is consistent between them
maxRes <- max(abs(c(binChi(lineBin.chi[[3]])$residuals,
                    binChi(randBin.chi[[3]])$residuals)))
## plot the perason residuals for both (Fig 4.5)
png("randomSplitResid.png", width = 3, height = 3, units = "in",
    res = 480)
narrowPlot(xgrid = seq(0, 1000, by = 250),
           ygrid = seq(0, 1000, by = 250),
           addGrid = FALSE, xlab = "x", ylab = "y")
plotBinning(randBin.chi[[3]], pch = 19, cex = 0.5, add = TRUE,
            fill = residualFill(randBin.chi[[3]], maxRes = maxRes),
            col = adjustcolor("black", 0.5))
dev.off()
png("lineSplitResid.png", width = 3, height = 3, units = "in",
    res = 480)
narrowPlot(xgrid = seq(0, 1000, by = 250),
           ygrid = seq(0, 1000, by = 250),
           addGrid = FALSE, xlab = "x", ylab = "y")
plotBinning(lineBin.chi[[3]], pch = 19, cex = 0.5, add = TRUE,
            fill = residualFill(lineBin.chi[[3]], maxRes = maxRes),
            col = adjustcolor("black", 0.5))
dev.off()


## INVESTIGATING THE NULL DISTRIBUTION ###############################
set.seed(506391)
n <- 1e2 # the sample size
nsim <- 1e4 # number of samples
depths <- 2:10 # range of depths
simDataSets <- replicate(nsim, data.frame(x = sample(1:n),
                                          y = sample(1:n)),
                         simplify = FALSE) # simulated null samples
## allocate storage space
depthSeq.chi <- array(NA, dim = c(4, length(depths), nsim))
depthSeq.mi <- array(NA, dim = c(4, length(depths), nsim))
depthSeq.rnd <- array(NA, dim = c(4, length(depths), nsim))
## set criteria/stop function
crits <- makeCriteria(depth >= dep, expn <= 10, n <= 1)
stopFn <- function(bns) stopper(bns, crits)
## and splitting functions
chiSplit <- function(bn) maxScoreSplit(bn, chiScores, minExp = 5)
miSplit <- function(bn) maxScoreSplit(bn, miScores, minExp = 5)
rndSplit <- function(bn) maxScoreSplit(bn, randScores, minExp = 5)

## iterate through the data and recursively bin
for (ii in 1:nsim) { # each simulation
    for (dep in depths) { # each depth
        depInd <- match(dep, depths) # get storage list index
        ## split using chi scores
        chtr <- binner(simDataSets[[ii]]$x, simDataSets[[ii]]$y,
                       stopper = stopFn, splitter = chiSplit)
        ## compute and store results
        depthSeq.chi[,depInd,ii] <-
            c(chi = binChi(chtr)$stat,
              mi = binMI(chtr)$stat,
              nbin = length(chtr),
              maxDep = max(sapply(chtr,function(bn) bn$depth)))
        ## repeat this for splitting based on the mutual information
        mitr <- binner(simDataSets[[ii]]$x, simDataSets[[ii]]$y,
                       stopper = stopFn, splitter = miSplit)
        depthSeq.mi[,depInd,ii] <- # store results
            c(chi = binChi(mitr)$stat,
              mi = binMI(mitr)$stat,
              nbin = length(mitr),
              maxDep = max(sapply(mitr, function(bn) bn$depth)))
        ## finally, random splits
        rntr <- binner(simDataSets[[ii]]$x, simDataSets[[ii]]$y,
                       stopper = stopFn, splitter = rndSplit)
        depthSeq.rnd[,depInd,ii] <-
            c(chi = binChi(rntr)$stat,
              mi = binMI(rntr)$stat,
              nbin = length(rntr),
              maxDep = max(sapply(rntr, function(bn) bn$depth)))
    }
    ## report back on progress
    if (ii %% 50 == 0) cat(paste("\r Done simulated data set:", ii))
}
## name the data held in each array slice
dimnames(depthSeq.chi) <- list(c("chi", "mi", "nbin", "maxDep"))
dimnames(depthSeq.mi) <- list(c("chi", "mi", "nbin", "maxDep"))
dimnames(depthSeq.rnd) <- list(c("chi", "mi", "nbin", "maxDep"))
## save the data
saveRDS(list(depths = depths, chiSplit = depthSeq.chi,
             miSplit = depthSeq.mi, randSplit = depthSeq.rnd),
        file = paste0("SplitsRandomDatan", n, ".Rds"))

## plot the paths: statistic by number of bins and depth for the three
## different splitting rules (Figs 4.6, 4.7, 4.8)
depthPal <- hcl.colors(9, "Dark 2")
## read in the simulation rather than run it every time
data <- readRDS("SplitsRandomDatan10000.Rds")
depths <- data$depths
for (spltr in c("chiSplit", "miSplit", "randSplit")) {
    png(paste0(spltr, "ChiDepth.png"), width = 4, height = 4,
               units = "in", res = 480)
    narrowPlot(xgrid = seq(0, 3, by = 0.5), # plot region
               xlab = expression(log[10]~"(Number of bins)"),
               ygrid = seq(-1, 4, by = 1),
               ylab = expression(log[10]~{"("~chi^2~statistic~")"}))
    points(log(data[[spltr]]["nbin",,],10), # points
           log(data[[spltr]]["chi",,],10),
           col = adjustcolor(depthPal, 0.1),
           pch = 20)
    ## mean points
    for (ii in 1:9) points(log(mean(data[[spltr]]["nbin",ii,]),10),
                           log(mean(data[[spltr]]["chi",ii,]),10),
                           col = "black", pch = 22,
                           bg = depthPal[ii])
    for (p in c(0.01)) { # chi quantile
        lines(log(2:600,10), log(qchisq(1-p, 1:599),10), lty = 2)
    }
    ## position the legend based on the splitting rule
    if (grepl("rand", spltr)) {
        legpos <- "topleft"
    } else legpos <- "bottomright"
    legend(x = legpos, legend = 2:10, title = "Max depth",
           bg = "white", cex = 0.8,
           fill = adjustcolor(depthPal, 0.6))
    dev.off()
}

## smaller plots to compare different sample sizes (Figs 4.10, 4.11,
## 4.12)
size <- 1.8
for (n in c(1e2, 1e3, 1e4)) {
    ## read in corresponding data
    data <- readRDS(paste0("SplitsRandomDatan", n, ".Rds"))
    depths <- data$depths
    if (n == 1e2) {
        mar <- c(2.1, 2.1, 0.1, 0.1)
        wid <- size + 0.1
    } else if (n == 1e3) {
        mar <- c(2.1, 1.1, 0.1, 0.1)
        wid <- size - 0.1
    } else {
        mar <- c(2.1, 1.1, 0.1, 3.1)
        wid <- size + 0.4
    }
    for (spltr in c("chiSplit", "miSplit", "randSplit")) {
        png(paste0(spltr, "ChiDepth", n, ".png"), width = wid,
            height = size, units = "in", res = 480)
        narrowPlot(xgrid = seq(0, 3, by = 1), ygrid = seq(-1, 4, by = 1),
                   xlab = expression(log[10]~"(Number of bins)"),
                   ylab = expression(log[10]~{"("~chi^2~statistic~")"}),
                   mars = mar)
        points(log(data[[spltr]]["nbin",,],10),
               log(data[[spltr]]["chi",,],10),
               col = adjustcolor(depthPal, 0.1),
               cex = 0.5, pch = 20)
        for (ii in 1:9) points(log(mean(data[[spltr]]["nbin",ii,]),10),
                               log(mean(data[[spltr]]["chi",ii,]),10),
                               col = "black", pch = 22, cex = 0.5,
                               bg = depthPal[ii])
        for (p in c(0.01)) {
            lines(log(2:600,10), log(qchisq(1-p, 1:599),10), lty = 2)
        }
        if (n == 1e4) {
            bnds <- par()$usr
            legend(x = bnds[2], y = bnds[4], legend = 2:10,
                   title = "Max depth", bg = "white",
                   cex = 0.6, xpd = NA,
                   fill = adjustcolor(depthPal, 0.6))
        }
        dev.off()
    }
}

## quantile regression, where the quantiles of the chi statistic are
## predicted based on the number of bins, makes for an easier
## comparison of the upper ends of the distributions
library(quantreg)
qnts <- c(0.95, 0.99, 0.999)
## choose the data with a sample size of 1000
data <- readRDS("SplitsRandomDatan10000.Rds")
depths <- data$depths
depthSeq.chi <- data$chiSplit
depthSeq.mi <- data$miSplit
depthSeq.rnd <- data$randSplit # split sections of the data
## random splitting
rndChiQnt <- rq(chi ~ nbin, tau = qnts,
                data = data.frame(nbin = c(depthSeq.rnd["nbin",,]),
                                  chi = c(depthSeq.rnd["chi",,])))
## chi splitting
chiChiQnt <- rq(chi ~ nbin, tau = qnts,
                data = data.frame(nbin = c(depthSeq.chi["nbin",,]),
                                  chi = c(depthSeq.chi["chi",,])))
## quantile of chi under mi splitting
chiMiQnt <- rq(chi ~ nbin, tau = qnts,
               data = data.frame(nbin = c(depthSeq.mi["nbin",,]),
                                 chi = c(depthSeq.mi["chi",,])))

## plot the quantiles (Fig 4.9)
xseq <- seq(4, 600, by = 1)
randQntPred <- predict(rndChiQnt, newdata = data.frame(nbin = xseq))
chiQntPred <- predict(chiChiQnt, newdata = data.frame(nbin = xseq))
miQntPred <- predict(chiMiQnt, newdata = data.frame(nbin = xseq))
png("binQuantileRegression.png", width = 3, height = 3, res = 480,
    units = "in")
narrowPlot(xgrid = seq(1, 3, by = 1), ylim = c(0, 3.5),
           xlab = expression(log[10]~"(Number of bins)"),
           ygrid = seq(0, 4, by = 1),
           ylab = expression(log[10]~{"("~chi^2~statistic~")"}),
           xlim = c(0.5, 3), mars = c(2.1, 2.1, 0.1, 0.1))
for (ii in 1:3) {
    lines(log(xseq, 10), log(randQntPred[,ii], 10),
          lty = ii, col = adjustcolor("steelblue", 0.8))
    lines(log(xseq, 10), log(chiQntPred[,ii], 10),
          lty = ii, col = adjustcolor("firebrick", 0.8))
    lines(log(xseq, 10), log(miQntPred[,ii], 10),
          lty = ii, col = adjustcolor("seagreen", 0.8))
    lines(log(xseq, 10), log(qchisq(qnts[ii], df = xseq), 10),
          lty = ii, col = adjustcolor("black", 0.8))
    text(log(xseq[1], 10), log(randQntPred[1,ii], 10),
         labels = paste0(100*qnts[ii], "%"), cex = 0.6,
         adj = c(0.5,-0.1))
}
legend(x = "bottomright",
       legend = c("Random bins", "Max MI bins",
                  expression("Max "~chi~bins),
                  expression(chi^2~" critical value")),
       cex = 0.8, col = adjustcolor(c("steelblue", "seagreen",
                                      "firebrick", "black"),
                                    0.8),
       lty = c(1,1,1,1), bg = "white")
dev.off()


## SIMULATED DATA PATTERNS ###########################################
## patterns from Newton (2009)
set.seed(70111238)
n <- 1000
xx <- matrix(NA,n,7)
yy <- matrix(NA,n,7)
## 1
i <- 1
x <- seq( -1, 1, length=n )
u <- x + runif(n)/3
v <-  4*( ( x^2 - 1/2 )^2 + runif(n)/500 )
xx[,i] <- u
yy[,i] <- v
## 2
i <- 2
x <- runif(n, min=(-1), max=1 )
y <- runif(n, min=(-1), max=1 )
theta <- -pi/8
rr <- rbind( c(cos(theta), -sin(theta) ),
             c( sin(theta), cos(theta) ) )
tmp <- cbind( x, y ) %*% rr
u <- tmp[,1]
v <-  tmp[,2]
xx[,i] <- u
yy[,i] <- v
## i=3
i <- 3
x <- runif(n, min=(-1), max=1 )
y <- runif(n, min=(-1), max=1 )
theta <- -pi/4
rr <- rbind( c(cos(theta), -sin(theta) ),
             c( sin(theta), cos(theta) ) )
tmp <- cbind( x, y ) %*% rr
u <- tmp[,1]
v <-  tmp[,2]
xx[,i] <- u
yy[,i] <- v
## 4
i <- 4
x <- seq(-1,1, length=n )
y <- (x ^2 + runif(n))/2
xx[,i] <- x
yy[,i] <- y
## 5
i <- 5
x <- seq(-1,1, length=n )
y <- (x ^2 + runif(n)/2 )*( sample( c(-1,1), size=n, replace=T ) )
xx[,i] <- x
yy[,i] <- y
## 6
i <- 6
x <- seq( -1, 1, length=n )
u <- sin( x*pi ) + rnorm( n )/8
v <- cos( x*pi ) + rnorm( n )/8
xx[,i] <- u
yy[,i] <- v
##
i <- 7
dx <- rnorm(n)/3
dy <- rnorm(n)/3
cx <- sample( c(-1,1), size=n, replace=T )
cy <- sample( c(-1,1), size=n, replace=T )
u <- cx + dx
v <- cy + dy
xx[,i] <- u
yy[,i] <- v

## plot the data (Fig 4.12)
m <- 1
pal <- c(RColorBrewer::brewer.pal(6, "Set2"), "black")
png(file="measurePatterns.png", height=m, width=6*m, units = "in",
    res = 480)
par(mfrow=c(1,7), mar=c(1,1,1,1)/2)
for(i in 1:7)
 {
     plot(xx[,i], yy[,i], xlab="", ylab="", axes=F, pch = 19,
          cex = 0.2, col = pal[i])
 }
dev.off()

## convert this data into pairwise ranks and plot it (Fig 4.13)
xxr <- apply(xx, 2, rank)
yyr <- apply(yy, 2, rank)
png(file="measurePatternsRank.png", height=m, width=6*m, units = "in",
    res = 480)
par(mfrow=c(1,7), mar=c(1,1,1,1)/2)
for(i in 1:7)
 {
     plot(xxr[,i], yyr[,i], xlab="", ylab="", axes=F, pch = 19,
          cex = 0.2, col = pal[i])
 }
dev.off()

## try the binning algorithm on these data
## again define the criteria dynamically
crits <- makeCriteria(depth >= ii, expn <= 10, n <= 1)
## define the stop function
stopFn <- function(bns) stopper(bns, crits)
## chi splitting on patterns
testChiBins <- vector("list", 10)
testMiBins <- vector("list", 9) # allocate storage
for (ii in 1:10) { # iterate through depths
    ## chi bins for each pattern
    testChiBins[[ii]] <- lapply(1:7, function(jj){
        binner(xxr[, jj], yyr[, jj], stopper = stopFn,
               splitter = chiSplit)
    })
    ## mi bins for each pattern
    testMiBins[[ii]] <- lapply(1:7, function(jj){
        binner(xxr[, jj], yyr[, jj], stopper = stopFn,
                  splitter = miSplit)
    })
}

## for random splitting, consider performing many repetitions due to
## the instability of any individual path
set.seed(7211)
testRndBins <- vector("list", 100) # allocate storage
for (jj in 1:100) { # 100 iterations
    testRndBins[[jj]] <- vector("list", 10) # allocate sub-storage
    for (ii in 1:10) { # iterate through depths
        ## randomly split the space for each pattern
        testRndBins[[jj]][[ii]] <- lapply(1:7, function(jj){
            binner(xxr[, jj], yyr[, jj],
                   stopper = stopFn, splitter = rndSplit)
       })
    }
}

## compute the chi square statistics from all
testChiChi <- lapply(testChiBins, function(lst) lapply(lst, binChi))
testMiChi <- lapply(testMiBins, function(lst) lapply(lst, binChi))
testRndChi <- lapply(testRndBins,
                     function(lst) {
                         lapply(lst,
                                function(el) lapply(el, binChi))
                     })

## plot the paths of every pattern under different splitting regimes
## compared to the null
data <- readRDS(paste0("SplitsRandomDatan1000.Rds")) # read in null
depths <- data$depths
depthSeq.chi <- data$chiSplit
depthSeq.mi <- data$miSplit
depthSeq.rnd <- data$randSplit # null data

## plot paths for an individual random split
png("simDataRandPath.png", width = 3, height = 3, units = "in",
    res = 480)
narrowPlot(xgrid = seq(0, 160, by = 40), xlab = "Number of bins",
           ygrid = seq(0, 1600, by = 400),
           ylab = expression(chi^2~statistic))
for (ii in 1:1e4) { # add the null lines
    lines(depthSeq.rnd["nbin",,ii],
          depthSeq.rnd["chi",,ii],
          col = adjustcolor("gray", 0.1))
}
## extract the statistic values from the split
paths <- sapply(testRndChi[[1]],
                function(lst) sapply(lst, function(el) el$stat))
## and the number of bins for each
nbin <- sapply(testRndChi[[1]], function(lst) {
    sapply(lst, function(el) length(el$residuals))
})
## add these as lines to the plot of null lines
for (jj in 1:7) {
    lines(nbin[jj,], paths[jj,], col = pal[jj])
    points(nbin[jj,], paths[jj,], col = pal[jj], pch = 19, cex = 0.5)
}
## add the 95% chi quantile
lines(1:160, qchisq(0.95, 1:160), lty = 2)
dev.off()

## make the same plot for paths from chi splitting
png("simDataMaxChiPath.png", width = 3, height = 3, units = "in",
    res = 480)
narrowPlot(xgrid = seq(0, 160, by = 40),
           xlab = "Number of bins",
           ygrid = seq(0, 1600, by = 400),
           ylab = expression(chi^2~statistic))
for (ii in 1:1e4) {
    lines(depthSeq.chi["nbin",,ii],
          depthSeq.chi["chi",,ii],
          col = adjustcolor("gray", 0.1))
}
paths <- sapply(testChiChi,
               function(lst) sapply(lst, function(el) el$stat))
nbin <- sapply(testChiChi, function(lst) {
    sapply(lst, function(el) length(el$residuals))
})
for (jj in 1:7) {
    lines(nbin[jj,], paths[jj,], col = pal[jj])
    points(nbin[jj,], paths[jj,], col = pal[jj], pch = 19, cex = 0.5)
}
dev.off()

## for the 100 random split repetitions, plot every single one
png("simDataRandAll.png", width = 3, height = 3, units = "in",
    res = 480)
narrowPlot(xgrid = seq(0, 160, by = 40),
           xlab = "Number of bins",
           ygrid = seq(0, 1200, by = 300), ylim = c(0, 1300),
           ylab = expression(chi^2~statistic))
for (ii in 1:1e4) {
    lines(depthSeq.rnd["nbin",,ii],
          depthSeq.rnd["chi",,ii],
          col = adjustcolor("gray", 0.1))
}
paths <- simplify2array(
    lapply(testRndChi, function(olst) {
        sapply(olst, function(lst) {
            sapply(lst, function(el) el$stat)
        })
    }))
nbin <- simplify2array(
    lapply(testRndChi, function(olst) {
        sapply(olst, function(lst) {
            sapply(lst, function(el) length(el$residuals))
        })
    }))
for (jj in 1:7) {
    for (ii in 1:100) {
        lines(nbin[jj,,ii], paths[jj,,ii],
              col = adjustcolor(pal[jj], 0.2))
    }
}
lines(1:160, qchisq(0.95, 1:160), lty = 2)
dev.off()

## let's see what the bins look like for every depth
## start by getting the maximum residual to make the shading constant
maxRes <- max(abs(unlist(sapply(unlist(testChiChi,
                                       recursive = FALSE),
                                function(el) el$residuals))))
## for every depth, display the binning for each pattern
for (depth in 2:10) {
    png(file = paste0("simDataBins", depth, ".png"), height = m,
        width= 6*m, units = "in", res = 480)
    par(mfrow=c(1,7), mar=c(1,1,1,1)/2)
    for(i in 1:7) {
        plotBinning(testChiBins[[depth]][[i]], xlab="", ylab="",
                    pch = 19, cex = 0.1, axes = F,
                    col = adjustcolor(pal[i], 0.5),
                    fill = residualFill(testChiBins[[depth]][[i]],
                                        maxRes = maxRes))
    }
    dev.off()
}

## do the same thing for the bins from MI splitting
maxRes <- max(abs(unlist(sapply(unlist(testMiChi,
                                       recursive = FALSE),
                                function(el) el$residuals))))
## for every depth, display the binning for each pattern
for (depth in 2:10) {
    png(file = paste0("simDataBins", depth, "MI.png"), height = m,
        width= 6*m, units = "in", res = 480)
    par(mfrow=c(1,7), mar=c(1,1,1,1)/2)
    for(i in 1:7) {
        plotBinning(testMiBins[[depth]][[i]], xlab="", ylab="",
                    pch = 19, cex = 0.1, axes = F,
                    col = adjustcolor(pal[i], 0.5),
                    fill = residualFill(testMiBins[[depth]][[i]],
                                        maxRes = maxRes))
    }
    dev.off()
}

## view the random bins a their maximum depth
## again, standardize the residuals
maxRes <- max(abs(unlist(sapply(unlist(testRndChi[[1]],
                                       recursive = FALSE),
                                function(el) el$residuals))))
depth <- 10
png(file = "simDataBinsRand.png", height = m, width = 6*m,
    units = "in", res = 480)
par(mfrow=c(1,7), mar=c(1,1,1,1)/2)
for(i in 1:7) {
    plotBinning(testRndBins[[10]][[depth]][[i]], xlab="", ylab="",
                pch = 19, cex = 0.1, axes = F,
                col = adjustcolor(pal[i], 0.5),
                fill = residualFill(testRndBins[[10]][[depth]][[i]],
                                    maxRes = maxRes))
}
dev.off()


## REAL DATA EXAMPLE #################################################
## S&P500 data: "SP500" demo in "zenplots" package
## abalone data
url <- "https://archive.ics.uci.edu/ml/machine-learning-databases/abalone/"
## loading here is locally, otherwise replace "~/Downloads/" with the url above
abalone <- read.table("~/Downloads/abalone.data", sep = ",")
abalone.names <- scan("~/Downloads/abalone.names", what = character(),
                      sep = "\n")[64:72]
names(abalone) <- regmatches(gsub(" ", ".", abalone.names),
                             regexpr("(?<=^\\t)[A-Za-z\\.]+(?=\\t)",
                                     gsub(" ", ".", abalone.names),
                                     perl = TRUE))
## take ranks, get all pairs
set.seed(2306201)
abPairs <- combn(ncol(abalone), 2)
abRanks <- lapply(abalone, rank, ties.method = "random")
abBins <- lapply(1:ncol(abPairs),
                 function(p) {
                     lapply(depths,
                            function(dep) {
                     makeTree(abRanks[[abPairs[1,p]]],
                              abRanks[[abPairs[2,p]]],
                              stopCriterion(depthLim = dep,
                                            areaLim = 10),
                              maxScoreSplit(chiScores), quickBin)
                 })})
names(abBins) <- paste(names(abRanks)[abPairs[1,]],
                       names(abRanks)[abPairs[2,]],
                       sep = "-")
abChis <- lapply(abBins, function(ds) lapply(ds, function(bnlst) binChi(unlist(bnlst))))
abMis <- lapply(abBins, function(ds) lapply(ds, function(bnlst) binMI(unlist(bnlst))))

## plot these against other points
plot(depthSeq.chi["nbin",,], depthSeq.chi["chi",,], xlab = "Number of bins",
     ylab = "", col = adjustcolor("gray", 0.5), pch = 20,
     ylim = c(0, max(sapply(abChis, function(lst) sapply(lst, function(scr) scr$stat)))),
     xlim = c(0, 275))
mtext(expression(chi^2~statistic), side = 2, line = 2)
lines(c(rbind(sapply(abChis[grepl("Rings", names(abChis))],
                     function(lst) sapply(lst, function(scr) length(scr$residuals))),NA)),
      c(rbind(sapply(abChis[grepl("Rings", names(abChis))],
                     function(lst) sapply(lst, function(scr) scr$stat)), NA)),
      pch = 20, col = adjustcolor("#fc8d62", 0.8), type = "b")
lines(c(rbind(sapply(abChis[grepl("Sex", names(abChis))],
                     function(lst) sapply(lst, function(scr) length(scr$residuals))),NA)),
      c(rbind(sapply(abChis[grepl("Sex", names(abChis))],
                     function(lst) sapply(lst, function(scr) scr$stat)), NA)),
      pch = 20, col = adjustcolor("#66c2a5", 0.8), type = "b")
lines(c(rbind(sapply(abChis[!grepl("Rings|Sex", names(abChis))],
                     function(lst) sapply(lst, function(scr) length(scr$residuals))),NA)),
      c(rbind(sapply(abChis[!grepl("Rings|Sex", names(abChis))],
                     function(lst) sapply(lst, function(scr) scr$stat)), NA)),
      pch = 20, col = adjustcolor("#8da0cb", 0.8), type = "b")
legend(x = "topleft", pch = 20, legend = c("Continuous", "Categorical", "Count", "Null"),
       col = adjustcolor(c("#8da0cb","#66c2a5","#fc8d62","gray"), 0.8),
       title = "Data type")

## ... also using MI
plot(depthSeq.chi["nbin",,], depthSeq.chi["mi",,], xlab = "Number of bins",
     ylab = "Mutual information", col = adjustcolor("gray", 0.5), pch = 20,
     ylim = c(0, max(sapply(abMis, function(lst) sapply(lst, function(scr) scr$stat)))),
     xlim = c(0, 275))
lines(c(rbind(sapply(abMis[grepl("Rings", names(abMis))],
                     function(lst) sapply(lst, function(scr) length(scr$residuals))),NA)),
      c(rbind(sapply(abMis[grepl("Rings", names(abMis))],
                     function(lst) sapply(lst, function(scr) scr$stat)), NA)),
      pch = 20, col = adjustcolor("#fc8d62", 0.8), type = "b")
lines(c(rbind(sapply(abMis[grepl("Sex", names(abMis))],
                     function(lst) sapply(lst, function(scr) length(scr$residuals))),NA)),
      c(rbind(sapply(abMis[grepl("Sex", names(abMis))],
                     function(lst) sapply(lst, function(scr) scr$stat)), NA)),
      pch = 20, col = adjustcolor("#66c2a5", 0.8), type = "b")
lines(c(rbind(sapply(abMis[!grepl("Rings|Sex", names(abMis))],
                     function(lst) sapply(lst, function(scr) length(scr$residuals))),NA)),
      c(rbind(sapply(abMis[!grepl("Rings|Sex", names(abMis))],
                     function(lst) sapply(lst, function(scr) scr$stat)), NA)),
      pch = 20, col = adjustcolor("#8da0cb", 0.8), type = "b")
legend(x = "topleft", pch = 20, legend = c("Continuous", "Categorical", "Count", "Null"),
       col = adjustcolor(c("#8da0cb","#66c2a5","#fc8d62","gray"), 0.8),
       title = "Data type")

## get some pairwise values
linex <- 1:(nrow(abalone))
liney <- 1:(nrow(abalone))
lineBin <- makeTree(linex, liney, # to scale values
                    stopCriterion(depthLim = 6, areaLim = 10),
                    maxScoreSplit(chiScores, ties = "x"), quickBin)
lineChi <- binChi(unlist(lineBin))
pairwiseAb <- matrix(lineChi$stat, ncol(abalone), ncol(abalone))
for (ii in 1:ncol(abPairs)) {
    tempC <- abPairs[,ii]
    tempChi <- abChis[[ii]]
    pairwiseAb[tempC[1], tempC[2]] <- pairwiseAb[tempC[2], tempC[1]] <- tempChi[[5]]$stat
}
pairwiseAb <- pairwiseAb/lineChi$stat


## definitely some stronger relationships than uniform...
## also, note the strong ordering for both, plot all of these
abOrders <- sapply(1:length(depths),
                   function(ii){
                       order(sapply(abChis, function(scr) scr[[ii]]$stat),
                             decreasing = TRUE)
                   })
maxRes <- sapply(1:length(depths),
                 function(ii) {
                     max(abs(unlist(sapply(abChis, function(scr) scr[[ii]]$residuals))))
                 })
ind <- 5
par(mfrow = c(6,6), mar = c(1.1,1.1,2.1,1.1))
for (ii in 1:length(abOrders[,ind])) {
    tempBins <- unlist(abBins[abOrders[,ind]][[ii]][[ind]])
    plotBinning(tempBins, pch = "", xaxt = "n", yaxt = "n",
                residualFill(tempBins, maxRes = maxRes[ind]),
                xaxt = "n", yaxt = "n",
                main = names(abBins)[abOrders[,ind]][ii])
}

## the iris data
data(iris)
wid <- rank(iris$Sepal.Width, ties.method = "random")
len <- rank(iris$Sepal.Length, ties.method = "random")
widLen <- makeTree(wid, len,
                   stopCriterion(depthLim = 5, areaLim = 5),
                   maxScoreSplit(chiScores, ties = "x"), quickBin)
widLen.mi <- makeTree(wid, len,
                   stopCriterion(depthLim = 5, areaLim = 5),
                   maxScoreSplit(miScores, ties = "x"), quickBin)

## some simulated structural data
data <- data.frame(x = c(sample(rnorm(2000, c(-10, 10), c(4,3))), rnorm(400)),
                   y = c(sample(rnorm(2000, c(-10, 10), c(1,5))), rnorm(400)))
fourCent <- makeTree(rank(data$x), rank(data$y),
                     stopCriterion(depthLim = 5, areaLim = 5),
                     maxScoreSplit(chiScores, ties = "x"), quickBin)

## upper corner
upCorn <- list()
upCorn$x <- runif(1000)
upCorn$x <- (1 - sqrt(upCorn$x))*1000
upCorn$y <- upCorn$x + runif(1000)*(1000 - upCorn$x)
upCornChi <- makeTree(upCorn$x, upCorn$y,
                     stopCriterion(depthLim = 5, areaLim = 5),
                     maxScoreSplit(chiScores, ties = "x"), quickBin)
