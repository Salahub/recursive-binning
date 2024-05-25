## build a prototype single split function to envision how the
## algorithm might work

## stop functions
makeCriteria <- function(...) {
    cl <- match.call() # capturing inputs
    crits <- as.list(cl) # change to a list
    ## remove self reference, collapse into single OR
    paste(sapply(crits[-1], deparse), collapse = " | ")
}
stopper <- function(binList, criteria) {
    sapply(binList,
           function(b) eval(parse(text = criteria), envir = b))
}

## functions to split the chosen bin
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
##' @describeIn marginalsplitters Splitting on y
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

## something to choose a split location
randomSplit <- function(bin, bound, N) {
  xsort <- order(bin$x)
  ysort <- order(bin$y) # get marginal ordering
  deltax <- diff(bin$bnds$x)
  deltay <- diff(bin$bnds$y)
  u <- runif(1)
  if (u <= 0.5) { # split on x
      lower <- bin$bnds$x[1] + ceiling(N*bound/deltay)
      upper <- bin$bnds$x[2] - ceiling(N*bound/deltay)
      if (upper <= lower) {
          list("STOP" = bin)
      } else{
          locs <- seq(from = lower, to = upper, by = 1)
          spos <- sample(locs, size = 1)
          above <- which(bin$x > spos)
          below <- which(bin$x <= spos)
          splitX(bin, spos, above, below)
      }
  } else {
      lower <- bin$bnds$y[1] + ceiling(N*bound/deltax)
      upper <- bin$bnds$y[2] - ceiling(N*bound/deltax)
      if (upper <= lower) {
          list("STOP" = bin)
      } else{
          locs <- seq(from = lower, to = upper, by = 1)
          spos <- sample(locs, size = 1)
          above <- which(bin$y > spos)
          below <- which(bin$y <= spos)
          splitY(bin, spos, above, below)
      }
  }
}

## assume we have a list of areas in decreasing order, choose a bin
## to split assuming splits must be strictly greater than the bound
## provided
chooseBin <- function(areas) {
    sample(seq_along(areas), size = 1, prob = areas)
    ## cute but not efficient for R
    #n <- length(areas)
    #cands <- numeric(n) # set storage
    #ii <- 1 # set indicator
    #parsum <- 0
    #while (areas[ii] >= 2*bound & ii <= n) {
    #    cands[ii] <- 1
    #    parsum <- parsum + areas[ii]
    #    ii <- ii + 1
    #}
    #u <- runif(1)
    #jj <- 0
    #part <- 0
    #while (u > part/parsum) {
    #    jj <- jj + 1
    #    part <- part + areas[jj]
    #}
    #return(jj)
}

## place an area in a sorted list
placeBin <- function(areas, newa, lower = 1) {
    upper <- length(areas)
    if (newa < areas[upper]) {
        length(areas) + 1
    } else if (newa > areas[lower]) {
        1
    } else {
        while (upper - lower > 1) {
            test <- floor((upper + lower)/2)
            if (newa < areas[test]) {
                lower <- test
            } else {
                upper <- test
            }
        }
        lower + 1
    }
}

## helper to place given an index
placeAtIndex <- function(vec, x, ind = 1) {
    len <- length(vec)
    if (ind == 1) {
        c(x, vec)
    } else if (ind > len) {
        c(vec, x)
    } else {
        c(vec[1:(ind-1)], x, vec[ind:len])
    }
}

## the outside wrapper function
randBinner <- function(x, y, stopper, minExp = 5, maxK = 5) {
    ## initialize bin with all the data contained
    bin <- list(x = x, y = y, # x and y points
                 bnds = list(x = c(0, max(x, na.rm = TRUE)),
                             y = c(0, max(y, na.rm = TRUE))),
                 expn = length(x), # default expectation is n
                n = length(x), depth = 0) # size, depth
    N <- length(x)
    binList <- randomSplit(bin, minExp, length(x))
    K <- 2
    areas <- sapply(binList, function(el) el$expn)
    binList <- binList[order(areas, decreasing = TRUE)]
    areas <- areas[order(areas, decreasing = TRUE)]
    stopStatus <- stopper(binList) # initialize logical vector

    while (any(!stopStatus) & K < maxK) { # check the stop criteria
        chosen <- chooseBin(areas[!stopStatus])
        newBins <- randomSplit(binList[!stopStatus][[chosen]],
                               minExp, N)
        if (identical(names(newBins), "STOP")) { # no valid splits
            stopStatus[chosen] <- TRUE
        } else {
            areas <- areas[-chosen] # remove split bin
            binList <- binList[-chosen]
            stopStatus <- stopStatus[-chosen]
            ord <- order(c(newBins[[1]]$expn, newBins[[2]]$expn,
                           areas))
            areas <- c(newBins[[1]]$expn, newBins[[2]]$expn,
                       areas)[ord]
            binList <- c(newBins, binList)[ord]
            stopStatus <- c(stopper(newBins), stopStatus)[ord]
            ## cute: but inefficient(?)
            ##for (bin in newBins) { # place new bins
            ##    ind <- placeBin(areas, bin$expn)
            ##    areas <- placeAtIndex(areas, bin$expn, ind)
            ##    binList <- placeAtIndex(binList, list(bin), ind)
            ##    stopStatus <- placeAtIndex(stopStatus,
            ##                               stopper(list(bin)), ind)
            ##}
            K <- K + 1
        }
    }

    binList # return the final list of bins
}


## SCRIPT ##########################################################
n <- 1e4
Ks <- seq(4, 400, by = 2)
nsim <- 100

## set criteria
stopCrits <- makeCriteria(expn <= 10, n < 1, depth == 10)
stopFun <- function(bns) stopper(bns, stopCrits)

## sample data
set.seed(24052024)
samples <- replicate(nsim*length(Ks),
                     list(x = sample(1:n), y = sample(1:n)),
                     simplify = FALSE)

## go through Ks and split
singleSplitBins <- vector(mode = "list", length = length(Ks)*nsim)
for (k in Ks) {
    for (ii in 1:nsim) {
        ind <- (k/2 - 2)*nsim + ii
        dat <- samples[[ind]]
        binning <- randBinner(dat$x, dat$y, stopper = stopFun,
                              minExp = 5, maxK = k)
        singleSplitBins[[ind]] <- binning
    }
    if ((k %% 50) == 0) cat("\r Done to", k)
}

## get stats
library(AssocBin)
## collect info
nbin <- sapply(singleSplitBins, length)
tests <- lapply(singleSplitBins, binChi)
stats <- sapply(tests, function(x) x$stat)
plot(log(nbin, 10), log(stats, 10), ylim = c(-1,3.5), pch = 19,
     col = adjustcolor("black", 0.2), cex = 0.5,
     main = "Splitting one bin at a time")
lines(log(seq(1, 400, by = 1), 10), lwd = 2,
      log(qchisq(0.99, df = seq(1, 400, by = 1)), 10),
      col = "firebrick", lty = 2)
