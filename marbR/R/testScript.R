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
chiScores <- function(vals) {
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
    scr
}
## mutual informations
miScores <- function(vals) {
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
    below + above
}
## random scoring for random splits
randScores <- function(vals) {
    runif(length(vals))+1
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

## splitter maximizing a score function
maxScoreSplit <- function(bin, scorer = diff) {
  xsort <- order(bin$x)
  ysort <- order(bin$y) # get marginal ordering
  xscore <- scorer(c(bin$bnds$x[1], bin$x[xsort], bin$bnds$x[2]))
  yscore <- scorer(c(bin$bnds$y[1], bin$y[ysort], bin$bnds$y[2]))
  xmax <- which.max(xscore)
  ymax <- which.max(yscore) # the score values
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
uniMaxScoreSplit <- function(bin, scorer = diff) {
  xsort <- order(bin$x)
  xscore <- scorer(c(bin$bnds$x[1], bin$x[xsort], bin$bnds$x[2]))
  xmax <- which.max(xscore)
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

## something to plot/visualize bins
plotBinning <- function(bins, xlab = "x", ylab = "y",
                        main = "Bins", ...) {
    nbins <- length(bins) # number of bins
    xbnds <- sapply(bins, function(bn) bn$bnds$x)
    ybnds <- sapply(bins, function(bn) bn$bnds$y)
    plot(NA, type = "n", xlim = range(xbnds), ylim = range(ybnds),
         xlab = xlab, ylab = ylab, main = main)
    for (ii in 1:nbins) {
        rect(xbnds[1,ii], ybnds[1,ii], xbnds[2,ii], ybnds[2,ii])
        points(bins[[ii]]$x, bins[[ii]]$y, ...)
    }
}


## TRY THESE OUT
criteria <- makeCriteria(expn <= 10, n <= 10, depth >= 10)
test <- binner(1:100, 1:100,
               stopper = function(bns) stopper(bns, criteria),
               splitter = function(bin) maxScoreSplit(bin, chiScores))
