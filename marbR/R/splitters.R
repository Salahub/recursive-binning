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

## halve a bin
halfSplit <- function(bin, margin = "x") {
    if (margin == "x") {
        xsort <- order(bin$x)
        hind <- floor(bin$n/2) # middle index
        newbnd <- bin$x[xsort][hind] # middle value
        above <- xsort[(hind+1):(bin$n)] # points above
        list(list(x = bin$x[-above], y = bin$y[-above], # new bins
                  bnds = list(x = c(bin$bnds$x[1], newbnd),
                              y = bin$bnds$y),
                  area = bin$area*(newbnd - bin$bnds$x[1])/
                      diff(bin$bnds$x),
                  n = bin$n-length(above), depth = bin$depth + 1),
             list(x = bin$x[above], y = bin$y[above],
                  bnds = list(x = c(newbnd, bin$bnds$x[2]),
                              y = bin$bnds$y),
                  area = bin$area*(bin$bnds$x[2] - newbnd)/
                      diff(bin$bnds$x),
                  n = length(above), depth = bin$depth + 1))
    } else if (margin == "y") {
        ysort <- order(bin$y)
        hind <- floor(bin$n/2) # middle index
        newbnd <- bin$y[ysort][hind] # middle value
        above <- ysort[(hind+1):(bin$n)] # points above
        list(list(x = bin$x[-above], y = bin$y[-above],
                  bnds = list(x = bin$bnds$x,
                              y = c(bin$bnds$y[1], newbnd)),
                  area = bin$area*(newbnd - bin$bnds$y[1])/
                      diff(bin$bnds$y),
                  n = bin$n - length(above), depth = bin$depth + 1),
             list(x = bin$x[above], y = bin$y[above],
                  bnds = list(x = bin$bnds$x,
                              y = c(newbnd, bin$bnds$y[2])),
                  area = bin$area*(bin$bnds$y[2] - newbnd)/
                      diff(bin$bnds$y),
                  n = length(above), depth = bin$depth + 1))
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
  if (xscore[xmax] >= yscore[ymax]) { # ties go to x
      xsplts <- bin$x[xsort]
      newbnd <- c(xsplts[1]-1, xsplts)[xmax] # new boundary
      below <- xsort[seq_len(xmax-1)] # get indices of points below
      above <- if (xmax == bin$n+1) integer(0) else xsort[xmax:bin$n]
      list(list(x = bin$x[below], y = bin$y[below], # new bins
                bnds = list(x = c(bin$bnds$x[1], newbnd),
                            y = bin$bnds$y),
                area = bin$area*(newbnd - bin$bnds$x[1])/
                    diff(bin$bnds$x),
                n = length(below), depth = bin$depth + 1),
           list(x = bin$x[above], y = bin$y[above],
                bnds = list(x = c(newbnd, bin$bnds$x[2]),
                            y = bin$bnds$y),
                area = bin$area*(bin$bnds$x[1] - newbnd)/
                    diff(bin$bnds$x),
                n = length(above), depth = bin$depth + 1))
  } else { # do the same on y
      ysplts <- bin$y[ysort]
      newbnd <- c(ysplts[1]-1, ysplts)[ymax]
      below <- ysort[seq_len(xmax-1)]
      above <- if (ymax == bin$n+1) integer(0) else ysort[ymax:bin$n]
      list(list(x = bin$x[below], y = bin$y[below],
                bnds = list(x = bin$bnds$x,
                            y = c(bin$bnds$y[1], newbnd)),
                area = bin$area*(newbnd - bin$bnds$y[1])/
                    diff(bin$bnds$y),
                n = length(below), depth = bin$depth + 1),
           list(x = bin$x[above], y = bin$y[above],
                bnds = list(x = bin$bnds$x,
                            y = c(newbnd, bin$bnds$y[2])),
                area = bin$area*(bin$bnds$y[2] - newbnd)/
                    diff(bin$bnds$y),
                n = length(above), depth = bin$depth + 1))
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
  list(list(x = bin$x[below], y = bin$y[below], # new bins
            bnds = list(x = c(bin$bnds$x[1], newbnd),
                        y = bin$bnds$y),
            area = bin$area*(newbnd - bin$bnds$x[1])/
                diff(bin$bnds$x),
            n = length(below), depth = bin$depth + 1),
       list(x = bin$x[above], y = bin$y[above],
            bnds = list(x = c(newbnd, bin$bnds$x[2]),
                        y = bin$bnds$y),
            area = bin$area*(bin$bnds$x[2] - newbnd)/
                diff(bin$bnds$x),
            n = length(above), depth = bin$depth + 1))
}

## NB: half/average splits don't really make sense here... we are in
## the rank space and so no points can possibly exist at half
