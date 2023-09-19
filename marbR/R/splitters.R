
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

## another function which halves a bin independent of the point
## locations to break ties and limit bin expecations
halfCutTie <- function(bin, xscore, yscore) {
    u <- as.numeric(yscore > xscore) # prefer to split on max score
    if (yscore == xscore) u <- runif(1)
    if (u < 0.5) { # y has a larger score, or random
        newbnd <- ceiling(mean(bin$bnds$x)) # split value
        abv <- bin$x > newbnd # which x values are above
        above <- which(abv) # indices above
        below <- which(!abv) # indices below
        splitX(bin, bd = newbnd, above = above, below = below)
    } else {
        newbnd <- ceiling(mean(bin$bnds$y)) # split value
        abv <- bin$y > newbnd # which y values are above
        above <- which(abv) # indices above
        below <- which(!abv) # indices below
        splitY(bin, bd = newbnd, above = above, below = below)
    }
}

## splitter maximizing a score function
maxScoreSplit <- function(bin, scorer, ties = halfCutTie,
                          pickMax = which.max, ...) {
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
  if (xallEq & yallEq) { # in the case of ties, use tie function
      ties(bin, xscore[1], yscore[1])
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
