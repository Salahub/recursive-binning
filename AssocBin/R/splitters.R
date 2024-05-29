##` Marginalsplitters
##' @title Helper functions for marginal splitting
##' @description These functions are helpers to safely split bins
##' along X or Y.
##' @details These unexported functions have been defined primarily
##' to clean up other code, but could be changed to obtain different
##' core functionality.
##' @param bin a bin to be split with elements `x`, `y`, `depth`,
##' `bnds` (list with elements `x` and `y`), `expn`, `n`
##' @param bd numeric split point within the bin bounds
##' @param above indices of `x` and `y` points in the bin above `bd`
##' @param below indices of `x` and `y` points in the bin below `bd`
##' @return A list of two bins resulting from the split of `bin` at
##' `bds`.
##' @author Chris Salahub
##' @describeIn marginalsplitters Splitting on x
splitX <- function(bin, bd, above, below) {
    belowfac <- (bd - bin$bnds$x[1])/diff(bin$bnds$x)
    abovefac <- (bin$bnds$x[2] - bd)/diff(bin$bnds$x)
    list(list(x = bin$x[below], y = bin$y[below],
              bnds = list(x = c(bin$bnds$x[1], bd),
                          y = bin$bnds$y),
              expn = bin$expn*belowfac,
              n = bin$n-length(above), depth = bin$depth + 1,
              stopped = FALSE),
         list(x = bin$x[above], y = bin$y[above],
              bnds = list(x = c(bd, bin$bnds$x[2]),
                          y = bin$bnds$y),
              expn = bin$expn*abovefac,
              n = length(above), depth = bin$depth + 1,
              stopped = FALSE))
}
##' @describeIn marginalsplitters Splitting on y
splitY <- function(bin, bd, above, below) {
    belowfac <- (bd - bin$bnds$y[1])/diff(bin$bnds$y)
    abovefac <- (bin$bnds$y[2] - bd)/diff(bin$bnds$y)
    list(list(x = bin$x[below], y = bin$y[below],
              bnds = list(x = bin$bnds$x,
                          y = c(bin$bnds$y[1], bd)),
              expn = bin$expn*belowfac,
              n = bin$n-length(above), depth = bin$depth + 1,
              stopped = FALSE),
         list(x = bin$x[above], y = bin$y[above],
              bnds = list(x = bin$bnds$x,
                          y = c(bd, bin$bnds$y[2])),
              expn = bin$expn*abovefac,
              n = length(above), depth = bin$depth + 1,
              stopped = FALSE))
}

##' @title Halve at an observed point
##' @description This function halves a bin under the restriction that
##' splits can only occur at observation coordinates.
##' @details Given a bin and a margin, this function splits the bin so
##' half the points are above the new split point and half are below.
##' @param bin a bin to be split with elements `x`, `y`, `depth`,
##' `bnds` (list with elements `x` and `y`), `expn`, `n`
##' @param margin string, one of `x` or `y`
##' @return A list of two bins resulting from the split of `bin` in
##' half along the specified margin
##' @examples
##' bin <- list(x = 1:10, y = sample(1:10),
##'             bnds = list(x = c(0, 10), y = c(0, 10)),
##'             expn = 10, n = 10, depth = 0)
##' halfSplit(bin)
##' halfSplit(bin, margin = "y")
##' @author Chris Salahub
halfSplit <- function(bin, margin = sample(c("x", "y"), 1)) {
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

##' @title Halve continuously to break ties
##' @description This function halves a bin based on the midpoint of
##' the bounds along whichever margin produces the larger score.
##' @details The goal of this function is to break ties within bin
##' splitting in a way which prevents very small or lopsided bins from
##' forming, a common problem with the `halfSplit` function
##' @param bin a bin to be split with elements `x`, `y`, `depth`,
##' `bnds` (list with elements `x` and `y`), `expn`, `n`
##' @param xscore numeric value giving the score for all splits along
##' x
##' @param yscore numeric value giving the score for all splits along
##' y
##' @param wider logical; is the bin wider than it is tall?
##' @param splitLongSide logical value, should we force splitting on
##' the longer side regardless of scores?
##' @return A list of two bins resulting from the split of `bin` in
##' half along the margin corresponding to the larger score.
##' @examples
##' bin <- list(x = 1:10, y = sample(1:10),
##'             bnds = list(x = c(0, 10), y = c(0, 10)),
##'             expn = 10, n = 10, depth = 0)
##' halfCutTie(bin, 1, 2, wider = FALSE) # splits on y
##' halfCutTie(bin, 2, 1, wider = FALSE) # splits on x
##' halfCutTie(bin, 1, 1, wider = FALSE) # ties are random
##' @author Chris Salahub
halfCutTie <- function(bin, xscore, yscore, wider,
                       splitLongSide = FALSE) {
    u <- as.numeric(yscore > xscore) # prefer to split on max score
    if (yscore == xscore) u <- runif(1)
    if (u < 0.5 | (splitLongSide & !wider)) { # y has a larger score
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

##' @title Bivariate score maximizing splitting
##' @description A function which splits a bin based on the location
##' maximizing a score function.
##' @details This function serves as a wrapper which manages the
##' interaction of a score function, marginal splitting functions,
##' tie breaking function, and a maximum selection function to split
##' a bin at the observation coordinate which maximizes the score
##' function.
##' @param bin a bin to be split with elements `x`, `y`, `depth`,
##' `bnds` (list with elements `x` and `y`), `expn`, `n`
##' @param scorer function which accepts a numeric vector of potential
##' split coordinates and the bounds of `bin` and returns a numeric
##' vector of scores for each
##' @param ties function which is called to break ties when all splits
##' generate the same score
##' @param pickMax function which accepts a list of scores and returns
##' the element of the largest score according to some rule
##' @param splitLongSide logical value, should we force splitting on
##' the longer side regardless of scores?
##' @param ... optional additional arguments to `scorer`
##' @return A list of two bins resulting from the split of `bin`
##' along the corresponding margin at the maximum location
##' @examples
##' bin <- list(x = 1:10, y = sample(1:10),
##'             bnds = list(x = c(0, 10), y = c(0, 10)),
##'             expn = 10, n = 10, depth = 0)
##' maxScoreSplit(bin, chiScores)
##' maxScoreSplit(bin, miScores) # pretty similar for both
##' maxScoreSplit(bin, randScores)
##' maxScoreSplit(bin, randScores) # different every time
##' @author Chris Salahub
maxScoreSplit <- function(bin, scorer, ties = halfCutTie,
                          pickMax = which.max, splitLongSide = FALSE,
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
  wider <- (bin$bnds$x[2] - bin$bnds$x[1]) >
      (bin$bnds$y[2] - bin$bnds$y[1])
  if (xallEq & yallEq) { # in the case of ties, use tie function
      ties(bin, xscore[1], yscore[1], wider = wider,
           splitLongSide = splitLongSide)
  } else if (xscore[xmax] >= yscore[ymax] |
             (splitLongSide & wider)) { # ties go to x
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

##' @title Random integer splitting
##' @description A function which splits a bin at a random integer
##' conforming to limits on minimum bin size.
##' @details This function serves as a wrapper which manages the
##' interaction of a score function, marginal splitting functions,
##' tie breaking function, and a maximum selection function to split
##' a bin at the observation coordinate which maximizes the score
##' function.
##' @param bin a bin to be split with elements `x`, `y`, `depth`,
##' `bnds` (list with elements `x` and `y`), `expn`, `n`
##' @param minExp numeric giving the minimum expected count allowed
##' in a bin
##' @param splitLongSide logical value, should we force splitting on
##' the longer side regardless of scores?
##' @param ... optional additional arguments (for compatibility)
##' @return A list of two bins resulting from the split of `bin`
##' along the corresponding margin at the maximum location
##' @examples
##' bin <- list(x = 1:10, y = sample(1:10),
##'             bnds = list(x = c(0, 10), y = c(0, 10)),
##'             expn = 10, n = 10, depth = 0)
##' rIntSplit(bin, minExp = 2)
##' @author Chris Salahub
rIntSplit <- function(bin, minExp = 5, splitLongSide = FALSE, ...) {
  expn <- bin$expn
  prop <- minExp/expn
  deltax <- diff(bin$bnds$x)
  deltay <- diff(bin$bnds$y)  
  wider <- (bin$bnds$x[2] - bin$bnds$x[1]) >
      (bin$bnds$y[2] - bin$bnds$y[1])
  u <- runif(1)
  if (u <= 0.5 | (splitLongSide & wider)) { # split on x
      lower <- bin$bnds$x[1] + ceiling(prop*deltax)
      upper <- bin$bnds$x[2] - ceiling(prop*deltax)
      if (upper <= lower) {
          bin$stopped <- TRUE
          list(bin)
      } else{
          locs <- seq(from = lower, to = upper, by = 1)
          spos <- sample(locs, size = 1)
          above <- which(bin$x > spos)
          below <- which(bin$x <= spos)
          splitX(bin, spos, above, below)
      }
  } else {
      lower <- bin$bnds$y[1] + ceiling(prop*deltay)
      upper <- bin$bnds$y[2] - ceiling(prop*deltay)
      if (upper <= lower) {
          bin$stopped <- TRUE
          list(bin)
      } else{
          locs <- seq(from = lower, to = upper, by = 1)
          spos <- sample(locs, size = 1)
          above <- which(bin$y > spos)
          below <- which(bin$y <= spos)
          splitY(bin, spos, above, below)
      }
  }
}

##' @title Random uniform splitting
##' @description Split bins randomly and uniformly
##' @details This function samples a coordinate uniformly along a
##' random margin and splits a bin at that coordinate. In contrast to
##' maxScoreSplit with randScores, this can introduce splits at
##' locations other than the points.
##' @param bin a bin to be split with elements `x`, `y`, `depth`,
##' `bnds` (list with elements `x` and `y`), `expn`, `n`
##' @param minExp numeric giving the minimum expected count allowed
##' in a bin
##' @param splitLongSide logical value, should we force splitting on
##' the longer side?
##' @param ... optional additional arguments (for compatibility)
##' @return A list of two bins resulting from the split of `bin`
##' at a random location on a random margin
##' @examples
##' bin <- list(x = 1:10, y = sample(1:10),
##'             bnds = list(x = c(0, 10), y = c(0, 10)),
##'             expn = 10, n = 10, depth = 0)
##' rUnifSplit(bin, minExp = 2)
##' @author Chris Salahub
rUnifSplit <- function (bin, minExp = 0, splitLongSide = FALSE, ...) {
    expn <- bin$expn
    prop <- minExp/expn
    xrng <- diff(bin$bnds$x)
    yrng <- diff(bin$bnds$y)
    xmax <- runif(1, min = bin$bnds$x[1] + prop*xrng,
                  max = bin$bnds$x[2] - prop*xrng)
    ymax <- runif(1, min = bin$bnds$y[1] + prop*yrng,
                  max = bin$bnds$y[2] - prop*yrng)
    wider <- (bin$bnds$x[2] - bin$bnds$x[1]) >
        (bin$bnds$y[2] - bin$bnds$y[1])
    u <- runif(1)
    if (u >= 0.5 | (splitLongSide & wider)) {
        newbnd <- xmax
        below <- which(bin$x <= xmax)
        above <- which(bin$x > xmax)
        splitX(bin, bd = newbnd, above = above, below = below)
    }
    else {
        newbnd <- ymax
        below <- which(bin$y <= ymax)
        above <- which(bin$y > ymax)
        splitY(bin, bd = newbnd, above = above, below = below)
    }
}

##' @title Univariate score maximizing splitting
##' @description A function which splits a bin based on the location
##' maximizing a score function.
##' @details This function is the univariate version of
##' `maxScoreSplit` and so is considerably simpler. It assumes the
##' variable to be split is named `x` in the bin, and the other
##' variable is to remain unsplit.
##' @param bin a bin to be split with elements `x`, `y`, `depth`,
##' `bnds` (list with elements `x` and `y`), `expn`, `n`
##' @param scorer function which accepts a numeric vector of potential
##' split coordinates and the bounds of `bin` and returns a numeric
##' vector of scores for each
##' @param pickMax function which accepts a list of scores and returns
##' the element of the largest score according to some rule
##' @param ... optional additional arguments to `scorer`
##' @return A list of two bins resulting from the split of `bin` at
##' the maximum split location along x
##' @author Chris Salahub
uniMaxScoreSplit <- function(bin, scorer = diff, pickMax = which.max,
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
