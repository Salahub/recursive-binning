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
##' @return A list of two bins resulting from the split of `bin` in
##' half along the margin corresponding to the larger score.
##' @examples
##' bin <- list(x = 1:10, y = sample(1:10),
##'             bnds = list(x = c(0, 10), y = c(0, 10)),
##'             expn = 10, n = 10, depth = 0)
##' halfCutTie(bin, 1, 2) # splits on y
##' halfCutTie(bin, 2, 1) # splits on x
##' halfCutTie(bin, 1, 1) # ties are random
##' @author Chris Salahub
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
