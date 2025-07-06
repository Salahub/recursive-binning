##` Marginal splitters
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
    list(makeBin(x = bin$x[below], y = bin$y[below],
                 bnds = list(x = c(bin$bnds$x[1], bd),
                             y = bin$bnds$y),
                 expn = bin$expn*belowfac,
                 n = bin$n-length(above), depth = bin$depth + 1,
                 stopped = FALSE),
         makeBin(x = bin$x[above], y = bin$y[above],
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
    list(makeBin(x = bin$x[below], y = bin$y[below],
                 bnds = list(x = bin$bnds$x,
                             y = c(bin$bnds$y[1], bd)),
                 expn = bin$expn*belowfac,
                 n = bin$n-length(above), depth = bin$depth + 1,
                 stopped = FALSE),
         makeBin(x = bin$x[above], y = bin$y[above],
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
##' @param squarify logical value, should we force splitting on
##' the longer side regardless of scores?
##' @return A list of two bins resulting from the split of `bin` in
##' half along the margin corresponding to the larger score.
##' @examples
##' bin <- makeBin(x = 1:10, y = sample(1:10))
##' halfCutTie(bin, 1, 2, wider = FALSE) # splits on y
##' halfCutTie(bin, 2, 1, wider = FALSE) # splits on x
##' halfCutTie(bin, 1, 1, wider = FALSE) # ties are random
##' @author Chris Salahub
halfCutTie <- function(bin, xscore, yscore, wider,
                       squarify = FALSE) {
    u <- as.numeric(yscore < xscore) # prefer to split on max score
    if (yscore == xscore) u <- runif(1)
    if (squarify) u <- as.numeric(wider)
    if (u >= 0.5) { # x has a larger score or bin is wider than tall
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

##' @title Size-restricted bivariate score maximizing splitting
##' @description Splits a bin based on the location maximizing a score
##' function with restrictions on minimum bin size.
##' @details This function serves as a wrapper which manages the
##' logic of splitting bins using a score function while maintaining
##' a minimum size and allowing forced splits along the wider edge.
##' @param bin a bin to be split with elements `x`, `y`, `depth`,
##' `bnds` (list with elements `x` and `y`), `expn`, `n`
##' @param scorer function which accepts a numeric vector of potential
##' split coordinates and the bounds of `bin` and returns a numeric
##' vector of scores for each
##' @param minExp value giving the smallest expected count allowed for
##' bin splits
##' @param squarify logical value, should we force splitting on
##' the longer side regardless of scores?
##' @return A list of two bins resulting from the split of `bin`
##' along the corresponding margin at the maximum location
##' @examples
##' bin <- makeBin(x = 1:10, y = sample(1:10))
##' maxScoreSplit(bin, chiScores)
##' maxScoreSplit(bin, miScores) # pretty similar for both
##' maxScoreSplit(bin, randScores)
##' maxScoreSplit(bin, randScores) # different every time
##' @author Chris Salahub
maxScoreSplit <- function(bin, scorer, minExp = 5,
                          squarify = FALSE) {
    expn <- bin$expn
    prop <- minExp/expn
    deltax <- diff(bin$bnds$x)
    deltay <- diff(bin$bnds$y)
    xfrm <- bin$bnds$x + prop*deltax*c(1, -1) # minExp bounds
    yfrm <- bin$bnds$y + prop*deltay*c(1, -1)
    if ((xfrm[2] - xfrm[1]) >= 0) {
        ## take advantage of R ordering: ensure frame includes pts
        xaug <- c(min(bin$x)-1, bin$x, xfrm) # augment with frame
        xsort <- order(xaug)
        xlowFrm <- which(xsort == (length(xaug) - 1)) # identify frame
        xupFrm <- which(xsort == length(xaug))
        xinfrm <- xlowFrm:xupFrm
        xcntblw <- cumsum(c(0, rep(1, bin$n), 0, 0)[xsort])[xinfrm]
        xscore <- scorer(bounds = c(bin$bnds$x[1],
                                    xaug[xsort][xinfrm],
                                    bin$bnds$x[2]),
                         nbelow = xcntblw, n = bin$n)
        xmax <- which.max(xscore)
        if (all(abs(xscore - xscore[1]) < sqrt(.Machine$double.eps))) {
            xmax <- ceiling((xupFrm - xlowFrm)/2)
        }  # halve on ties
    } else xmax <- NA # indicate inability to split
    if ((yfrm[2] - yfrm[1]) >= 0) {
        yaug <- c(min(bin$y)-1, bin$y, yfrm)
        ysort <- order(yaug)
        ylowFrm <- which(ysort == (length(yaug) - 1))
        yupFrm <- which(ysort == length(yaug))
        yinfrm <- ylowFrm:yupFrm
        ycntblw <- cumsum(c(0, rep(1, bin$n), 0, 0)[ysort])[yinfrm]
        yscore <- scorer(bounds = c(bin$bnds$y[1],
                                    yaug[ysort][yinfrm],
                                    bin$bnds$y[2]),
                         nbelow = ycntblw, n = bin$n)
        ymax <- which.max(yscore)
        if (all(abs(yscore - yscore[1]) < sqrt(.Machine$double.eps))) {
            ymax <- ceiling((yupFrm - ylowFrm)/2)
        }
    } else ymax <- NA
    ## control the bin being split
    if (is.na(xmax) & is.na(ymax)) {
        swtch <- NA
    } else {
        if (is.na(xmax)) {
            swtch <- 1
        } else if (is.na(ymax)) {
            swtch <- 0
        } else {
            wider <- deltax >= deltay
            if (yscore[ymax] > xscore[xmax]) {
                swtch <- 1
            } else if (yscore[ymax] == xscore[xmax]) {
                swtch <- runif(1)
            } else swtch <- 0
            if (squarify) {
                if (wider) swtch <- 0 else swtch <- 1
            }
        }   
    }
    ## use swtch variable to control behaviour
    if (is.na(swtch)) {
        bin$stopped <- TRUE
        list(bin)
    } else if (swtch > 0.5) {
        newbnd <- yaug[ysort][yinfrm[ymax]]
        ## remove frame, remove new lowest min-1, adjust indices 
        yclean <- ysort[-c(yupFrm, ylowFrm)][-1]-1 # remove extras
        below <- yclean[seq_len(ycntblw[ymax])]
        if (newbnd >= bin$y[yclean[bin$n]]) { # the complement of below
            above <- integer(0)
        } else above <- yclean[(ycntblw[ymax]+1):bin$n]
        splitY(bin, bd = newbnd, above = above, below = below)
    } else {
        newbnd <- xaug[xsort][xinfrm[xmax]]
        xclean <- xsort[-c(xupFrm, xlowFrm)][-1] - 1
        below <- xclean[seq_len(xcntblw[xmax])]
        if (newbnd >= bin$x[xclean[bin$n]]) {
            above <- integer(0)
        } else above <- xclean[(xcntblw[xmax]+1):bin$n]
        splitX(bin, bd = newbnd, above = above, below = below)
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
##' @param minExp numeric giving the minimum expected count allowed
##' in a bin
##' @param on one of "x" or "y": the margin to split
##' @return A list of two bins resulting from the split of `bin` at
##' the maximum split location along y
##' @author Chris Salahub
uniMaxScoreSplit <- function(bin, scorer, minExp = 5, on = "y") {
    expn <- bin$expn
    prop <- minExp/expn
    if (on == "y") {
        deltay <- diff(bin$bnds$y)
        yfrm <- bin$bnds$y + ceiling(prop*deltay)*c(1, -1)
        if ((yfrm[2] - yfrm[1]) >= 0) {
            yaug <- c(min(bin$y)-1, bin$y, yfrm)
            ysort <- order(yaug)
            ylowFrm <- which(ysort == (length(yaug) - 1))
            yupFrm <- which(ysort == length(yaug))
            yinfrm <- ylowFrm:yupFrm
            ycntblw <- cumsum(c(0, rep(1, bin$n), 0, 0)[ysort])[yinfrm]
            yscore <- scorer(bounds = c(bin$bnds$y[1],
                                        yaug[ysort][yinfrm],
                                        bin$bnds$y[2]),
                             nbelow = ycntblw, n = bin$n)
            ymax <- which.max(yscore)
            if (all(abs(yscore - yscore[1]) < sqrt(.Machine$double.eps))) {
                ymax <- ceiling((yupFrm - ylowFrm)/2)
            }
        } else ymax <- NA
        if (is.na(ymax)) {
            bin$stopped <- TRUE
            list(bin)
        } else {
            newbnd <- yaug[ysort][yinfrm[ymax]]
            yclean <- ysort[-c(yupFrm, ylowFrm)][-1] - 1
            below <- yclean[seq_len(ycntblw[ymax])]
            if (newbnd >= bin$y[yclean[bin$n]]) {
                above <- integer(0)
            } else above <- yclean[(ycntblw[ymax]+1):bin$n]
            splitY(bin, bd = newbnd, above = above, below = below)
        }
    } else if (on == "x") {
        deltax <- diff(bin$bnds$x)
        xfrm <- bin$bnds$x + ceiling(prop*deltax)*c(1, -1)
        if ((xfrm[2] - xfrm[1]) >= 0) {
            xaug <- c(min(bin$x)-1, bin$x, xfrm)
            xsort <- order(xaug)
            xlowFrm <- which(xsort == (length(xaug) - 1))
            xupFrm <- which(xsort == length(xaug))
            xinfrm <- xlowFrm:xupFrm
            xcntblw <- cumsum(c(0, rep(1, bin$n), 0, 0)[xsort])[xinfrm]
            xscore <- scorer(bounds = c(bin$bnds$x[1],
                                        xaug[xsort][xinfrm],
                                        bin$bnds$x[2]),
                             nbelow = xcntblw, n = bin$n)
            xmax <- which.max(xscore)
            if (all(abs(xscore - xscore[1]) < sqrt(.Machine$double.eps))) {
                xmax <- ceiling((xupFrm - xlowFrm)/2)
            }
        } else xmax <- NA
        if (is.na(xmax)) {
            bin$stopped <- TRUE
            list(bin)
        } else {
            newbnd <- xaug[xsort][xinfrm[xmax]]
            xclean <- xsort[-c(xupFrm, xlowFrm)][-1] - 1
            below <- xclean[seq_len(xcntblw[xmax])]
            if (newbnd >= bin$x[xclean[bin$n]]) {
                above <- integer(0)
            } else above <- xclean[(xcntblw[xmax]+1):bin$n]
            splitY(bin, bd = newbnd, above = above, below = below)
        }
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
##' @param minExp value giving the smallest expected count allowed for
##' bin splits
##' @param pickMax function which accepts a list of scores and returns
##' the element of the largest score according to some rule
##' @param ... optional additional arguments to `scorer`
##' @return A list of two bins resulting from the split of `bin`
##' along the corresponding margin at the maximum location
##' @author Chris Salahub
sandboxMaxSplit <- function(bin, scorer, ties = halfCutTie,
                            minExp = 5, pickMax = which.max, ...) {
    xsort <- order(bin$x)
    ysort <- order(bin$y) # get marginal ordering
    xscore <- scorer(bin, xsort, ...)
    yscore <- scorer(bin, ysort, ...) # score splits based on ordering
    xmax <- pickMax(xscore)
    ymax <- pickMax(yscore) # the score values
    xallEq <- all(abs(xscore - xscore[1]) < sqrt(.Machine$double.eps))
    yallEq <- all(abs(yscore - yscore[1]) < sqrt(.Machine$double.eps))
    u <- as.numeric(xscore[xmax] >= yscore[ymax]) # ties go to x
    if (xallEq & yallEq) { # in the case of ties, use tie function
        ties(bin, xscore[1], yscore[1])
    } else if (u >= 0.5) {
        xsplts <- bin$x[xsort]
        newbnd <- c(xsplts[1]-1, xsplts)[xmax] # new boundary
        below <- xsort[seq_len(xmax-1)] # get indices of points below
        above <- if (xmax == bin$n+1) integer(0) else xsort[xmax:bin$n]
        splitX(bin, bd = newbnd, above = above, below = below)
    } else {
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
##' @param squarify logical value, should we force splitting on
##' the longer side?
##' @return A list of two bins resulting from the split of `bin`
##' along the corresponding margin at the maximum location
##' @examples
##' bin <- makeBin(x = 1:10, y = sample(1:10))
##' rIntSplit(bin, minExp = 2)
##' @author Chris Salahub
rIntSplit <- function(bin, minExp = 5, squarify = TRUE) {
  expn <- bin$expn
  prop <- minExp/expn
  deltax <- diff(bin$bnds$x)
  deltay <- diff(bin$bnds$y)  
  wider <- deltax > deltay
  if (squarify) u <- as.numeric(wider) else u <- runif(1)
  if (u >= 0.5) { # split on x
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

##' @title Univariate random integer splitting
##' @description A function which splits a bin along x at a random
##' integer conforming to limits on minimum bin size.
##' @details This function serves as a wrapper which manages the
##' interaction of a score function, marginal splitting functions,
##' tie breaking function, and a maximum selection function to split
##' a bin along a single margin at the observation coordinate which
##' maximizes the score function.
##' @param bin a bin to be split with elements `x`, `y`, `depth`,
##' `bnds` (list with elements `x` and `y`), `expn`, `n`
##' @param minExp numeric giving the minimum expected count allowed
##' in a bin
##' @param on one of "x" or "y": the margin to split
##' @return A list of two bins resulting from the split of `bin`
##' along the corresponding margin at the maximum location
##' @examples
##' bin <- makeBin(x = 1:10, y = sample(1:10))
##' rIntSplit(bin, minExp = 2)
##' @author Chris Salahub
uniRIntSplit <- function(bin, minExp = 5, on = "y") {
  expn <- bin$expn
  prop <- minExp/expn
  if (on == "y") {
      deltay <- diff(bin$bnds$y)
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
  } else if (on == "x") {
      deltax <- diff(bin$bnds$x)
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
      stop('Argument "on" must be one of "x" or "y"')
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
##' @param squarify logical value, should we force splitting on
##' the longer side?
##' @return A list of two bins resulting from the split of `bin`
##' at a random location on a random margin
##' @examples
##' bin <- makeBin(x = 1:10, y = sample(1:10))
##' rUnifSplit(bin, minExp = 2)
##' @author Chris Salahub
rUnifSplit <- function (bin, minExp = 0, squarify = FALSE) {
    expn <- bin$expn
    prop <- minExp/expn
    deltax <- diff(bin$bnds$x)
    deltay <- diff(bin$bnds$y)  
    wider <- deltax > deltay
    if (squarify) u <- as.numeric(wider) else u <- runif(1)
    if (u >= 0.5) {
        lower <- bin$bnds$x[1] + prop*deltax
        upper <- bin$bnds$x[2] - prop*deltax
        if (upper <= lower) {
            bin$stopped <- TRUE
            list(bin)
        } else{
            spos <- runif(1, min = lower, max = upper)
            above <- which(bin$x > spos)
            below <- which(bin$x <= spos)
            splitX(bin, spos, above, below)
        }
    }
    else {
        lower <- bin$bnds$y[1] + prop*deltay
        upper <- bin$bnds$y[2] - prop*deltay
        if (upper <= lower) {
            bin$stopped <- TRUE
            list(bin)
        } else{
            spos <- runif(1, min = lower, max = upper)
            above <- which(bin$y > spos)
            below <- which(bin$y <= spos)
            splitY(bin, spos, above, below)
        }
    }
}
