## max gap splitter
maxGapSplit <- function(bin) {
  xsort <- order(bin$x)
  ysort <- order(bin$y) # get marginal ordering
  xdiffs <- diff(bin$x[xsort])
  ydiffs <- diff(bin$y[ysort]) # get differences
  xmax <- which.max(xdiffs)
  ymax <- which.max(ydiffs) # the maximum differences
  if (xdiffs[xmax] >= ydiffs[ymax]) { # ties go to x
    newbnd <- mean(bin$x[xsort][xmax:(xmax+1)]) # new bin boundary
    above <- xsort[(xmax+1):(bin$n)] # get indices of points above
    list(list(x = bin$x[-above], y = bin$y[-above],
              bnds = list(x = c(bin$bnds$x[1], newbnd),
                          y = bin$bnds$y),
              n = bin$n-length(above), depth = bin$depth + 1),
         list(x = bin$x[above], y = bin$y[above],
              bnds = list(x = c(newbnd, bin$bnds$x[2]),
                          y = bin$bnds$y),
              n = length(above), depth = bin$depth + 1)) # split bins
  } else {
    newbnd <- mean(bin$y[ysort][ymax:(ymax+1)]) # new bin boundary
    above <- ysort[(ymax+1):(bin$n)] # get indices of points above
    list(list(x = bin$x[-above], y = bin$y[-above],
              bnds = list(x = bin$bnds$x,
                          y = c(bin$bnds$y[1], newbnd)),
              n = bin$n - length(above), depth = bin$depth + 1),
         list(x = bin$x[above], y = bin$y[above],
              bnds = list(x = bin$bnds$x,
                          y = c(newbnd, bin$bnds$y[2])),
              n = length(above), depth = bin$depth + 1)) # split bins
  }
}

## score-maximizing split
maxScoreSplit <- function(scorefn, ties = "random", eTol = 10) {
    if (!(ties %in% c("random", "x", "y"))) { # check ties
        stop("ties must be one of 'random', 'x', or 'y'")
    }
    ## closure which splits based on scorefn
    function(inds, n, bnds, x, y, xord, yord, diffTol = 1.5e-8) {
        nbin <- length(inds) # size of the bin to split
        xsort <- getSorted(x, inds, xord, n)
        ysort <- getSorted(y, inds, yord, n) # elements in sorted order
        xsort <- c(xsort[1] - 1, xsort)
        ysort <- c(ysort[1] - 1, ysort) # split below bottom
        xdiff <- diff(bnds$x[1], xsort, bnds$x[2])
        ydiff <- diff(bnds$y[1], ysort, bnds$y[2]) # lengths
        ## split variable
        if (ties == "random") {
            splitVar <- sample(c("x", "y"), size = 1)
        } else {
            splitVar <- ties
        } # overwritten to split on max later
        ## splitting logic (get maximum)
        xmax <- getMax(nbin, xdiffs, scorefn, eTol)
        ymax <- getMax(nbin, ydiffs, scorefn, eTol) # index and score
        if (xmax[2] > ymax[2]) { # x maximizes
            splitVar <- "x"
            if (is.na(xmax[1])) {
                splitPos <- floor(nbin/2) + 1 # halve for no max
            } else {
                splitPos <- xmax[1] # take max index
            }
        } else if (ymax[2] > xmax[2]) {
            splitVar <- "y"
            if (is.na(ymax[1])) {
                splitPos <- floor(nbin/2) + 1
            } else {
                splitPos <- xmax[1]
            }
        } else {
            splitPos <- floor(nbin/2) + 1
        }
        ## apply the splits
        if (splitVar == "x") {
            splitVal <- xsort[splitPos]
            splitInds <- x[inds] <= splitVal
        } else {
            splitVal <- ysort[splitPos]
            splitInds <- y[inds] <= splitVal
        }
        allZero <- if (nbin >= 1) { # degenerate case
                       all(c(xsort[2:(nbin+1)]-xsort[2],
                             ysort[2:(nbin+1)]-ysort[2]) < diffTol)
                   } else TRUE
        list(var = splitVar, at = splitVal, allZero = allZero,
             below = inds[splitInds], above = inds[!splitInds])
    }
}
