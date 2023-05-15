## consider the general form of a binner for the data
binner <- function(data, stopper, splitter) {
  # initialize bin with all the data contained
  binList <- list(list(x = data$x, y = data$y, # x and y points
                       bnds = list(x = range(data$x, na.rm = TRUE),
                                   y = range(data$y, na.rm = TRUE)), # boundaries
                       n = length(data$x), depth = 1)) # size of bin, depth
  stopStatus <- stopper(binList) # initialize logical vector

  while (any(!stopStatus)) { # check the stop criteria
    newBins <- binList[stopStatus] # stopped bins
    for (bin in binList[!stopStatus]) { # split all other bins
      newBins <- c(newBins, splitter(bin)) # split and add
    }
    binList <- newBins # update binList
    stopStatus <- stopper(binList) # check criteria
  }

  binList # return the final list of bins
}

## now define the splitter and stopper
## bin size stopper
sizeStop <- function(binList, lim = 100) {
  lens <- sapply(binList, function(bin) bin$n) # bin sizes
  lens < lim # check against limit
}

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
