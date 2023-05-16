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
