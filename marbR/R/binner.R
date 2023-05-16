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
