## consider the general form of a binner for the data
binner <- function(x, y, stopper, splitter) {
  # initialize bin with all the data contained
  bins <- list(list(x = x, y = y, # x and y points
                    bnds = list(x = range(data$x, na.rm = TRUE),
                                y = range(data$y, na.rm = TRUE)),
                    expn = length(data$x), # default expectation is n
                    n = length(data$x), depth = 1)) # size of bin, depth
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
