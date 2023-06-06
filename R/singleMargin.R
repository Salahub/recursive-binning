## functionality that re-thinks previous work: write a function to
## to split on a given margin, one to choose a margin, and only
## univariate logic
checkFn <- function(depth, size, area) {


## start with the simplest case: a univariate recursive function
univarRecurse <- function(x, criterion, splitter, leaf) {
    xsort <- sort(x) # sort the data
    n <- length(x) # sample size
    ## the recursive tree function
    branch <- function(inds, bounds, depth = 0, size = length(inds)) {
        check <- criterion(depth, size, # check stop
                           area = prod(sapply(bounds, diff))/2)
        if (check) { # apply leaf if stopping
            leaf(size = size, bounds = bounds, x = x[inds], depth = depth)
        } else if (size == 0) { # hard limit regardless of criterion
            leaf(size = 0, bounds = bounds, x = x[inds], depth = depth)
        } else {
            splt <- splitter(inds, bounds, x) # split the bin
            nms <- with(splt, c(paste0(var, "<=", at),
                                paste0(var, ">", at))) # nice names
            bnds.ab <- bnds.be <- bounds # new bin bounds
            bnds.ab[1] <- splt$at
            bnds.be[2] <- splt$at # fill approriate values
            setNames(list(branch(splt$below, bounds = bnds.be,
                                 depth = depth + 1,
                                 size = length(splt$below)),
                          branch(splt$above, bounds = bnds.ab,
                                 depth = depth + 1,
                                 size = length(splt$above))),
                     nms) # recursive step
        }
    }
    branch(1:n, bounds = c(0, xsort[n])) # initialize
}

## the bivariate recursion needs to choose a margin to split on
bivarRecurse <- function(x, criterion, splitter, leaf) {
    xsort <- sort(x) # sort the data
    n <- length(x) # sample size
    ## the recursive tree function
    branch <- function(inds, bounds, depth = 0, size = length(inds)) {
        check <- criterion(depth, size, # check stop
                           area = prod(sapply(bounds, diff))/2)
        if (check) { # apply leaf if stopping
            leaf(size = size, bounds = bounds, x = x[inds], depth = depth)
        } else if (size == 0) { # hard limit regardless of criterion
            leaf(size = 0, bounds = bounds, x = x[inds], depth = depth)
        } else {
            splt <- splitter(inds, bounds, x) # split the bin
            nms <- with(splt, c(paste0(var, "<=", at),
                                paste0(var, ">", at))) # nice names
            bnds.ab <- bnds.be <- bounds # new bin bounds
            bnds.ab[1] <- splt$at
            bnds.be[2] <- splt$at # fill approriate values
            setNames(list(branch(splt$below, bounds = bnds.be,
                                 depth = depth + 1,
                                 size = length(splt$below)),
                          branch(splt$above, bounds = bnds.ab,
                                 depth = depth + 1,
                                 size = length(splt$above))),
                     nms) # recursive step
        }
    }
    branch(1:n, bounds = c(0, xsort[n])) # initialize
}

