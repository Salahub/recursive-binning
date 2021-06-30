## a small helper, extracts the indices of x in order
getSorted <- function(x, inds, ord = order(x), n = length(x)) {
    toGet <- rep(FALSE, n)
    toGet[inds] <- TRUE
    x[ord][toGet[ord]]
}

## simple possible criteria
stopCriterion <- function(depthLim = 4, sizeLim = 1, areaLim = 5) {
    function(depth, size, area) {
        list(Stop = depth >= depthLim | size <= sizeLim | area <= 2*areaLim)
    }
}

## scoring functions
chiScores <- function(n, sorted, bnds) {
    diffs <- diff(c(bnds[1], sorted, bnds[2]))
    total <- cumsum(diffs)
    h1 <- total[1:(n+1)] # length below
    h2 <- total[n+2] - h1 # length above
    d <- n/total[n+2] # density
    i <- 0:n # number below split i
    ni <- n - i # number above i
    (i - d*h1)^2/(h1*d) + (ni - d*h2)^2/(h2*d)
}
miScores <- function(n, sorted, bnds) {
    diffs <- diff(c(bnds[1], sorted, bnds[2]))
    total <- cumsum(diffs)
    h1 <- total[1:(n+1)] # length below
    h2 <- total[n+2] - h1 # length above
    d <- n/total[n+2] # density
    i <- 0:n # number below point i
    ni <- n - i # number above i
    below <- (i/n)*log(i/(d*h1))
    above <- (ni/n)*log(ni/(d*h2)) # split expectation
    below[1] <- 0
    above[n+1] <- 0 # handle known zeros
    below + above
}
difScores <- function(n, sorted, bnds) {
    diffs <- diff(c(bnds[1], sorted, bnds[2]))
    diffs/(diff(bnds))
}

## split to maximize a score
maxScoreSplit <- function(scores, ties = sample(c("x","y"), 1)) {
    function(inds, n, bnds, x, y, xord, yord, diffTol = 1.5e-8) {
        nbin <- length(inds) # size of bin to split
        xsort <- getSorted(x, inds, xord, n)
        ysort <- getSorted(y, inds, yord, n) # sorted order elements
        xsort <- c(xsort[1]-1, xsort)
        ysort <- c(ysort[1]-1, ysort) # split below bottom point
        xscore <- scores(nbin, xsort, bnds$x)
        yscore <- scores(nbin, ysort, bnds$y) # marginal split scores
        xscore[!is.finite(xscore)] <- 0
        yscore[!is.finite(yscore)] <- 0 # no null splits
        if (all(xscore == xscore[1])) { # perfect agreement
            if (all(yscore == yscore[1])) { # in both dims
                if (yscore[1] > xscore[1]) {
                    splitVar <- "y"
                } else if (xscore[1] > yscore[1]) {
                    splitVar <- "x"
                } else {
                    splitVar <- ties # follow tie convention
                }
                splitPos <- floor(nbin/2)+1 # set index
            } else {
                splitVar <- "y"
                splitPos <- which.max(yscore)
            }
        } else if (all(yscore == yscore[1])) {
            splitVar <- "x"
            splitPos <- which.max(xscore)
        } else { # split on the maximizing value
            ymax <- which.max(yscore)
            xmax <- which.max(xscore)
            if (xscore[xmax] > yscore[ymax]) {
                splitVar <- "x"
                splitPos <- xmax
            } else if (xscore[xmax] < yscore[ymax]) {
                splitVar <- "y"
                splitPos <- ymax
            } else {
                splitPos <- c(x = xmax, y = ymax)[ties]
                splitVar <- ties
            }
        }
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

## wrapper to make a tree
makeTree <- function(x, y, criterion, splitter, leaf) {
    xord <- order(x)
    yord <- order(y) # sorted data
    n <- length(x) # sample size
    ## recursive tree function
    branch <- function(inds, bounds, depth = 0, size = length(inds)) {
        check <- criterion(depth, size = length(inds),
                           area = prod(sapply(bounds, diff))/n)
        if (check$Stop) { # stop and apply leaf
            leaf(size = size, bounds = bounds, x = x[inds], y = y[inds],
                 depth = depth)
        } else if (length(inds) == 0) {
            leaf(size = 0, bounds = bounds, x = numeric(0), y = numeric(0),
                 depth = depth)
        } else { # split
            splitted <- splitter(inds, n, bounds, x, y,
                                 xord = xord, yord = yord)
            if (splitted$allZero) { # denegerate case
                leaf(size = length(inds),
                     bounds = bounds, x = x[inds], y = y[inds],
                     depth = depth)
            } else {
                names <- c(paste0(splitted$var, "<=", splitted$at),
                           paste0(splitted$var, ">", splitted$at))
                aboveBnds <- belowBnds <- bounds # fill new bounds
                aboveBnds[[splitted$var]][1] <- splitted$at
                belowBnds[[splitted$var]][2] <- splitted$at
                setNames(list(branch(splitted$below, bounds = belowBnds,
                                     depth = depth + 1,
                                     size = length(splitted$below)),
                              branch(splitted$above, bounds = aboveBnds,
                                     depth = depth + 1,
                                     size = length(splitted$above))),
                         names) # pretty names
            }
        }
    }
    branch(1:length(x), bounds = list(x = c(0, x[xord[n]]),
                                      y = c(0, y[yord[n]])))
}
