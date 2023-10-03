## PACKAGES ##########################################################
## recursive binning package
library(marbR)



## CONSTANTS/RUN PARAMETERS ##########################################
writeout <- FALSE # should simulations be run and output written?
depthPal <- hcl.colors(9, "Dark 2") # colours by depth



## FUNCTIONS #########################################################
## custom plotting function with narrow margins
narrowPlot <- function(xgrid, ygrid, main = "", xlab = "", ylab = "",
                       xticks = xgrid, yticks = ygrid,
                       mars = c(2.1, 2.1, 1.1, 1.1),
                       xlim = range(xgrid), ylim = range(ygrid),
                       addGrid = TRUE, ...) {
    par(mar = mars) # set narrow margins
    plot(NA, ylim = ylim, xlim = xlim, xaxt = 'n', xlab = "",
         yaxt = 'n', ylab = "", main = "", ...)
    ## add labels
    mtext(main, side = 3, line = 0, cex = 0.8) # main
    mtext(ylab, side = 2, line = 1, cex = 0.8) # ylab
    mtext(xlab, side = 1, line = 1, padj = 0, cex = 0.8) # xlab
    ## add grid lines
    if (addGrid) {
        abline(h = ygrid, v = xgrid, lty = 1,
               col = adjustcolor("gray", alpha.f = 0.4))
    }
    ## and ticks
    mtext(side = 1, at = xgrid, text = "|", line = 0, cex = 0.5,
          padj = -2)
    mtext(text = xticks, at = xgrid, side = 1, cex = 0.8)
    mtext(side = 2, at = ygrid, text = "|", line = 0, cex = 0.5,
          padj = 1)
    mtext(text = yticks, at = ygrid, side = 2, cex = 0.8)
}

## reduce size of binning in storage by dropping points
dropBinPoints <- function(bins) {
    lapply(bins, function(bn) {
        bn$x <- NULL; bn$y <- NULL; bn
    })
}



## SIMULATED DATA PATTERNS ###########################################
## patterns from Newton (2009) provided in a list of functions
patFns <- list(
    wave = function(n) {
        x <- seq(-1, 1, length=n)
        u <- x + runif(n)/3; v <- 4*((x^2 - 1/2)^2 + runif(n)/500)
        cbind(x = u, y = v)
    },
    rotatedSquare = function(n) {
        x <- runif(n, min = -1, max = 1)
        y <- runif(n, min = -1, max = 1)
        theta <--pi/8
        rr <- rbind(c(cos(theta), -sin(theta)),
                    c(sin(theta), cos(theta)))
        tmp <- cbind(x, y) %*% rr
        colnames(tmp) <- c("x",  "y")
        tmp
    },
    circle = function(n) {
        x <- runif(n, min = -1, max = 1)
        y <- runif(n, min = -1, max = 1)
        theta <- -pi/4
        rr <- rbind(c(cos(theta), -sin(theta)),
                    c(sin(theta), cos(theta)))
        tmp <- cbind(x, y) %*% rr
        colnames(tmp) <- c("x",  "y")
        tmp
    },
    valley = function(n) {
        x <- seq(-1,1, length=n )
        y <- (x ^2 + runif(n))/2
        cbind(x = x, y = y)
    },
    cross = function(n) {
        x <- seq(-1, 1, length = n)
        y <- (x^2 + runif(n)/2)*(sample(c(-1,1), size=n, replace = T))
        cbind(x = x, y = y)
    },
    ring = function(n) {
        x <- seq(-1, 1, length = n)
        u <- sin(x*pi) + rnorm(n)/8
        v <- cos(x*pi) + rnorm(n)/8
        cbind(x = u, y = v)
    },
    noise = function(n) {
        dx <- rnorm(n)/3
        dy <- rnorm(n)/3
        cx <- sample(c(-1, 1), size=n, replace = T)
        cy <- sample(c(-1, 1), size=n, replace = T)
        u <- cx + dx
        v <- cy + dy
        cbind(x = u, y = v)
    })

## write a wrapper for these patterns to generate an array of all
generatePatterns <- function(n) {
    simplify2array(lapply(patFns, function(fn) fn(n)))
}

## generate many repetitions of each to bin
set.seed(70111238)
n <- 1000
nsim <- 100
simData <- replicate(nsim, generatePatterns(n))

## plot the first data realization
m <- 1
pal <- c(RColorBrewer::brewer.pal(6, "Set2"), "black")
par(mfrow=c(1,7), mar=c(1,1,1,1)/2)
for(i in 1:7) {
    plot(simData[, "x", i, 1], simData[, "y", i, 1], xlab = "",
         ylab = "", axes = F, pch = 19, cex = 0.2, col = pal[i])
}

## convert this data into pairwise ranks and plot it
simXr <- apply(simData[, "x", , ], c(2, 3), rank)
simYr <- apply(simData[, "y", , ], c(2, 3), rank)
for(i in 1:7) {
     plot(simXr[, i, 1], simYr[, i, 1], xlab = "", ylab = "",
          axes= F, pch = 19, cex = 0.2, col = pal[i])
}

## try the binning algorithm on these data
## define a range of depths
depths <- 1:10
## define the criteria dynamically, works due to R's lexical scoping
crits <- makeCriteria(depth >= ii, expn <= 10, n == 0)
## define the stop function
stopFn <- function(bns) stopper(bns, crits)
## and splitting functions
chiSplit <- function(bn) maxScoreSplit(bn, chiScores, minExp = 5)
miSplit <- function(bn) maxScoreSplit(bn, miScores, minExp = 5)
rndSplit <- function(bn) maxScoreSplit(bn, randScores, minExp = 5)
## allocate storage for every split method
testChiBins <- vector("list", nsim)
testMiBins <- vector("list", nsim)
testRndBins <- vector("list", nsim)

## bin each realization (takes a few minutes)
for (jj in 1:nsim) {
    ## each list element is also a list for each
    testChiBins[[jj]] <- vector("list", length(depths))
    testMiBins[[jj]] <- vector("list", length(depths))
    testRndBins[[jj]] <- vector("list", length(depths))
    for (ii in seq_along(depths)) { # iterate through depths
        ## chi bins for each pattern
        testChiBins[[jj]][[ii]] <- lapply(1:7, function(kk) {
            binner(simXr[, kk, jj], simYr[, kk, jj],
                   stopper = stopFn, splitter = chiSplit)
        })
        ## mi bins for each pattern
        testMiBins[[jj]][[ii]] <- lapply(1:7, function(kk) {
            binner(simXr[, kk, jj], simYr[, kk, jj],
                   stopper = stopFn, splitter = miSplit)
        })
        ## finally, random bins for each pattern
        testRndBins[[jj]][[ii]] <- lapply(1:7, function(kk) {
            binner(simXr[, kk, jj], simYr[, kk, jj],
                   stopper = stopFn, splitter = rndSplit)
        })
    }
}

## compute the chi square statistics for each split method
testChiChi <- lapply(testChiBins, # nested list  makes it ugly
                     function(lst) {
                         lapply(lst,
                                function(el) lapply(el, binChi))
                     })
testMiChi <- lapply(testMiBins, ## same thing for mi...
                     function(lst) {
                         lapply(lst,
                                function(el) lapply(el, binChi))
                     })
testRndChi <- lapply(testRndBins, ## ... and random splitting
                     function(lst) {
                         lapply(lst,
                                function(el) lapply(el, binChi))
                     })

## for ease of plotting, convert these tests to statistic values,
## final bin counts
## define some helpers to make this cleaner...
## wrapper to apply function to a nested list and return an array
deNest <- function(nstdLst, fn) {
    lapply(nstdLst, function(olst) {
        sapply(olst, function(lst) {
            sapply(lst, fn)
        })
    })
}
## the internal functions to work with deNest
getStat <- function(x) x$stat
getnBin <- function(x) length(x$residuals)
getMaxRes <- function(x) max(abs(x$residuals))

## apply this to everything else
chiPaths <- deNest(testChiChi, getStat)
chiNbin <- deNest(testChiChi, getnBin)
miPaths <- deNest(testMiChi, getStat)
miNbin <- deNest(testMiChi, getnBin)
rndPaths <- deNest(testRndChi, getStat)
rndNbin <- deNest(testRndChi, getnBin)

## plot the paths of every pattern under different splitting regimes
## compared to the null
data(null1000) # read in null
data <- get("null1000")
depths <- data$depths
depthSeq.chi <- data$chiSplit
depthSeq.mi <- data$miSplit
depthSeq.rnd <- data$randSplit # null data

par(mfrow = c(1,1))
## plot paths for an individual random split
narrowPlot(xgrid = seq(0, 160, by = 40), xlab = "Number of bins",
           ygrid = seq(0, 1600, by = 400),
           ylab = expression(chi^2~statistic))
for (ii in 1:1e4) { # add the null lines
    lines(depthSeq.rnd["nbin",,ii],
          depthSeq.rnd["chi",,ii],
          col = adjustcolor("gray", 0.1))
}
## add these as lines to the plot of null lines
for (jj in 1:7) {
    lines(rndNbin[[1]][jj,], rndPaths[[1]][jj,], col = pal[jj])
    points(rndNbin[[1]][jj,], rndPaths[[1]][jj,], col = pal[jj],
           pch = 19, cex = 0.5)
}
## add the 95% chi quantile
lines(1:160, qchisq(0.95, 1:160), lty = 2)

## make the same plot for paths from chi splitting
narrowPlot(xgrid = seq(0, 160, by = 40),
           xlab = "Number of bins",
           ygrid = seq(0, 1600, by = 400),
           ylab = expression(chi^2~statistic))
for (ii in 1:1e4) {
    lines(depthSeq.chi["nbin",,ii],
          depthSeq.chi["chi",,ii],
          col = adjustcolor("gray", 0.1))
}
for (jj in 1:7) {
    lines(chiNbin[[1]][jj,], chiPaths[[1]][jj,], col = pal[jj])
    points(chiNbin[[1]][jj,], chiPaths[[1]][jj,], col = pal[jj],
           pch = 19, cex = 0.5)
}

## for the random split repetitions, plot every one
narrowPlot(xgrid = seq(0, 160, by = 40),
           xlab = "Number of bins",
           ygrid = seq(0, 1200, by = 300), ylim = c(0, 1300),
           ylab = expression(chi^2~statistic))
for (ii in 1:1e4) {
    lines(depthSeq.rnd["nbin",,ii],
          depthSeq.rnd["chi",,ii],
          col = adjustcolor("gray", 0.1))
}
for (jj in 1:7) {
    for (ii in 1:100) {
        lines(rndNbin[[ii]][jj,], rndPaths[[ii]][jj,],
              col = adjustcolor(pal[jj], 0.2))
    }
}
lines(1:160, qchisq(0.95, 1:160), lty = 2)

## do the same for the chi splits
narrowPlot(xgrid = seq(0, 160, by = 40),
           xlab = "Number of bins",
           ygrid = seq(0, 1200, by = 300), ylim = c(0, 1300),
           ylab = expression(chi^2~statistic))
for (ii in 1:1e4) {
    lines(depthSeq.chi["nbin",,ii],
          depthSeq.chi["chi",,ii],
          col = adjustcolor("gray", 0.1))
}
for (jj in 1:7) {
    for (ii in 1:100) {
        lines(chiNbin[[ii]][jj,], chiPaths[[ii]][jj,],
              col = adjustcolor(pal[jj], 0.2))
    }
}
lines(1:160, qchisq(0.95, 1:160), lty = 2)

## next, check the bins for every depth
## start by getting the maximum residual to make the shading constant
maxRes <- max(sapply(unlist(testChiChi[[1]],
                            recursive = FALSE),
                     getMaxRes))
## for every depth, display the binning for each pattern
for (depth in 2:10) {
    par(mfrow=c(1,7), mar=c(1,1,1,1)/2)
    for(i in 1:7) {
        plot(NA, ylim = c(1, n), xlim = c(1, n), # remove axes
             axes = F, xlab = "", ylab = "", main = "")
        plotBinning(testChiBins[[1]][[depth]][[i]], pch = 19,
                    cex = 0.1, add = TRUE,
                    col = adjustcolor("grey", 0.8),
                    fill = residualFill(testChiBins[[1]][[depth]][[i]],
                                        maxRes = maxRes))
    }
}

## do the same thing for the bins from MI splitting
maxRes <- max(sapply(unlist(testMiChi[[1]],
                            recursive = FALSE),
                     getMaxRes))
for (depth in 2:10) {
    par(mfrow=c(1,7), mar=c(1,1,1,1)/2)
    for(i in 1:7) {
        plot(NA, ylim = c(1, n), xlim = c(1, n),
             axes = F, xlab = "", ylab = "", main = "")
        plotBinning(testMiBins[[1]][[depth]][[i]], pch = 19,
                    cex = 0.1, add = TRUE,
                    col = adjustcolor("grey", 0.8),
                    fill = residualFill(testMiBins[[1]][[depth]][[i]],
                                        maxRes = maxRes))
    }
}

## finally, repeat this for the random splitting
maxRes <- max(sapply(unlist(testRndChi[[10]],
                            recursive = FALSE),
                     getMaxRes))
for (depth in 2:10) {
    par(mfrow=c(1,7), mar=c(1,1,1,1)/2)
    for(i in 1:7) {
        plot(NA, ylim = c(1, n), xlim = c(1, n),
             axes = F, xlab = "", ylab = "", main = "")
        plotBinning(testRndBins[[10]][[depth]][[i]],
                    pch = 19, cex = 0.1, add = TRUE,
                    col = adjustcolor("grey", 0.8),
                    fill = residualFill(testRndBins[[10]][[depth]][[i]],
                                        maxRes = maxRes))
    }
}
