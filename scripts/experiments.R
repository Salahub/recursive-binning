## PACKAGES ##########################################################
## recursive binning package
library(AssocBin)



## CONSTANTS/RUN PARAMETERS ##########################################
writeout <- FALSE # should simulations be run and output written?
depthPal <- hcl.colors(9, "Dark 2") # colours by depth



## FUNCTIONS #########################################################
## add marginal histograms to a scatterplot
addMarHists <- function(x, y, xcuts, ycuts) {
    bds <- par()$usr
    rowDist <- table(cut(x, xcuts))
    colDist <- table(cut(y, ycuts)) # marginal distributions
    vboxBds <- c(bds[2], bds[2] + 0.1*(bds[2] - bds[1]), bds[3:4])
    hboxBds <- c(bds[1:2], bds[4], bds[4] + 0.1*(bds[4] - bds[3]))
    ## density boxes
    rect(vboxBds[1], vboxBds[3], vboxBds[2], vboxBds[4], xpd = NA)
    rect(hboxBds[1], hboxBds[3], hboxBds[2], hboxBds[4], xpd = NA)
    ## add marginal histograms
    vseq <- ycuts
    rect(vboxBds[1], vseq[1:length(colDist)],
         vboxBds[1] + 0.9*diff(vboxBds[1:2])*(colDist/max(colDist)),
         vseq[2:(length(colDist) + 1)], xpd = NA,
         col = adjustcolor("firebrick", 0.5))
    hseq <- xcuts
    rect(hseq[1:length(rowDist)], hboxBds[3],
         hseq[2:(length(rowDist) + 1)], xpd = NA,
         hboxBds[3] + 0.9*diff(hboxBds[3:4])*(rowDist/max(rowDist)),
         col = adjustcolor("firebrick", 0.5))
}

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



## SIMPLE EXAMPLES ###################################################
## generate some random data
set.seed(16062021)
randx <- sample(1:1e3)
randy <- sample(1:1e3) # random ranks
## plots the progression of the algorithm for several different depths
dep <- c(1, 2, 10)
## R's lexical scoping allows for dynamic criteria construction: d is
## not defined, but will be within a later loop
criteria <- makeCriteria(expn <= 10, n == 0, depth >= d)
## define a stopper based on this
stopFn <- function(bns) stopper(bns, criteria)
## and some splitters for different score functions
chiSplit <- function(bn) maxScoreSplit(bn, chiScores, minExp = 5)
miSplit <- function(bn) maxScoreSplit(bn, miScores, minExp = 5)
rndSplit <- function(bn) maxScoreSplit(bn, randScores, minExp = 5)

## random data plot (Fig 4.1(a))
png("randomData.png", width = 2, height = 2, units = "in", res = 480)
narrowPlot(xgrid = seq(0, 1000, by = 500),
           ygrid = seq(0, 1000, by = 500),
           addGrid = FALSE, xlab = "x", ylab = "y")
points(randx, randy, pch = 19, cex = 0.25,
       col = adjustcolor("gray50", 0.5))
dev.off()

## allocate storage
randBin.chi <- vector(mode = "list", length = length(dep))
## split random data using different rules and depths
for (ii in seq_along(dep)) {
    d <- dep[ii]
    randBin.chi[[ii]] <- binner(randx, randy, stopper = stopFn,
                                splitter = chiSplit)
}

## the early progression of splitting (Figs 4.2(a), 4.3(a))
png("randomSplit1.png", width = 2, height = 2, units = "in", res = 480)
narrowPlot(xgrid = seq(0, 1000, by = 500),
           ygrid = seq(0, 1000, by = 500),
           addGrid = FALSE, xlab = "x", ylab = "y")
plotBinning(randBin.chi[[1]], pch = 19, cex = 0.25, add = TRUE,
            col = adjustcolor("gray50", 0.5))
dev.off()
png("randomSplit2.png", width = 2, height = 2, units = "in", res = 480)
narrowPlot(xgrid = seq(0, 1000, by = 500),
           ygrid = seq(0, 1000, by = 500),
           addGrid = FALSE, xlab = "x", ylab = "y")
plotBinning(randBin.chi[[2]], pch = 19, cex = 0.25, add = TRUE,
            col = adjustcolor("gray50", 0.5))
dev.off()

## plot bins coloured by depth (Fig 4.4(a))
png("randomSplitDepth.png", width = 2, height = 2, units = "in",
    res = 480)
narrowPlot(xgrid = seq(0, 1000, by = 500),
           ygrid = seq(0, 1000, by = 500),
           addGrid = FALSE, xlab = "x", ylab = "y")
plotBinning(randBin.chi[[3]], pch = 19, cex = 0.25, add = TRUE,
            fill = depthFill(randBin.chi[[3]]),
            col = adjustcolor("black", 0.5))
dev.off()

## data with perfect rank agreement
linex <- 1:1e3
liney <- 1:1e3

## allocate storage
lineBin.chi <- vector(mode = "list", length = length(dep))
## split random data using different rules and depths
for (ii in seq_along(dep)) {
    d <- dep[ii]
    lineBin.chi[[ii]] <- binner(linex, liney, stopper = stopFn,
                                splitter = chiSplit)
}

## plot the line data (Fig 4.1(b))
png("lineData.png", width = 2, height = 2, units = "in", res = 480)
narrowPlot(xgrid = seq(0, 1000, by = 500),
           ygrid = seq(0, 1000, by = 500),
           addGrid = FALSE, xlab = "x", ylab = "y")
points(linex, liney, pch = 19, cex = 0.25,
       col = adjustcolor("gray50", 0.5))
dev.off()

## the early progression of splitting (Figs 4.2(b), 4.3(b))
png("lineSplit1.png", width = 2, height = 2, units = "in", res = 480)
narrowPlot(xgrid = seq(0, 1000, by = 500),
           ygrid = seq(0, 1000, by = 500),
           addGrid = FALSE, xlab = "x", ylab = "y")
plotBinning(lineBin.chi[[1]], pch = 19, cex = 0.25, add = TRUE,
            col = adjustcolor("gray50", 0.5))
dev.off()
png("lineSplit2.png", width = 2, height = 2, units = "in", res = 480)
narrowPlot(xgrid = seq(0, 1000, by = 500),
           ygrid = seq(0, 1000, by = 500),
           addGrid = FALSE, xlab = "x", ylab = "y")
plotBinning(lineBin.chi[[2]], pch = 19, cex = 0.25, add = TRUE,
            col = adjustcolor("gray50", 0.5))
dev.off()

## plot bins coloured by depth (Fig 4.4(b))
png("lineSplitDepth.png", width = 2, height = 2, units = "in",
    res = 480)
narrowPlot(xgrid = seq(0, 1000, by = 500),
           ygrid = seq(0, 1000, by = 500),
           addGrid = FALSE, xlab = "x", ylab = "y")
plotBinning(lineBin.chi[[3]], pch = 19, cex = 0.25, add = TRUE,
            fill = depthFill(lineBin.chi[[3]]),
            col = adjustcolor("black", 0.5))
dev.off()

## compute pearson residuals for both and take the maximum magnitude
## to ensure the shading is consistent between them
maxRes <- max(abs(c(binChi(lineBin.chi[[3]])$residuals,
                    binChi(randBin.chi[[3]])$residuals)))
## plot the perason residuals for both (Fig 4.5)
png("randomSplitResid.png", width = 2, height = 2, units = "in",
    res = 480)
narrowPlot(xgrid = seq(0, 1000, by = 500),
           ygrid = seq(0, 1000, by = 500),
           addGrid = FALSE, xlab = "x", ylab = "y")
plotBinning(randBin.chi[[3]], pch = 19, cex = 0.25, add = TRUE,
            fill = residualFill(randBin.chi[[3]], maxRes = maxRes),
            col = adjustcolor("black", 0.5))
dev.off()
png("lineSplitResid.png", width = 2, height = 2, units = "in",
    res = 480)
narrowPlot(xgrid = seq(0, 1000, by = 500),
           ygrid = seq(0, 1000, by = 500),
           addGrid = FALSE, xlab = "x", ylab = "y")
plotBinning(lineBin.chi[[3]], pch = 19, cex = 0.25, add = TRUE,
            fill = residualFill(lineBin.chi[[3]], maxRes = maxRes),
            col = adjustcolor("black", 0.5))
dev.off()



## INVESTIGATING THE NULL DISTRIBUTION ###############################
set.seed(506391)
n <- 1e4 # the sample size
nsim <- 1e4 # number of samples
depths <- 2:10 # range of depths
simDataSets <- replicate(nsim, data.frame(x = sample(1:n),
                                          y = sample(1:n)),
                         simplify = FALSE) # simulated null samples
## allocate storage space
depthSeq.chi <- array(NA, dim = c(4, length(depths), nsim))
depthSeq.mi <- array(NA, dim = c(4, length(depths), nsim))
depthSeq.rnd <- array(NA, dim = c(4, length(depths), nsim))
depthSeq.runif <- array(NA, dim = c(4, length(depths), nsim))
## set criteria/stop function
crits <- makeCriteria(depth >= dep, expn <= 10, n == 0)
stopFn <- function(bns) stopper(bns, crits)
## and splitting functions
chiSplit <- function(bn) maxScoreSplit(bn, chiScores, minExp = 5)
miSplit <- function(bn) maxScoreSplit(bn, miScores, minExp = 5)
rndSplit <- function(bn) maxScoreSplit(bn, randScores, minExp = 5)
runfSplit <- function(bn) rUnifSplit(bn, minExp = 5)

## section which simulates null distribution
if (writeout) { # run simulation and write it out
    ## iterate through the data and recursively bin
    for (ii in 1:nsim) { # each simulation
        for (dep in depths) { # each depth
            depInd <- match(dep, depths) # get storage list index
            ## split using chi scores
            chtr <- binner(simDataSets[[ii]]$x, simDataSets[[ii]]$y,
                           stopper = stopFn, splitter = chiSplit)
            ## compute and store results
            depthSeq.chi[,depInd,ii] <-
                c(chi = binChi(chtr)$stat,
                  mi = binMI(chtr)$stat,
                  nbin = length(chtr),
                  maxDep = max(sapply(chtr,function(bn) bn$depth)))
            ## repeat this for splitting based on the mut. inf.
            mitr <- binner(simDataSets[[ii]]$x, simDataSets[[ii]]$y,
                           stopper = stopFn, splitter = miSplit)
            depthSeq.mi[,depInd,ii] <- # store results
                c(chi = binChi(mitr)$stat,
                  mi = binMI(mitr)$stat,
                  nbin = length(mitr),
                  maxDep = max(sapply(mitr, function(bn) bn$depth)))
            ## random splits
            rntr <- binner(simDataSets[[ii]]$x, simDataSets[[ii]]$y,
                           stopper = stopFn, splitter = rndSplit)
            depthSeq.rnd[,depInd,ii] <-
                c(chi = binChi(rntr)$stat,
                  mi = binMI(rntr)$stat,
                  nbin = length(rntr),
                  maxDep = max(sapply(rntr, function(bn) bn$depth)))
            ## finally, random uniform splits
            rntr <- binner(simDataSets[[ii]]$x, simDataSets[[ii]]$y,
                           stopper = stopFn, splitter = runfSplit)
            depthSeq.runif[,depInd,ii] <-
                c(chi = binChi(rntr)$stat,
                  mi = binMI(rntr)$stat,
                  nbin = length(rntr),
                  maxDep = max(sapply(rntr, function(bn) bn$depth)))
        }
        ## report back on progress
        if (ii %% 50 == 0) cat(paste("\r Done simulated data set:",
                                     ii))
    }
    ## name the data held in each array slice
    dimnames(depthSeq.chi) <- list(c("chi", "mi", "nbin", "maxDep"))
    dimnames(depthSeq.mi) <- list(c("chi", "mi", "nbin", "maxDep"))
    dimnames(depthSeq.rnd) <- list(c("chi", "mi", "nbin", "maxDep"))
    dimnames(depthSeq.runif) <- list(c("chi", "mi", "nbin", "maxDep"))
    ## save the data
    saveRDS(list(depths = depths, chiSplit = depthSeq.chi,
                 miSplit = depthSeq.mi, randSplit = depthSeq.rnd,
                 runifSplit = depthSeq.runif),
            file = paste0("null", n, ".Rds"))
}

## plot the paths: statistic by number of bins and depth for the three
## different splitting rules (Figs 4.6, 4.7, 4.8)
## read in the simulation rather than run it every time
#load("null10000.rda")
data <- readRDS("null10000.Rds")
depths <- data$depths
for (spltr in c("chiSplit", "miSplit", "randSplit", "runifSplit")) {
    png(paste0(spltr, "ChiDepth.png"), width = 4, height = 4,
               units = "in", res = 480, bg = "transparent")
    narrowPlot(xgrid = seq(0, 3, by = 0.5), # plot region
               xlab = expression(log[10]~"K"),
               ygrid = seq(-1, 4, by = 1), bg = "transparent",
               ylab = expression(log[10]~{"("~X^2~")"}))
    points(log(data[[spltr]]["nbin",,],10), # points
           log(data[[spltr]]["chi",,],10),
           col = adjustcolor(depthPal, 0.1),
           pch = 20)
    ## mean points
    for (ii in 1:9) {
        points(log(mean(data[[spltr]]["nbin",ii,]),10),
               log(mean(data[[spltr]]["chi",ii,]),10),
               pch = 3, cex = 0.5, lwd = 1.2, col = "black")
        points(log(mean(data[[spltr]]["nbin",ii,]),10),
               log(mean(data[[spltr]]["chi",ii,]),10),
               pch = 3, cex = 0.5, lwd = 0.8, col = depthPal[ii])
    }
    for (p in c(0.01)) { # chi quantiles
        lines(log(2:600,10), log(qchisq(1-p, 1:599),10),
              lty = 2)
    }
    ## position the legend based on the splitting rule
    if (grepl("rand", "runif", spltr)) {
        legpos <- "topleft"
    } else legpos <- "bottomright"
    legend(x = legpos, legend = 2:10, title = "Max depth",
           bg = "transparent", cex = 0.8,
           fill = adjustcolor(depthPal, 0.6))
    dev.off()
}

## smaller plots to compare different sample sizes (Figs 4.10, 4.11,
## 4.12)
size <- 1.8
for (n in c(1e2, 1e3, 1e4)) {
    ## load corresponding data
    #data <- load(paste0("null", n, ".rda"))
    data <- readRDS(paste0("null", n, ".Rds"))
    depths <- data$depths
    if (n == 1e2) {
        mar <- c(2.1, 2.1, 0.1, 0.1)
        wid <- size + 0.1
    } else if (n == 1e3) {
        mar <- c(2.1, 1.1, 0.1, 0.1)
        wid <- size - 0.1
    } else {
        mar <- c(2.1, 1.1, 0.1, 3.1)
        wid <- size + 0.4
    }
    for (spltr in c("chiSplit", "miSplit", "randSplit",
                    "runifSplit")) {
        png(paste0(spltr, "ChiDepth", n, ".png"), width = wid,
            height = size, units = "in", res = 480, bg = "transparent")
        narrowPlot(xgrid = seq(0, 3, by = 1), bg = "transparent",
                   ygrid = seq(-1, 4, by = 1),
                   xlab = expression(log[10]~"K"),
                   ylab = expression(log[10]~{"("~X^2~")"}),
                   mars = mar)
        points(log(data[[spltr]]["nbin",,],10),
               log(data[[spltr]]["chi",,],10),
               col = adjustcolor(depthPal, 0.1),
               cex = 0.5, pch = 20)
        for (ii in 1:9) {
            points(log(mean(data[[spltr]]["nbin",ii,]),10),
                   log(mean(data[[spltr]]["chi",ii,]),10),
                   pch = 3, cex = 0.5, lwd = 1.2, col = "black")
            points(log(mean(data[[spltr]]["nbin",ii,]),10),
                   log(mean(data[[spltr]]["chi",ii,]),10),
                   pch = 3, cex = 0.5, lwd = 0.8, col = depthPal[ii])
        }
        for (p in c(0.01)) {
            lines(log(2:600,10), log(qchisq(1-p, 1:599),10), lty = 2)
        }
        if (n == 1e4) {
            bnds <- par()$usr
            legend(x = bnds[2], y = bnds[4], legend = 2:10,
                   title = "Max depth", # bg = "white",
                   cex = 0.6, xpd = NA, bg = "transparent",
                   fill = adjustcolor(depthPal, 0.6))
        }
        dev.off()
    }
}

## quantile regression, where the quantiles of the chi statistic are
## predicted based on the number of bins, makes for an easier
## comparison of the upper ends of the distributions
library(quantreg)
qnts <- c(0.95, 0.99, 0.999)
## choose the data with a sample size of 1000
load("null10000.rda")
data <- get("null10000")
depths <- data$depths
depthSeq.chi <- data$chiSplit
depthSeq.mi <- data$miSplit
depthSeq.rnd <- data$randSplit # split sections of the data
## random splitting
rndChiQnt <- rq(chi ~ nbin, tau = qnts,
                data = data.frame(nbin = c(depthSeq.rnd["nbin",,]),
                                  chi = c(depthSeq.rnd["chi",,])))
## chi splitting
chiChiQnt <- rq(chi ~ nbin, tau = qnts,
                data = data.frame(nbin = c(depthSeq.chi["nbin",,]),
                                  chi = c(depthSeq.chi["chi",,])))
## quantile of chi under mi splitting
chiMiQnt <- rq(chi ~ nbin, tau = qnts,
               data = data.frame(nbin = c(depthSeq.mi["nbin",,]),
                                 chi = c(depthSeq.mi["chi",,])))

## plot the quantiles (Fig 4.9)
xseq <- seq(4, 600, by = 1)
randQntPred <- predict(rndChiQnt, newdata = data.frame(nbin = xseq))
chiQntPred <- predict(chiChiQnt, newdata = data.frame(nbin = xseq))
miQntPred <- predict(chiMiQnt, newdata = data.frame(nbin = xseq))
png("binQuantileRegression.png", width = 3, height = 3, res = 480,
    units = "in")
narrowPlot(xgrid = seq(1, 3, by = 1), ylim = c(0, 3.5),
           xlab = expression(log[10]~"(Number of bins)"),
           ygrid = seq(0, 4, by = 1),
           ylab = expression(log[10]~{"("~chi^2~statistic~")"}),
           xlim = c(0.5, 3), mars = c(2.1, 2.1, 0.1, 0.1))
for (ii in 1:3) {
    lines(log(xseq, 10), log(randQntPred[,ii], 10),
          lty = ii, col = adjustcolor("steelblue", 0.8))
    lines(log(xseq, 10), log(chiQntPred[,ii], 10),
          lty = ii, col = adjustcolor("firebrick", 0.8))
    lines(log(xseq, 10), log(miQntPred[,ii], 10),
          lty = ii, col = adjustcolor("seagreen", 0.8))
    lines(log(xseq, 10), log(qchisq(qnts[ii], df = xseq), 10),
          lty = ii, col = adjustcolor("black", 0.8))
    text(log(xseq[1], 10), log(randQntPred[1,ii], 10),
         labels = paste0(100*qnts[ii], "%"), cex = 0.6,
         adj = c(0.5,-0.1))
}
legend(x = "bottomright",
       legend = c("Random bins", "Max MI bins",
                  expression("Max "~chi~bins),
                  expression(chi^2~" critical value")),
       cex = 0.8, col = adjustcolor(c("steelblue", "seagreen",
                                      "firebrick", "black"),
                                    0.8),
       lty = c(1,1,1,1), bg = "white")
dev.off()



## SIMULATED DATA PATTERNS ###########################################
## patterns from Newton (2009) provided in a list of functions
patFns <- list(
    wave = function(n) {
        x <- seq(-1, 1, length=n)
        u <- x + runif(n)/3; v <- 4*((x^2 - 1/2)^2 + runif(n)/500)
        cbind(x = u, y = v)
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
    saddle = function(n) {
        x <- seq(-1,1, length=n )
        y <- (x ^2 + runif(n))/2
        cbind(x = x, y = y)
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
    rotSquare = function(n) {
        x <- runif(n, min = -1, max = 1)
        y <- runif(n, min = -1, max = 1)
        theta <--pi/8
        rr <- rbind(c(cos(theta), -sin(theta)),
                    c(sin(theta), cos(theta)))
        tmp <- cbind(x, y) %*% rr
        colnames(tmp) <- c("x",  "y")
        tmp
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

## plot the first data realization (Fig 4.12)
m <- 1
pal <- c(RColorBrewer::brewer.pal(6, "Pastel2"), "gray50")
for(i in 1:7) {
    png(file=paste0("patterns-", names(patFns)[i], ".png"), height = m,
        width = m, units = "in", res = 480, bg = "transparent")
    par(mar=c(1,1,1,1)/2)
    plot(simData[, "x", i, 1], simData[, "y", i, 1], xlab = "",
         ylab = "", axes = F, pch = 19, cex = 0.2, col = pal[i],
         bg = "transparent")
    dev.off()
}

## convert this data into pairwise ranks and plot it (Fig 4.13)
simXr <- apply(simData[, "x", , ], c(2, 3), rank)
simYr <- apply(simData[, "y", , ], c(2, 3), rank)
for(i in 1:7) {
    png(file=paste0("ranks-", names(patFns)[i], ".png"), height = m,
        width = m, units = "in", res = 480, bg = "transparent")
    par(mar=c(1,1,1,1)/2)
     plot(simXr[, i, 1], simYr[, i, 1], xlab = "", ylab = "",
          axes= F, pch = 19, cex = 0.2, col = pal[i],
          bg = "transparent")
    dev.off()
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

## bin each realization (takes some time)
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
load("null1000.rda") # read in null
data <- get("null1000")
depths <- data$depths
depthSeq.chi <- data$chiSplit
depthSeq.mi <- data$miSplit
depthSeq.rnd <- data$randSplit # null data

## plot paths for an individual random split (Figure 4.14(b))
png("simDataRandPath.png", width = 3, height = 3, units = "in",
    res = 480, bg = "transparent")
narrowPlot(xgrid = seq(0, 160, by = 40),
           xlab = expression(K),
           ygrid = seq(0, 1600, by = 400),
           ylab = expression(X^2), bg = "transparent")
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
dev.off()

## make the same plot for paths from chi splitting (Figure 4.14(a))
png("simDataMaxChiPath.png", width = 3, height = 3, units = "in",
    res = 480, bg = "transparent")
narrowPlot(xgrid = seq(0, 160, by = 40),
           xlab = expression(K),
           ygrid = seq(0, 1600, by = 400),
           ylab = expression(X^2), bg = "transparent")
for (ii in 1:1e4) { # null lines
    lines(depthSeq.chi["nbin",,ii],
          depthSeq.chi["chi",,ii],
          col = adjustcolor("gray", 0.1))
}
for (jj in 1:7) { # observed path
    lines(chiNbin[[1]][jj,], chiPaths[[1]][jj,], col = pal[jj])
    points(chiNbin[[1]][jj,], chiPaths[[1]][jj,], col = pal[jj],
           pch = 19, cex = 0.5)
}
dev.off()

## for the random split repetitions, plot every one (Fig 4.15(b))
png("simDataRandAll.png", width = 5, height = 3, units = "in",
    res = 480, bg = "transparent")
narrowPlot(xgrid = seq(0, 160, by = 40),
           xlab = expression(K),
           ygrid = seq(0, 1600, by = 400),
           ylab = expression(X^2), bg = "transparent")
for (ii in 1:1e4) { # null lines
    lines(depthSeq.rnd["nbin",,ii],
          depthSeq.rnd["chi",,ii],
          col = adjustcolor("gray", 0.1))
}
for (jj in 1:7) { # observed paths
    for (ii in 1:nsim) {
        lines(rndNbin[[ii]][jj,], rndPaths[[ii]][jj,],
              col = adjustcolor(pal[jj], 0.2))
    }
}
## get median lines
medianNbin <- apply(simplify2array(rndNbin), c(1,2), median)
medianPaths <-  apply(simplify2array(rndPaths), c(1,2), median)
for (jj in 1:7) { # add mean lines
    lines(medianNbin[jj,], medianPaths[jj,], col = "gray30",
          lwd = 3)
    lines(medianNbin[jj,], medianPaths[jj,], col = pal[jj],
          lwd = 2)
}
lines(1:160, qchisq(0.95, 1:160), lty = 2)
dev.off()

## do the same for the chi splits (Fig 4.15(a))
png("simDataMaxChiAll.png", width = 5, height = 3, units = "in",
    res = 480, bg = "transparent")
narrowPlot(xgrid = seq(0, 160, by = 40), bg = "transparent",
           xlab = expression(K),
           ygrid = seq(0, 1600, by = 400),
           ylab = expression(X^2))
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
## get median lines
medianNbin <- apply(simplify2array(chiNbin), c(1,2), median)
medianPaths <-  apply(simplify2array(chiPaths), c(1,2), median)
for (jj in 1:7) { # add mean lines
    lines(medianNbin[jj,], medianPaths[jj,], col = "gray30",
          lwd = 3)
    lines(medianNbin[jj,], medianPaths[jj,], col = pal[jj],
          lwd = 2)
}
lines(1:160, qchisq(0.95, 1:160), lty = 2)
dev.off()

## plotting bins for every depth: first remind of patterns
for(i in 1:7) {
    png(file=paste0("red-", names(patFns)[i], ".png"), height = m,
        width = m, units = "in", res = 480, bg = "transparent")
    par(mar=c(1,1,1,1)/2)
    plot(simXr[, i, 1], simYr[, i, 1], xlab = "", ylab = "",
          axes= F, pch = 19, cex = 0.2, col = "firebrick",
         bg = "transparent")
    dev.off()
}

## next, check the bins for every depth (Fig 4.16)
## start by getting the maximum residual to make the shading constant
maxRes <- max(sapply(unlist(testChiChi[[1]],
                            recursive = FALSE),
                     getMaxRes))
## for every depth, display the binning for each pattern
for (depth in 2:10) {
    png(file = paste0("simDataBins", depth, ".png"), height = m,
        width= 6*m, units = "in", res = 480, bg = "transparent")
    par(mfrow=c(1,7), mar=c(1,1,1,1)/2)
    for(i in 1:7) {
        plot(NA, ylim = c(1, n), xlim = c(1, n), # remove axes
             axes = F, xlab = "", ylab = "", main = "",
             bg = "transparent")
        plotBinning(testChiBins[[1]][[depth]][[i]], pch = 19,
                    cex = 0.1, add = TRUE,
                    col = NA,
                    fill = residualFill(testChiBins[[1]][[depth]][[i]],
                                        colrng = c("steelblue", "white",
                                                   "firebrick"),
                                        maxRes = maxRes))
    }
    dev.off()
}

## do the same thing for the bins from MI splitting
maxRes <- max(sapply(unlist(testMiChi[[1]],
                            recursive = FALSE),
                     getMaxRes))
for (depth in 2:10) {
    png(file = paste0("simDataBins", depth, "MI.png"), height = m,
        width= 6*m, units = "in", res = 480, bg = "transparent")
    par(mfrow=c(1,7), mar=c(1,1,1,1)/2)
    for(i in 1:7) {
        plot(NA, ylim = c(1, n), xlim = c(1, n),
             axes = F, xlab = "", ylab = "", main = "",
             bg = "transparent")
        plotBinning(testMiBins[[1]][[depth]][[i]], pch = 19,
                    cex = 0.1, add = TRUE,
                    col = adjustcolor("grey", 0.8),
                    fill = residualFill(testMiBins[[1]][[depth]][[i]],
                                        maxRes = maxRes))
    }
    dev.off()
}

## and the random bins for every depth (Figure 4.17 for depth of 10)
## again, standardize the residuals
maxRes <- max(sapply(unlist(testRndChi[[10]],
                            recursive = FALSE),
                     getMaxRes))
for (depth in 2:10) {
    png(file = paste0("simDataBinsRand", depth, ".png"), height = m,
        width= 6*m, units = "in", res = 480, bg = "transparent")
    par(mfrow=c(1,7), mar=c(1,1,1,1)/2)
    for(i in 1:7) {
        plot(NA, ylim = c(1, n), xlim = c(1, n),
             axes = F, xlab = "", ylab = "", main = "",
             bg = "transparent")
        plotBinning(testRndBins[[10]][[depth]][[i]],
                    pch = 19, cex = 0.1, add = TRUE,
                    col = NA,
                    fill = residualFill(testRndBins[[10]][[depth]][[i]],
                                        maxRes = maxRes))
    }
    dev.off()
}


## REAL DATA EXAMPLE #################################################
## "SP500" demo in "zenplots" package, code from Marius
## Hofert, produces a set of pseudo-observations that are uniform
## these are loaded here and converted to ranks
data(sp500pseudo)
spRanks <- apply(sp500pseudo, 2, rank, ties.method = "random")
rownames(spRanks) <- NULL
spPairs <- combn(ncol(spRanks), 2) # all possible pairs
## next, we iterate through all pairs and bin to a maximum depth of 6
## define the criteria to used
crits <- makeCriteria(depth >= 6, expn <= 10, n == 0)
stopFn <- function(bns) stopper(bns, crits)
## and potential splitting functions
chiSplit <- function(bn) maxScoreSplit(bn, chiScores, minExp = 5)
miSplit <- function(bn) maxScoreSplit(bn, miScores, minExp = 5)
rndSplit <- function(bn) maxScoreSplit(bn, randScores, minExp = 5)
## allocate storage
spBins <- vector("list", ncol(spPairs))
msgInd <- ((1:ncol(spPairs)) %% 1000) == 0
## iterate through all pairs
## ~ 57 mins
set.seed(85912024)
system.time({for (ii in seq_len(ncol(spPairs))) { ## ~57 mins
    pair <- spPairs[, ii] # indices of pairs
    spBins[[ii]] <- binner(spRanks[, pair[1]], spRanks[, pair[2]],
                           stopper = stopFn,
                           splitter = rndSplit)
    if (msgInd[ii]) {
        cat("\r Completed ", ii, " pairs")
    }
             }})
## drop points for smaller storage size
spBinsNP <- lapply(spBins, dropBinPoints)
## save binnings
saveRDS(spBinsNP, file = paste0("sp500binsRnd", "NoPts.Rds"))

## plot observed statistics on the S&P500 data for the random or
## maximized splits
random <- FALSE
if (random) { ## load pre-processed data
    spBinsNP <- readRDS("sp500binsRndNoPts.Rds")
} else {
    spBinsNP <- readRDS("sp500binsNoPts.Rds")
}
## get chi statistics across the bins
spChis <- lapply(spBinsNP, function(bns) binChi(bns))
spChiStats <- sapply(spChis, function(x) x$stat)
spChiResid <- sapply(spChis, function(x) x$residuals)
spChiNbin <- sapply(spChiResid, length)
## order by most interesting
spOrd <- order(spChiStats, decreasing = TRUE)
spMaxRes <- 17.35 # hard coded: largest residual across both

## plot the statistic quantiles
qnts <- c(0.95, 0.99, 0.999) # chosen quantiles
if (random) {
    nullP <- pchisq(spChiStats, spChiNbin - 1)
    spRGB <- colorRamp(c("steelblue", "firebrick"),
                       bias = 10)(nullP^2)/255
} else {
    library(quantreg)
    ## the null density
    nulls <- readRDS("SplitsRandomDatan1000.Rds")
    ## splitting the null statistics by number of bins allows us to more
    ## easily compute the empirical quantiles of the sp500 data
    splitNulls <- split(c(nulls$chiSplit["chi",c(5,6),]),
                        c(nulls$chiSplit["nbin",c(5,6),]))
    binVals <- as.numeric(names(splitNulls))
    ## compute the empirical p-values by comparing each of the observed
    ## statistics to the corresponding split bin of the null distribution
    empP <- sapply(seq_along(spChiStats),
                   function(ii) {
                       nb <- as.character(spChiNbin[ii])
                       sum(spChiStats[ii] >
                           splitNulls[[nb]])/length(splitNulls[[nb]])
                   })
    ## convert these to a hue for plotting of points
    spRGB <- colorRamp(c("steelblue", "firebrick"),
                       bias = 10)(empP^2)/255
    ## perform quantile regression for the null data as well
    modQnt <- rq(chi ~ nbin, tau = c(0.95, 0.99, 0.999),
                 data = data.frame(nbin = c(nulls$chiSplit["nbin",,]),
                                   chi = c(nulls$chiSplit["chi",,])))
    predQnt <- predict(modQnt,
                       newdata = data.frame(nbin = binVals))
}
## convert p-values to a colour
spCol <- rgb(spRGB[,1], spRGB[,2], spRGB[,3])

## plot the sp500 point cloud coloured by empirical p-value with the
## quantile regression lines alongside
if (random) {
    nme <- "sp500empPColourRnd.png"
} else {
    nme <- "sp500empPColour.png"
}
png(nme, width = 3.5, height = 3.5, units = "in",
    res = 480)
narrowPlot(xgrid = seq(1, 2, by = 0.5),
           ygrid = seq(1, 3, by = 1), ylim = c(1, 3.2),
           xlab = expression(log[10]~"("~n[bin]~")"),
           ylab = expression(log[10]~{"("~X^2~")"}),
           mars = c(2.1, 2.1, 2.1, 2.1))
points(log(spChiNbin, 10), log(spChiStats, 10), cex = 0.5,
       col = adjustcolor(spCol, 0.2), pch = 20)
for (ii in 1:3) {
    yadj <- 1 - 0.5*(ii - 1)
    if (random) {
        lines(log(seq(10^1.15, 10^2, 1), 10),
              log(qchisq(qnts[ii], df = seq(10^1.15, 10^2, 1) - 1), 10),
              lty = ii)
        text(labels = paste0(100*qnts[ii], "%"),
             x = 1.15, y = log(qchisq(qnts[ii], df = 10^1.15-1), 10),
             adj = c(1, yadj), cex = 0.6)
    } else{
        lines(log(binVals[-c(1,2)], 10),
              log(predQnt[-c(1,2),ii], 10), lty = ii)
        text(labels = paste0(100*qnts[ii], "%"),
             x = log(binVals[3], 10), y = log(predQnt[3, ii], 10),
             adj = c(1, yadj), cex = 0.6)
    }
}
addMarHists(log(spChiNbin, 10), log(spChiStats, 10),
            xcuts = seq(1, 2, by = 0.03125),
            ycuts = seq(1, 3, by = 0.0625))
dev.off()

## compare the values between the two
spChiStatsBoth <- readRDS("sp500Stats.Rds")
spChiPBoth <- readRDS("sp500Pvals.Rds")
spChiStatRanks <- apply(spChiStatsBoth, 2, rank,
                        ties.method = "random")
rejCol <- c("black", "firebrick")[((spChiPBoth[, "max"] > 0.95) !=
                                   (spChiPBoth[, "rand"] > 0.95)) + 1]

## plot it all
png("spRandvsMaxRanks.png", width = 3.5, height = 3.5, units = "in",
    res = 480)
narrowPlot(xgrid = seq(0, 100000, by = 50000),
           ygrid = seq(0, 100000, by = 50000),
           ylim = c(0, 106030), xlim = c(0, 106030),
           xlab = expression({"Random"~X^2~"rank"}),
           ylab = expression({"Max"~X^2~"rank"}),
           mars = c(2.1, 2.1, 2.1, 2.1))
points(spChiStatRanks[, "rand"], spChiStatRanks[, "max"],
       pch = ".", col = adjustcolor("black", 0.3))
points(spChiStatRanks[c(22451, 65548), "rand"],
       spChiStatRanks[c(22451, 65548), "max"],
       pch = 1, col = "firebrick")
text(labels = c("SNA:PH", "MCO:AEE"), x = spChiStatRanks[c(22451, 65548), "rand"],
     y = spChiStatRanks[c(22451, 65548), "max"], adj = c(-0.1,0.5),
     cex = 0.6)
                                        #col = rejCol)
dev.off()

## which pairs maximize the difference in either direction?
rankDiffs <- spChiStatRanks[, "max"] - spChiStatRanks[, "rand"]
rankDiffOrd <- order(rankDiffs, decreasing = TRUE)
topMax <- rankDiffOrd[1]
topRand <- rankDiffOrd[length(rankDiffOrd)]
## plot these two
png("spBiggestDiffs.png", width = 3.5, height = 3.5, units = "in",
    res = 480)
par(mfrow = c(2, 2), mar = c(0.1, 0.55, 1.1, 0.55))
for (prInd in c(topMax, topRand)) {
    pr <- spPairs[, prInd] # pair indices
    tempBin <- binner(spRanks[, pr[1]], spRanks[, pr[2]],
                      stopper = stopFn, splitter = chiSplit)
    bns <- list("Max" = tempBin, "Random" = spBinsNP[[prInd]])
    for (ii in seq_along(bns)) {
        plot(NA, xlim = c(1, (nrow(spRanks))),
             ylim = c(1, (nrow(spRanks))), axes = "F",
             xlab = colnames(spRanks)[pr[1]],
             ylab = colnames(spRanks)[pr[2]],
             main = "")
        mtext(paste0(paste(colnames(spRanks)[pr], collapse = ":"),
                     " (", names(bns)[ii], ")"),
              cex = 0.6)
        plotBinning(bns[[ii]], pch = NA,
                    fill = residualFill(bns[[ii]],
                                        maxRes = spMaxRes,
                                        colrng = c("steelblue",
                                                   "white",
                                                   "firebrick")),
                    add = TRUE)
        points(spRanks[, pr[1]], spRanks[, pr[2]], pch = ".",
               col = adjustcolor("gray50"))
    }
}
dev.off()

## these seem to have some association... is the low rank in the
## random case a fluke? look at the most egregious difference
testStat <- numeric(200)
for (ii in seq_along(testStat)) {
    tempBin <- binner(spRanks[, spPairs[1, 22451]],
                      spRanks[, spPairs[2, 22451]],
                      stopper = stopFn,
                      splitter = rndSplit)
    testStat[ii] <- binChi(tempBin)$stat
}


## use this to plot the top pairs and their binnings
if (random) {
    nme <- "sp500top16Rnd.png"
} else {
    nme <- "sp500top16.png"
}
png(nme, width = 3.5, height = 3.5, units = "in",
    res = 480)
par(mfrow = c(4, 4), mar = c(0.1, 0.55, 1.1, 0.55))
for (prInd in spOrd[1:16]) {
    pr <- spPairs[, prInd] # pair indices
    plot(NA, xlim = c(1, (nrow(spRanks))),
         ylim = c(1, (nrow(spRanks))), axes = "F",
         xlab = colnames(spRanks)[pr[1]],
         ylab = colnames(spRanks)[pr[2]],
         main = "")
    mtext(paste(colnames(spRanks)[pr], collapse = ":"),
          cex = 0.6)
    tempBin <- binner(spRanks[, pr[1]], spRanks[, pr[2]],
                      stopper = stopFn,
                      splitter = chiSplit)
    plotBinning(tempBin, pch = NA,
                fill = residualFill(tempBin,
                                    maxRes = spMaxRes,
                                    colrng = c("steelblue",
                                               "white",
                                               "firebrick")),
                add = TRUE)
    points(spRanks[, pr[1]], spRanks[, pr[2]], pch = ".",
           col = adjustcolor("gray50"))
}
dev.off()

## do the same for pairs in the middle of the distribution
png("sp500mid36.png", width = 5, height = 5, units = "in",
    res = 480)
par(mfrow = c(6, 6), mar = c(0.1, 0.55, 1.1, 0.55))
for (prInd in spOrd[seq(length(spOrd)/2 - 17, by = 1,
                        length.out = 36)]) {
    pr <- spPairs[, prInd] # pair indices
    plot(NA, xlim = c(1, (nrow(spRanks))),
         ylim = c(1, (nrow(spRanks))), axes = "F",
         xlab = colnames(spRanks)[pr[1]],
         ylab = colnames(spRanks)[pr[2]],
         main = "")
    mtext(paste(colnames(spRanks)[pr], collapse = ":"),
          cex = 0.6)
    plotBinning(spBinsNP[[prInd]],
                fill = residualFill(spBinsNP[[prInd]],
                                    maxRes = spMaxRes),
                add = TRUE)
    points(spRanks[, pr[1]], spRanks[, pr[2]], pch = ".",
           col = adjustcolor("gray50"))
}
dev.off()

## finally view the weakest associations
png("sp500last36.png", width = 5, height = 5, units = "in",
    res = 480)
par(mfrow = c(6, 6), mar = c(0.1, 0.55, 1.1, 0.55))
for (prInd in spOrd[seq(length(spOrd)-35, by = 1,
                        length.out = 36)]) {
    pr <- spPairs[, prInd] # pair indices
    plot(NA, xlim = c(1, (nrow(spRanks))),
         ylim = c(1, (nrow(spRanks))), axes = "F",
         xlab = colnames(spRanks)[pr[1]],
         ylab = colnames(spRanks)[pr[2]],
         main = "")
    mtext(paste(colnames(spRanks)[pr], collapse = ":"),
          cex = 0.6)
    plotBinning(spBinsNP[[prInd]],
                fill = residualFill(spBinsNP[[prInd]],
                                    maxRes = spMaxRes),
                add = TRUE)
    points(spRanks[, pr[1]], spRanks[, pr[2]], pch = ".",
           col = adjustcolor("gray50"))
}
dev.off()

## a twist: try this on unprocessed S&P500 data
spRaw <- readRDS("sp500raw.Rds")
spRawRanks <- apply(spRaw, 2, rank, ties.method = "random")
rownames(spRawRanks) <- NULL
spPairs <- combn(ncol(spRawRanks), 2)
## set up recursive binning
crits <- makeCriteria(depth >= 6, expn <= 10, n == 0)
stopFn <- function(bns) stopper(bns, crits)
## and potential splitting functions
chiSplit <- function(bn) maxScoreSplit(bn, chiScores, minExp = 5)
rndSplit <- function(bn) maxScoreSplit(bn, randScores, minExp = 5)
## set up storage, messages
spBins <- vector("list", ncol(spPairs))
msgInd <- ((1:ncol(spPairs)) %% 1000) == 0
## iterate through all pairs
## ~ 57 mins
set.seed(8190525)
system.time({for (ii in seq_len(ncol(spPairs))) { ## ~57 mins
    pair <- spPairs[, ii] # indices of pairs
    spBins[[ii]] <- binner(spRawRanks[, pair[1]],
                           spRawRanks[, pair[2]],
                           stopper = stopFn,
                           splitter = chiSplit)
    if (msgInd[ii]) {
        cat("\r Completed ", ii, " pairs")
    }
             }})
## drop points for smaller storage size
spBinsNP <- lapply(spBins, dropBinPoints)
## save binnings
saveRDS(spBinsNP, file = paste0("sp500binsChi-raw", "NoPts.Rds"))

## plot observed statistics on the S&P500 data for the random or
## maximized splits
random <- FALSE
if (random) { ## load pre-processed data
    spBinsNPRaw <- readRDS("sp500binsRnd-rawNoPts.Rds")
} else {
    spBinsNPRaw <- readRDS("sp500binsChi-rawNoPts.Rds")
}
## get chi statistics across the bins
spChisRaw <- lapply(spBinsNPRaw, function(bns) binChi(bns))
spChiStatsRaw <- sapply(spChisRaw, function(x) x$stat)
spChiResidRaw <- sapply(spChisRaw, function(x) x$residuals)
spChiNbinRaw <- sapply(spChiResidRaw, length)
## order by interestingness
spOrdRaw <- order(spChiStatsRaw, decreasing = TRUE)
## try plotting this
plot(spOrdRaw, spOrd, pch = ".")
