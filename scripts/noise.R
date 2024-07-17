library(AssocBin)
library(RColorBrewer)

## increasing noise simulation
n <- 1000
nsim <- 250
sds <- seq(0, 20, by = 0.2)

## bin settings
stopCrits <- makeCriteria(depth > 6, n < 1, expn <= 10)
stopFn <- function(bns) stopper(bns, stopCrits)
sqrfySplit <- function(bn) rIntSplit(bn, squarify = TRUE)

## a function to run this simulation for a given y function on x
simPattern <- function(dataGen, n = 1000, nsim = 500, 
                       sds = seq(0, 8, by = 0.2), sFn = stopFn,
                       splitter = rIntSplit) {
    simSettings <- expand.grid(sim = 1:nsim, sd = sds)
    pval <- stat <- snr <- rsq <- numeric(nsim)
    for (ii in seq_len(nrow(simSettings))) {
        dat <- dataGen(n)
        x <- dat$x
        y <- dat$y
        yv <- y + rnorm(n, sd = simSettings$sd[ii])
        xv <- x + rnorm(n, sd = simSettings$sd[ii])
        xr <- rank(xv, ties.method = "random")
        yr <- rank(yv, ties.method = "random")
        snr[ii] <- (var(y) + mean(y)^2)/(simSettings$sd[ii])^2
        rsq[ii] <- var(y)/var(yv)
        binning <- binner(xr, yr, splitter = splitter, stopper = sFn)
        stat[ii] <- binChi(binning)$stat
        pval[ii] <- pchisq(stat[ii], df = length(binning) - 1,
                           lower.tail = FALSE, log.p = TRUE)
    }
    return(list(rsq = rsq, stat = stat, pval = pval, sd = simSettings$sd,
                snr = snr))
}

## a slightly different one for functional patterns
simFunc <- function(dataGen, n = 1000, nsim = 500, 
                    sds = seq(0, 8, by = 0.2), sFn = stopFn,
                    splitter = rIntSplit) {
    simSettings <- expand.grid(sim = 1:nsim, sd = sds)
    pval <- stat <- snr <- rsq <- numeric(nsim)
    for (ii in seq_len(nrow(simSettings))) {
        dat <- dataGen(n)
        x <- dat$x
        y <- dat$y
        yv <- y + rnorm(n, sd = simSettings$sd[ii])
        xr <- rank(x, ties.method = "random")
        yr <- rank(yv, ties.method = "random")
        rsq[ii] <- var(y)/var(yv)
        binning <- binner(xr, yr, splitter = splitter, stopper = sFn)
        stat[ii] <- binChi(binning)$stat
        pval[ii] <- pchisq(stat[ii], df = length(binning) - 1,
                           lower.tail = FALSE, log.p = TRUE)
    }
    return(list(rsq = rsq, stat = stat, pval = pval, sd = simSettings$sd))
}

## median plot
medianPlot <- function(x, y, breaks = 100, add = FALSE, ...) {
    xbr <- cut(x, breaks = breaks)
    ymed <- tapply(y, xbr, median)
    xmed <- sapply(strsplit(gsub("\\(|\\]", "", names(ymed)), ","), 
                   function(x) mean(as.numeric(x)))
    if (add) {
        points(x = xmed, y = ymed, ...)
    } else {
        plot(x = xmed, y = ymed, ...)
    }
}

## plot example patterns
plotEx <- function(dataGen, dat, n = 1000, rsq = 0.2,
                   tol = 0.005, ...) {
    ## get noise sds
    sd <- mean(dat$sd[which(abs(dat$rsq - rsq) <= tol)])
    ex <- dataGen(n)
    cury <- ex$y + rnorm(n, sd = sd)
    curRsq <- var(ex$y)/var(cury)
    while (abs(rsq - curRsq) > 0.001) {
        ex <- dataGen(n)
        cury <- ex$y + rnorm(n, sd = sd)
        curRsq <- var(ex$y)/var(cury)
    }
    plot(x = ex$x + rnorm(n, sd = sd), y = cury, ...)
}

## different patterns ##########################################
## functional
## sine wave
genSine <- function(n, rng = c(0, 6*pi)) {
    x <- runif(n, min = rng[1], max = rng[2])
    y <- sin(x)*diff(rng)/2 + rng[1]
    list(x = x, y = y)
}

## parabola
genPar <- function(n, rng = c(0, 6*pi)) {
    x <- runif(n, min = rng[1], max = rng[2])
    y <- diff(rng)*(x - mean(rng))^2/(mean(rng)^2)
    list(x = x, y = y)
}

## cubic
genCubic <- function(n, rng = c(0, 6*pi)) {
    x <- runif(n, min = rng[1], max = rng[2])
    xadj <- x - mean(rng)
    scl <- diff(rng)/2
    y <- diff(rng)*(20*xadj^3/scl^3 - 6*xadj^2/scl^2
        + xadj/scl)/42
    list(x = x, y = y)
}

## exponential
genExponential <- function(n, rng = c(0, 6*pi)) {
    x <- runif(n, min = rng[1], max = rng[2])
    y <- diff(rng)*exp(x)/exp(diff(rng))
    list(x = x, y = y)
}

## straight line
genLine <- function(n, rng = c(0, 6*pi)) {
    x <- runif(n, min = rng[1], max = rng[2])
    y <- x
    list(x = x, y = y)
}

## non-coexistence
genNonco <- function(n, rng = c(0, 6*pi)) {
    x <- c(rep(0, n/2), runif(n/2, min = rng[1], max = rng[2]))
    y <- runif(n, min = rng[1], max = rng[2])*(x == 0)
    list(x = x, y = y)
}

## non-functional
## two lines
genTwoLines <- function(n, rng = c(0, 6*pi)) {
    x <- runif(n, min = rng[1], max = rng[2])
    u <- as.numeric(runif(n) >= 0.5) + 1
    y <- x*c(0.2, 1)[u]
    list(x = x, y = y)
}

## X
genX <- function(n, rng = c(0, 6*pi)) {
    x <- runif(n, min = rng[1], max = rng[2])
    u <- as.numeric(runif(n) >= 0.5) + 1
    y <- x*c(1, -1)[u] + c(0, max(x))[u]
    list(x = x, y = y)
}

## ellipse
genEllipse <- function(n, rng = c(0, 6*pi)) {
    x <- runif(n, min = rng[1], max = rng[2])
    u <- as.numeric(runif(n) >= 0.5) + 1
    scl <- (rng[2] - rng[1])/2
    y <- c(1, -1)[u]*scl*sqrt((1 - ((x - mean(rng))/(scl))^2))
    list(x = x, y = y)
}

## die four pattern
genFour <- function(n, rng = c(0, 6*pi)) {
    pos1 <- rng[1] + 1/4*diff(rng)
    pos2 <- rng[1] + 3/4*diff(rng)
    x <- c(rep(pos1, n/2), rep(pos2, n/2))
    y <- c(rep(pos1, n/4), rep(pos2, n/4), rep(pos1, n/4),
           rep(pos2, n/4))
    list(x = x, y = y)
}

## die five pattern
genFive <- function(n, rng = c(0, 6*pi)) {
    pos1 <- rng[1] + 1/4*diff(rng)
    pos2 <- rng[1] + 3/4*diff(rng)
    x <- c(rep(pos1, 2*n/5), rep(mean(rng), n/5), rep(pos2, 2*n/5))
    y <- c(rep(pos1, n/5), rep(pos2, n/5), rep(mean(rng), n/5),
           rep(pos1, n/5), rep(pos2, n/5))
    list(x = x, y = y)
}

## all together
exFuns <- list(sine = genSine, nonco = genNonco, lines = genTwoLines,
               x = genX, ellipse = genEllipse, four = genFour,
               five = genFive)
funcFuns <- list(sine = genSine, line = genLine, parabola = genPar,
                 exponential = genExponential, cubic = genCubic)
                
## simulate #########################################################
## examples
set.seed(540192)
examples <- lapply(exFuns, simPattern, sds = sds, n = n, nsim = nsim,
                   splitter = sqrfySplit)
## the functional ones
set.seed(302101293)
functions <- lapply(funcFuns, simFunc, sds = sds, n = n, nsim = nsim,
                    splitter = sqrfySplit)

## plot 
br <- 100
pal <- brewer.pal(7, "Pastel2")
pchs <- 0:7
names(pal) <- names(pchs) <- c("sine", "lines", "nonco", "x", "ellipse",
                               "four", "five")
png("performance.png", width = 6, height = 5, units = "in",
    bg = "transparent", res = 480)
par(mar = c(5.1, 4.1, 1.1, 5.1))
with(examples[[names(pal)[1]]], 
     medianPlot(1 - rsq, exp(pval), xlab = expression(1~`-`~R^2), 
                breaks = br, col = pal["sine"], bg = "transparent",
                ylab = expression(Median~p), ylim = c(0,1),
                xlim = c(0,1), pch = pchs["sine"]))
for (nm in names(pal)[-1]) {
    with(examples[[nm]],
         medianPlot(1 - rsq, exp(pval), breaks = br, add = TRUE,
                    col = pal[nm], pch = pchs[nm]))
}
legend(x = "topright", legend = names(pal), col = pal, pch = pchs,
       bty = "n", inset = c(-0.2, 0), xpd = TRUE)
dev.off()

## the same thing on the log scale with a different range
br <- seq(0, 0.2, length.out = 50)
png("logPperformance.png", width = 6, height = 5, units = "in",
    bg = "transparent", res = 480)
par(mar = c(5.1, 4.1, 1.1, 5.1))
with(examples[[names(pal)[1]]],
     medianPlot(1 - rsq, pval/log(10), xlab = expression(1~`-`~R^2), 
                breaks = br, col = pal["sine"], bg = "transparent",
                ylab = expression(Median~log[10]~p), xlim = c(0, 0.2),
                pch = pchs["sine"], ylim = c(-250, 0)))
for (nm in names(pal)[-1]) {
    with(examples[[nm]], 
         medianPlot(1 - rsq, pval/log(10), breaks = br, add = TRUE,
                    col = pal[nm], pch = pchs[nm]))
}
legend(x = "topright", legend = names(pal), col = pal, pch = pchs,
       bty = "n", inset = c(-0.2, 0), xpd = TRUE)
dev.off()

## plot some example patterns
genFns <- list(sine = genSine, lines = genTwoLines, nonco = genNonco,
               x = genX, ellipse = genEllipse, four = genFour, 
               five = genFive)
png("patternExamples.png", width = 6, height = 9, units = "in",
    bg = "transparent", res = 480)
par(mfrow = c(7,5))
for (ii in seq_along(genFns)) {
  for (jj in rev(c(0.2, 0.5, 0.8, 0.9, 0.98))) {
    margs <- c(0.1, 0.1, 0.1, 0.1)
    if (ii == 1) margs[3] <- 1.1
    if (jj == 0.98) margs[2] <- 1.1
    par(mar = margs)
    plotEx(genFns[[ii]], examples[[names(genFns)[ii]]], rsq = jj,
           xaxt = "n", yaxt = "n", col = pal[ii], pch = pchs[ii],
           bg = "transparent")
    if (ii == 1) mtext(side = 3, line = 0, text = 1 - jj)
    if (jj == 0.98) mtext(side = 2, line = 0, text = names(genFns)[ii])
  }
}
dev.off()

## give the same treatment to the functional patterns
br <- 100
pal <- brewer.pal(4, "Pastel1")
pchs <- 8:11
names(pal) <- names(pchs) <- c("line", "parabola",
                               "exponential", "cubic")
png("funcPerf.png", width = 6, height = 5, units = "in",
    bg = "transparent", res = 480)
par(mar = c(5.1, 4.1, 1.1, 6.1))
with(functions[[names(pal)[1]]], 
     medianPlot(1 - rsq, exp(pval), xlab = expression(1~`-`~R^2), 
                breaks = br, col = pal["line"], bg = "transparent",
                ylab = expression(Median~p), ylim = c(0,1),
                xlim = c(0,1), pch = pchs["line"]))
for (nm in names(pal)[-1]) {
    with(functions[[nm]],
         medianPlot(1 - rsq, exp(pval), breaks = br, add = TRUE,
                    col = pal[nm], pch = pchs[nm]))
}
legend(x = "topright", legend = names(pal), col = pal, pch = pchs,
       bty = "n", inset = c(-0.3, 0), xpd = TRUE)
dev.off()

br <- seq(0, 0.2, length.out = 51)
png("logPperfFun.png", width = 6, height = 5, units = "in",
    bg = "transparent", res = 480)
par(mar = c(5.1, 4.1, 1.1, 6.1))
with(functions[[names(pal)[1]]],
     medianPlot(1 - rsq, pval/log(10), xlab = expression(1~`-`~R^2), 
                breaks = br, col = pal[1], bg = "transparent",
                ylab = expression(Median~log[10]~p), xlim = c(0, 0.2),
                pch = pchs[1], ylim = c(-700, 0)))
for (nm in names(pal)[-1]) {
    with(functions[[nm]], 
         medianPlot(1 - rsq, pval/log(10), breaks = br, add = TRUE,
                    col = pal[nm], pch = pchs[nm]))
}
legend(x = "topright", legend = names(pal), col = pal, pch = pchs,
       bty = "n", inset = c(-0.3, 0), xpd = TRUE)
dev.off()

genFns <- list(line = genLine, parabola = genPar,
               cubic = genCubic, exponential = genExponential)
png("patternFunctions.png", width = 6, height = 5, units = "in",
    bg = "transparent", res = 480)
par(mfrow = c(4,5))
for (ii in seq_along(genFns)) {
  for (jj in rev(c(0.2, 0.5, 0.8, 0.9, 0.98))) {
    margs <- c(0.1, 0.1, 0.1, 0.1)
    if (ii == 1) margs[3] <- 1.1
    if (jj == 0.98) margs[2] <- 1.1
    par(mar = margs)
    plotEx(genFns[[ii]], functions[[names(genFns)[ii]]], rsq = jj,
           xaxt = "n", yaxt = "n", col = pal[names(genFns)[ii]],
           pch = pchs[names(genFns)[ii]], bg = "transparent")
    if (ii == 1) mtext(side = 3, line = 0, text = 1 - jj)
    if (jj == 0.98) mtext(side = 2, line = 0, text = names(genFns)[ii])
  }
}
dev.off()
