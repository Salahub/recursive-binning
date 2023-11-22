## PACKAGES ##########################################################
## recursive binning package
library(AssocBin)


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

## a custom plot with thinner margins
narrowPlot <- function(xgrid, ygrid, main = "", xlab = "", ylab = "",
                       xticks = xgrid, yticks = ygrid,
                       xlim = range(xgrid), ylim = range(ygrid),
                       addGrid = TRUE, ...) {
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


## S&P DATA EXAMPLE ##################################################
## S&P500 data: "SP500" demo in "zenplots" package, code from Marius
## Hofert, produces a set of pseudo-observations that are uniform
## these are loaded here and converted to ranks
data(sp500pseudo)
spRanks <- apply(sp500pseudo, 2, rank, ties.method = "random")
rownames(spRanks) <- NULL
## take a subset for inspection pairwise
set.seed(63212023)
spRanks <- spRanks[, sample(1:ncol(spRanks), 50)]
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
## iterate through all pairs
## ~ 20s (grows roughly linearly)
system.time({for (ii in seq_len(ncol(spPairs))) { ## ~20s
        pair <- spPairs[, ii] # indices of pairs
        spBins[[ii]] <- binner(spRanks[, pair[1]], spRanks[, pair[2]],
                               stopper = stopFn, splitter = chiSplit)
	}
})

## get chi statistics across the bins
spChis <- lapply(spBins, function(bns) binChi(bns))
spChiStats <- sapply(spChis, function(x) x$stat)
spChiResid <- sapply(spChis, function(x) x$residuals)
spChiNbin <- sapply(spChiResid, length)
## order by most interesting
spOrd <- order(spChiStats, decreasing = TRUE)
spMaxRes <- max(abs(unlist(spChiResid)))

## compute some null examples
nVar <- ncol(spRanks)
null <- sapply(1:nVar, function(ii) sample(1:nrow(spRanks)))
## the same number of columns are present in the null data here, so
## use the same pair indices
nullBins <- vector("list", ncol(spPairs))
for (ii in seq_len(ncol(spPairs))) {
        pair <- spPairs[, ii] # indices of pairs
        nullBins[[ii]] <- binner(null[, pair[1]], null[, pair[2]],
                                 stopper = stopFn, splitter = chiSplit)
}

## get chi statistics for the null examples
nullChis <- lapply(nullBins, function(bns) binChi(bns))
nullChiStats <- sapply(nullChis, function(x) x$stat)
nullChiResid <- sapply(nullChis, function(x) x$residuals)
nullChiNbin <- sapply(nullChiResid, length)
nullMaxRes <- max(abs(unlist(nullChiResid)))

## plot the sp500 point cloud alongside the null
oldPar <- par(mar = c(2.1, 2.1, 3.1, 3.1))
narrowPlot(xgrid = seq(1, 2, by = 0.25),
           ygrid = seq(1, 3, by = 0.5), ylim = c(1, 3.2),
           xlab = expression(log[10]~"(Number of bins)"),
           ylab = expression(log[10]~{"("~chi^2~statistic~")"}))
points(log(nullChiNbin, 10), log(nullChiStats, 10), cex = 1,
       pch = 20, col = adjustcolor("steelblue", 0.2))
points(log(spChiNbin, 10), log(spChiStats, 10),
       col = adjustcolor("firebrick", 0.2),
       pch = 20)
legend(x = "bottomright", cex = 0.8, legend = c("Null", "S&P500"),
       pch = 20, col = c("steelblue", "firebrick"))
addMarHists(log(spChiNbin, 10), log(spChiStats, 10),
            xcuts = seq(1, 2, by = 0.03125),
            ycuts = seq(1.5, 3, by = 0.0625))

## use this to plot the top pairs and their binnings
par(mfrow = c(6, 6), mar = c(0.1, 0.55, 1.1, 0.55))
for (prInd in spOrd[1:36]) {
    pr <- spPairs[, prInd] # pair indices
    plot(NA, xlim = c(1, (nrow(spRanks))),
         ylim = c(1, (nrow(spRanks))), axes = "F",
         xlab = colnames(spRanks)[pr[1]],
         ylab = colnames(spRanks)[pr[2]],
         main = "")
    mtext(paste(colnames(spRanks)[pr], collapse = ":"),
          cex = 0.6)
    plotBinning(spBins[[prInd]],
                fill = residualFill(spBins[[prInd]],
                                    maxRes = spMaxRes),
                add = TRUE, pch = "")
    points(spRanks[, pr[1]], spRanks[, pr[2]], pch = ".",
           col = adjustcolor("gray50"))
}
## fix parameters
par(oldPar)
