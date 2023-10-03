## PACKAGES ##########################################################
## recursive binning package
library(marbR)
library(quantreg)



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



## REAL DATA EXAMPLE #################################################
## S&P500 data: "SP500" demo in "zenplots" package, code from Marius
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
system.time({for (ii in seq_len(ncol(spPairs))) { ## ~57 mins
        pair <- spPairs[, ii] # indices of pairs
        spBins[[ii]] <- binner(spRanks[, pair[1]], spRanks[, pair[2]],
                               stopper = stopFn, splitter = chiSplit)
        if (msgInd[ii]) cat("\r Completed ", ii, " pairs")
	}
})
## drop points for smaller storage size
spBinsNP <- lapply(spBins, dropBinPoints)
## save binnings
if (writeout) saveRDS(spBinsNP,
                      file = paste0("sp500bins", "NoPts.Rds"))

## load pre-processed data
spBinsNP <- readRDS("sp500binsNoPts.Rds")
## get chi statistics across the bins
spChis <- lapply(spBinsNP, function(bns) binChi(bns))
spChiStats <- sapply(spChis, function(x) x$stat)
spChiResid <- sapply(spChis, function(x) x$residuals)
spChiNbin <- sapply(spChiResid, length)
## order by most interesting
spOrd <- order(spChiStats, decreasing = TRUE)
spMaxRes <- max(abs(unlist(spChiResid)))

## quantiles of the statistics to regress
qnts <- c(0.95, 0.99, 0.999)
## add the null density
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
spCol <- rgb(spRGB[,1], spRGB[,2], spRGB[,3])
## perform quantile regression for the null data as well
modQnt <- rq(chi ~ nbin, tau = c(0.95, 0.99, 0.999),
             data = data.frame(nbin = c(nulls$chiSplit["nbin",,]),
                               chi = c(nulls$chiSplit["chi",,])))
predQnt <- predict(modQnt,
                   newdata = data.frame(nbin = binVals))

## plot the sp500 point cloud alongside the null
png("sp500vsNullPoints.png", width = 3.5, height = 3.5, units = "in",
    res = 480)
narrowPlot(xgrid = seq(1, 2, by = 0.25),
           ygrid = seq(1, 3, by = 0.5), ylim = c(1, 3.2),
           xlab = expression(log[10]~"(Number of bins)"),
           ylab = expression(log[10]~{"("~chi^2~statistic~")"}))
points(log(nulls$chiSplit["nbin", c(5, 6),], 10),
       log(nulls$chiSplit["chi", c(5, 6),], 10), cex = 1,
       pch = 20, col = adjustcolor("steelblue", 0.03))
points(log(spChiNbin, 10), log(spChiStats, 10),
       col = adjustcolor("firebrick", 0.03),
       pch = 20)
for (ii in 1:3) {
    lines(log(binVals[-c(1,2)], 10),
          log(predQnt[-c(1,2),ii], 10), lty = ii)
    yadj <- 1 - 0.5*(ii - 1)
    text(labels = paste0(100*qnts[ii], "%"),
         x = log(binVals[3], 10), y = log(predQnt[3, ii], 10),
         adj = c(1, yadj), cex = 0.6)
}
legend(x = "bottomright", cex = 0.8, legend = c("Null", "S&P500"),
       pch = 20, col = c("steelblue", "firebrick"))
dev.off()

## plot the sp500 point cloud coloured by empirical p-value with the
## quantile regression lines alongside
png("sp500empPColour.png", width = 3.5, height = 3.5, units = "in",
    res = 480)
narrowPlot(xgrid = seq(1, 2, by = 0.25),
           ygrid = seq(1.5, 3, by = 0.5), ylim = c(1.4, 3.2),
           xlab = expression(log[10]~"(Number of bins)"),
           ylab = expression(log[10]~{"("~chi^2~statistic~")"}),
           mars = c(2.1, 2.1, 2.1, 2.1))
points(log(spChiNbin, 10), log(spChiStats, 10), cex = 0.5,
       col = adjustcolor(spCol, 0.2), pch = 20)
for (ii in 1:3) {
    lines(log(binVals[-c(1,2)], 10),
          log(predQnt[-c(1,2),ii], 10), lty = ii)
    yadj <- 1 - 0.5*(ii - 1)
    text(labels = paste0(100*qnts[ii], "%"),
         x = log(binVals[3], 10), y = log(predQnt[3, ii], 10),
         adj = c(1, yadj), cex = 0.6)
}
addMarHists(log(spChiNbin, 10), log(spChiStats, 10),
            xcuts = seq(1, 2, by = 0.03125),
            ycuts = seq(1.5, 3, by = 0.0625))
dev.off()

## use this to plot the top pairs and their binnings
png("sp500top36.png", width = 5, height = 5, units = "in",
    res = 480)
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
    plotBinning(spBinsNP[[prInd]],
                fill = residualFill(spBinsNP[[prInd]],
                                    maxRes = spMaxRes),
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
