source("binMethods.R")
source("recursiveBins.R")

## go against some random data
set.seed(16062021)
randx <- sample(1:1e4)
randy <- sample(1:1e4)
randBin <- makeTree(randx, randy,
                    stopCriterion(depthLim = 10, areaLim = 5),
                    maxScoreSplit(chiScores, ties = "x"), quickBin)
randBin.mi <- makeTree(randx, randy,
                    stopCriterion(depthLim = 4, areaLim = 5),
                    maxScoreSplit(miScores, ties = "x"), quickBin)
## plot these
maxRes <- max(c(binChi(unlist(randBin))$residuals,
                binChi(unlist(randBin.mi))$residuals))
plotBinning(unlist(randBin), pch = ".",
            fill = residualFill(unlist(randBin), maxRes = maxRes))
dev.new()
plotBinning(unlist(randBin.mi), pch = ".",
            fill = residualFill(unlist(randBin.mi), maxRes = maxRes))

## a straight line
linex <- 1:1e4
liney <- 1:1e4
lineBin <- makeTree(linex, liney,
                    stopCriterion(depthLim = 10, areaLim = 5),
                    maxScoreSplit(chiScores, ties = "x"), quickBin)
lineBin.mi <- makeTree(linex, liney,
                    stopCriterion(depthLim = 1, areaLim = 5),
                    maxScoreSplit(miScores, ties = "x"), quickBin)
## plot these
maxRes <- max(c(binChi(unlist(lineBin))$residuals,
                binChi(unlist(lineBin.mi))$residuals))
plotBinning(unlist(lineBin), pch = ".",
            fill = residualFill(unlist(lineBin), maxRes = maxRes))
dev.new()
plotBinning(unlist(lineBin.mi), pch = ".",
            fill = residualFill(unlist(lineBin.mi), maxRes = maxRes))

## how does changing the depth limit impact this?
set.seed(506391)
nsim <- 1e3
depths <- 2:10
simDataSets <- replicate(nsim, data.frame(x = sample(1:1e4),
                                          y = sample(1:1e4)),
                         simplify = FALSE)
depthSeq.chi <- array(NA, dim = c(4,length(depths),nsim))
depthSeq.mi <- array(NA, dim = c(4,length(depths),nsim))
for (ii in 1:nsim) {
    depthSeq.chi[,,ii] <- sapply(2:10,
                                 function(depth) {
                                     tr <- unlist(makeTree(simDataSets[[ii]]$x,
                                                           simDataSets[[ii]]$y,
                                                           stopCriterion(depthLim = depth,
                                                                         areaLim = 10),
                                                           maxScoreSplit(chiScores, ties = "x"),
                                                           quickBin))
                                     c(chi = binChi(tr)$stat, mi = binMI(tr)$stat,
                                       nbin = length(tr), maxDep = max(sapply(tr, attr, "depth")))
                                 })
    depthSeq.mi[,,ii] <- sapply(2:10,
                                function(depth) {
                                    tr <- unlist(makeTree(simDataSets[[ii]]$x,
                                                          simDataSets[[ii]]$y,
                                                          stopCriterion(depthLim = depth,
                                                                        areaLim = 10),
                                                          maxScoreSplit(miScores, ties = "x"),
                                                          quickBin))
                                    c(chi = binChi(tr)$stat, mi = binMI(tr)$stat,
                                      nbin = length(tr), maxDep = max(sapply(tr, attr, "depth")))
                                })
    if (ii %% 50 == 0) cat(paste("\r Simulated data set:", ii))
}
dimnames(depthSeq.chi) <- list(c("chi", "mi", "nbin", "maxDep"))
dimnames(depthSeq.mi) <- list(c("chi", "mi", "nbin", "maxDep"))

## plot densities
nbinSeq <- sort(unique(c(depthSeq.chi["nbin",,])))
depthDens <- lapply(1:length(depths),
                    function(d) density(depthSeq.chi["chi",d,]))
plot(NA, xlim = range(depthSeq.chi["chi",,]),
     ylim = range(0,1), xlab = expression(chi^2~"statistic"),
     ylab = "Scaled density")
for (ii in c(1,3,5,7,9)) {
    tempDens <- depthDens[[ii]]
    tempDens$y <- tempDens$y/max(tempDens$y, na.rm = TRUE)
    polygon(tempDens, col = adjustcolor(hcl.colors(9, "Dark 2"), 0.2)[ii])
}
for (ii in c(1,3,5,7,9)) {
    lines(x = rep(mean(depthSeq.chi["nbin",ii,]),2), y = c(0,0.2), lwd = 2,
          col = hcl.colors(9, "Dark 2")[ii])
}

## increased depth increases statistic values
plot(depthSeq.chi["nbin",,], depthSeq.chi["chi",,], xlab = "Number of bins",
     ylab = "", col = adjustcolor(hcl.colors(9, "Dark 2"), 0.5),
     pch = 20)
mtext(expression(chi^2~statistic), side = 2, line = 2)
abline(v = seq(0, 200, by = 50), col = "gray50", lty = 2)
abline(h = seq(0, 1500, by = 500), col = "gray50", lty = 2)
for (ii in 1:9) points(mean(depthSeq.chi["nbin",ii,]),
                       mean(depthSeq.chi["chi",ii,]), col = "black",
                       bg = hcl.colors(9, "Dark 2")[ii], pch = 22)
for (p in c(0.01)) {
    lines(2:251, qchisq(1-p, 1:250), lty = 3)
}
legend(x = "bottomright", legend = 2:10, title = "Max. Depth",
       bg = "white", fill = adjustcolor(hcl.colors(9, "Dark 2"), 0.5))

## try p-values?
plot(rbind(sqrt(depthSeq.chi["nbin",,]),NA),
     rbind(pchisq(depthSeq.chi["chi",,], depthSeq.chi["nbin",,]-1,
                  lower.tail = FALSE),NA),
     type = "b", col = adjustcolor("steelblue", 0.5))

## for both chi square and MIC
plot(depthSeq.chi["nbin",,], depthSeq.chi["mi",,], xlab = "Number of bins",
     ylab = "", col = adjustcolor(hcl.colors(9, "Dark 2"), 0.5), pch = 20)
mtext("Mutual information", side = 2, line = 3)
for (ii in 1:9) points(mean(depthSeq.chi["nbin",ii,]),
                       mean(depthSeq.chi["mi",ii,]), col = "black",
                       bg = hcl.colors(9, "Dark 2")[ii], pch = 22)
abline(v = seq(0, 200, by = 50), col = "gray50", lty = 2)
abline(h = seq(0, 0.05, by = 0.01), col = "gray50", lty = 2)
legend(x = "bottomright", legend = 2:10, title = "Max. Depth",
       bg = "white", fill = adjustcolor(hcl.colors(9, "Dark 2"), 0.5))

## square root scaling
plot(sqrt(depthSeq.chi["nbin",,]), sqrt(depthSeq.chi["chi",,]),
     xlab = "Square root of number of bins",
     ylab = "", col = adjustcolor(hcl.colors(9, "Dark 2"), 0.5), pch = 20)
mtext(expression(sqrt(chi^2)~statistic), side = 2, line = 2)
for (ii in 1:9) points(sqrt(mean(depthSeq.chi["nbin",ii,])),
                       sqrt(mean(depthSeq.chi["chi",ii,])), col = "black",
                       bg = hcl.colors(9, "Dark 2")[ii], pch = 22)
abline(v = seq(2, 14, by = 2), col = "gray50", lty = 2)
abline(h = seq(0, 40, by = 10), col = "gray50", lty = 2)
legend(x = "bottomright", legend = 2:10, title = "Max. Depth",
       bg = "white", fill = adjustcolor(hcl.colors(9, "Dark 2"), 0.5))
## robust linear fit for asymmetric errors
library(MASS)
simpleFit <- rlm(chi ~ nb, data = data.frame(nb = sqrt(c(depthSeq.chi["nbin",,])),
                                             chi = sqrt(c(depthSeq.chi["chi",,]))),
                 psi = psi.bisquare)
abline(simpleFit$coef[1], simpleFit$coef[2])

## of mutual information
plot(sqrt(depthSeq.chi["nbin",,]), sqrt(depthSeq.chi["mi",,]),
     xlab = "Square root of number of bins",
     ylab = "Square root of mutual information",
     col = adjustcolor(hcl.colors(9, "Dark 2"), 0.5), pch = 20)
for (ii in 1:9) points(sqrt(mean(depthSeq.chi["nbin",ii,])),
                       sqrt(mean(depthSeq.chi["mi",ii,])), col = "black",
                       bg = hcl.colors(9, "Dark 2")[ii], pch = 22)
abline(v = seq(2, 14, by = 2), col = "gray50", lty = 2)
abline(h = seq(0, 0.2, by = 0.05), col = "gray50", lty = 2)

## compare the mutual information between chi/mi splits
plot(depthSeq.mi["mi",,], depthSeq.chi["mi",,],
     xlab = "Mutual information score function",
     ylab = "",
     col = adjustcolor(hcl.colors(9, "Dark 2"), 0.5), pch = 20)
mtext(expression(chi^2~"score function"), side = 2, line = 2)
legend(x = "bottomright", legend = 2:10, title = "Max. Depth",
       bg = "white", col = adjustcolor(hcl.colors(9, "Dark 2"), 0.5),
       pch = 20)
abline(0, 1, lty = 2)
plot(depthSeq.mi["chi",,], depthSeq.chi["chi",,],
     xlab = "Mutual information score function",
     ylab = "",
     col = adjustcolor(hcl.colors(9, "Dark 2"), 0.5), pch = 20)
mtext(expression(chi^2~"score function"), side = 2, line = 2)
legend(x = "bottomright", legend = 2:10, title = "Max. Depth",
       bg = "white", col = adjustcolor(hcl.colors(9, "Dark 2"), 0.5),
       pch = 20)
abline(0, 1, lty = 2)

## p-values
plot(pchisq(depthSeq.mi["chi",,], depthSeq.chi["nbin",,]-1,
            lower.tail = FALSE),
     pchisq(depthSeq.chi["chi",,], depthSeq.chi["nbin",,]-1,
              lower.tail = FALSE),
     xlab = "Mutual information score function",
     ylab = "",
     col = adjustcolor(hcl.colors(9, "Dark 2"), 0.5), pch = 20)
mtext(expression(chi^2~"score function"), side = 2, line = 2)
legend(x = "bottomright", legend = 2:10, title = "Max. Depth",
       bg = "white", col = adjustcolor(hcl.colors(9, "Dark 2"), 0.5),
       pch = 20)
abline(0, 1, lty = 2)


## plot sequential position coloured by initial quintile
plot(c(rbind(sqrt(depthSeq.chi["nbin",,]), NA)),
     c(rbind(sqrt(depthSeq.chi["chi",,]), NA)),
     type = 'n', xlab = "Square root of number of bins", ylab = "")
mtext(expression(sqrt(chi^2)~statistic), side = 2, line = 2)
lineCols <-  adjustcolor(hcl.colors(5, "Dark 3"),
                         0.5)[as.numeric(cut(sqrt(depthSeq.chi["chi",1,]),
                                             quantile(sqrt(depthSeq.chi["chi",1,]),
                                                      c(0,0.2,0.4,0.6,0.8,1))))]
for (ii in 1:nsim) {
    lines(sqrt(depthSeq.chi["nbin",,ii]),
          sqrt(depthSeq.chi["chi",,ii]),
          col = lineCols[ii], type = "b")
}

## ... and for the MIC
plot(c(rbind(sqrt(depthSeq.chi["nbin",,]), NA)),
     c(rbind(sqrt(depthSeq.chi["mi",,]), NA)),
     type = 'n', xlab = "Square root of number of bins",
     ylab = "Square root of mutual information")
lineCols <-  adjustcolor(hcl.colors(5, "Dark 3"),
                         0.5)[as.numeric(cut(sqrt(depthSeq.chi["mi",1,]),
                                             quantile(sqrt(depthSeq.chi["mi",1,]),
                                                      c(0,0.2,0.4,0.6,0.8,1))))]
for (ii in 1:nsim) {
    lines(sqrt(depthSeq.chi["nbin",,ii]),
          sqrt(depthSeq.chi["mi",,ii]),
          col = lineCols[ii])
}


## abalone data
url <- "https://archive.ics.uci.edu/ml/machine-learning-databases/abalone/"
## loading here is locally, otherwise replace "~/Downloads/" with the url above
abalone <- read.table("~/Downloads/abalone.data", sep = ",")
abalone.names <- scan("~/Downloads/abalone.names", what = character(),
                      sep = "\n")[64:72]
names(abalone) <- regmatches(gsub(" ", ".", abalone.names),
                             regexpr("(?<=^\\t)[A-Za-z\\.]+(?=\\t)",
                                     gsub(" ", ".", abalone.names),
                                     perl = TRUE))
## take ranks, get all pairs
set.seed(2306201)
abPairs <- combn(ncol(abalone), 2)
abRanks <- lapply(abalone, rank, ties.method = "random")
abBins <- lapply(1:ncol(abPairs),
                 function(p) {
                     lapply(depths,
                            function(dep) {
                     makeTree(abRanks[[abPairs[1,p]]],
                              abRanks[[abPairs[2,p]]],
                              stopCriterion(depthLim = dep,
                                            areaLim = 10),
                              maxScoreSplit(chiScores), quickBin)
                 })})
names(abBins) <- paste(names(abRanks)[abPairs[1,]],
                       names(abRanks)[abPairs[2,]],
                       sep = "-")
abChis <- lapply(abBins, function(ds) lapply(ds, function(bnlst) binChi(unlist(bnlst))))
abMis <- lapply(abBins, function(ds) lapply(ds, function(bnlst) binMI(unlist(bnlst))))

## plot these against other points
plot(depthSeq.chi["nbin",,], depthSeq.chi["chi",,], xlab = "Number of bins",
     ylab = "", col = adjustcolor("gray", 0.5), pch = 20,
     ylim = c(0, max(sapply(abChis, function(lst) sapply(lst, function(scr) scr$stat)))),
     xlim = c(0, 275))
mtext(expression(chi^2~statistic), side = 2, line = 2)
lines(c(rbind(sapply(abChis[grepl("Rings", names(abChis))],
                     function(lst) sapply(lst, function(scr) length(scr$residuals))),NA)),
      c(rbind(sapply(abChis[grepl("Rings", names(abChis))],
                     function(lst) sapply(lst, function(scr) scr$stat)), NA)),
      pch = 20, col = adjustcolor("#fc8d62", 0.8), type = "b")
lines(c(rbind(sapply(abChis[grepl("Sex", names(abChis))],
                     function(lst) sapply(lst, function(scr) length(scr$residuals))),NA)),
      c(rbind(sapply(abChis[grepl("Sex", names(abChis))],
                     function(lst) sapply(lst, function(scr) scr$stat)), NA)),
      pch = 20, col = adjustcolor("#66c2a5", 0.8), type = "b")
lines(c(rbind(sapply(abChis[!grepl("Rings|Sex", names(abChis))],
                     function(lst) sapply(lst, function(scr) length(scr$residuals))),NA)),
      c(rbind(sapply(abChis[!grepl("Rings|Sex", names(abChis))],
                     function(lst) sapply(lst, function(scr) scr$stat)), NA)),
      pch = 20, col = adjustcolor("#8da0cb", 0.8), type = "b")
legend(x = "topleft", pch = 20, legend = c("Continuous", "Categorical", "Count", "Null"),
       col = adjustcolor(c("#8da0cb","#66c2a5","#fc8d62","gray"), 0.8),
       title = "Data type")

## ... also using MI
plot(depthSeq.chi["nbin",,], depthSeq.chi["mi",,], xlab = "Number of bins",
     ylab = "Mutual information", col = adjustcolor("gray", 0.5), pch = 20,
     ylim = c(0, max(sapply(abMis, function(lst) sapply(lst, function(scr) scr$stat)))),
     xlim = c(0, 275))
lines(c(rbind(sapply(abMis[grepl("Rings", names(abMis))],
                     function(lst) sapply(lst, function(scr) length(scr$residuals))),NA)),
      c(rbind(sapply(abMis[grepl("Rings", names(abMis))],
                     function(lst) sapply(lst, function(scr) scr$stat)), NA)),
      pch = 20, col = adjustcolor("#fc8d62", 0.8), type = "b")
lines(c(rbind(sapply(abMis[grepl("Sex", names(abMis))],
                     function(lst) sapply(lst, function(scr) length(scr$residuals))),NA)),
      c(rbind(sapply(abMis[grepl("Sex", names(abMis))],
                     function(lst) sapply(lst, function(scr) scr$stat)), NA)),
      pch = 20, col = adjustcolor("#66c2a5", 0.8), type = "b")
lines(c(rbind(sapply(abMis[!grepl("Rings|Sex", names(abMis))],
                     function(lst) sapply(lst, function(scr) length(scr$residuals))),NA)),
      c(rbind(sapply(abMis[!grepl("Rings|Sex", names(abMis))],
                     function(lst) sapply(lst, function(scr) scr$stat)), NA)),
      pch = 20, col = adjustcolor("#8da0cb", 0.8), type = "b")
legend(x = "topleft", pch = 20, legend = c("Continuous", "Categorical", "Count", "Null"),
       col = adjustcolor(c("#8da0cb","#66c2a5","#fc8d62","gray"), 0.8),
       title = "Data type")

## get some pairwise values
linex <- 1:(nrow(abalone))
liney <- 1:(nrow(abalone))
lineBin <- makeTree(linex, liney, # to scale values
                    stopCriterion(depthLim = 6, areaLim = 10),
                    maxScoreSplit(chiScores, ties = "x"), quickBin)
lineChi <- binChi(unlist(lineBin))
pairwiseAb <- matrix(lineChi$stat, ncol(abalone), ncol(abalone))
for (ii in 1:ncol(abPairs)) {
    tempC <- abPairs[,ii]
    tempChi <- abChis[[ii]]
    pairwiseAb[tempC[1], tempC[2]] <- pairwiseAb[tempC[2], tempC[1]] <- tempChi[[5]]$stat
}
pairwiseAb <- pairwiseAb/lineChi$stat


## definitely some stronger relationships than uniform...
## also, note the strong ordering for both, plot all of these
abOrders <- sapply(1:length(depths),
                   function(ii){
                       order(sapply(abChis, function(scr) scr[[ii]]$stat),
                             decreasing = TRUE)
                   })
maxRes <- sapply(1:length(depths),
                 function(ii) {
                     max(abs(unlist(sapply(abChis, function(scr) scr[[ii]]$residuals))))
                 })
ind <- 5
par(mfrow = c(6,6), mar = c(1.1,1.1,2.1,1.1))
for (ii in 1:length(abOrders[,ind])) {
    tempBins <- unlist(abBins[abOrders[,ind]][[ii]][[ind]])
    plotBinning(tempBins, pch = "", xaxt = "n", yaxt = "n",
                residualFill(tempBins, maxRes = maxRes[ind]),
                xaxt = "n", yaxt = "n",
                main = names(abBins)[abOrders[,ind]][ii])
}

## the iris data
data(iris)
wid <- rank(iris$Sepal.Width, ties.method = "random")
len <- rank(iris$Sepal.Length, ties.method = "random")
widLen <- makeTree(wid, len,
                   stopCriterion(depthLim = 5, areaLim = 5),
                   maxScoreSplit(chiScores, ties = "x"), quickBin)
widLen.mi <- makeTree(wid, len,
                   stopCriterion(depthLim = 5, areaLim = 5),
                   maxScoreSplit(miScores, ties = "x"), quickBin)

## some simulated structural data
data <- data.frame(x = c(sample(rnorm(2000, c(-10, 10), c(4,3))), rnorm(400)),
                   y = c(sample(rnorm(2000, c(-10, 10), c(1,5))), rnorm(400)))
fourCent <- makeTree(rank(data$x), rank(data$y),
                     stopCriterion(depthLim = 5, areaLim = 5),
                     maxScoreSplit(chiScores, ties = "x"), quickBin)

## upper corner
upCorn <- list()
upCorn$x <- runif(1000)
upCorn$x <- (1 - sqrt(upCorn$x))*1000
upCorn$y <- upCorn$x + runif(1000)*(1000 - upCorn$x)
upCornChi <- makeTree(upCorn$x, upCorn$y,
                     stopCriterion(depthLim = 5, areaLim = 5),
                     maxScoreSplit(chiScores, ties = "x"), quickBin)
