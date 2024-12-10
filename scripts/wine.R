## TODO: load other packages to compare?
library(AssocBin)
library(BET)

## read in wine data
vars <- scan("./winequality-red.csv", what = character(), nlines = 1,
             sep = ";")
reds <- matrix(scan("./winequality-red.csv", what = numeric(),
                    sep = ";", skip = 1),
               ncol = length(vars), byrow = TRUE)
whites <- matrix(scan("./winequality-white.csv", what = numeric(),
                      sep = ";", skip = 1), ncol = length(vars),
                 byrow = TRUE)

## process into a data frame
wines <- data.frame(rbind(reds, whites),
                    type = factor(c(rep("red", nrow(reds)),
                                    rep("white", nrow(whites)))))
colnames(wines) <- c("fixed.acid", "volatile.acid", "citric.acid",
                     "resid.sug", "chlorides", "free.SO2",
                     "total.SO2", "density", "pH", "sulphate",
                     "alcohol", "quality", "type")
## introduce a new categorical variable
wines$alcohol <- cut(wines$alcohol, breaks = c(0, 9, 12, 20),
                     labels = c("low", "med", "hi"),
                     include.lowest = TRUE)
## address small counts in the tails
wines$quality <- cut(wines$quality, breaks = c(0, 4.5, 5.5, 6.5,
                                               7.5, 10),
                     labels = c("4-","5","6","7","8+"))

## look at dependence
set.seed(7042024)
stopCrits <- makeCriteria(depth > 10, n < 1, expn <= 10)
sqrRSplt <- function(bn) rIntSplit(bn, squarify = TRUE)
wineDep <- inDep(wines, stopCriteria = stopCrits, conCon = sqrRSplt)
summary(wineDep)
png("wineTriplets.png", width = 4, height = 6, units = "in", res = 480,
    bg = "transparent")
plot(wineDep, pch = ".", border = NA, bg = "transparent",
     which = c(1, 2, 3, 4), buffer = 0.02)
dev.off()

sqrMaxSplit <- function(bn) maxScoreSplit(bn, chiScores, squarify = TRUE)
uniMaxSplit <- function(bn) uniMaxScoreSplit(bn, chiScores)
wineDepMax <- inDep(wines, stopCriteria = stopCrits, conCon = sqrMaxSplit,
                    catCon = uniMaxSplit)
summary(wineDepMax)
plot(wineDepMax, pch = ".", border = NA)

## compare between the two
prs <- wineDep$pairs
plot(wineDep$logps[prs], wineDepMax$logps[prs])
plot(rank(wineDep$logps[prs]), rank(wineDepMax$logps[prs]))
abline(a = 0, b = 1)
