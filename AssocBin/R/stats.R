##' Binstatistics
##' @title Statistics for bins
##' @description These functions compute statistics based on observed
##' and expected counts for a list of bins.
##' @details Three functions are provided by default, `binChi`
##' computes the chi-squared statistic by taking the squared
##' difference between observed and expected counts and dividing this
##' by the expected counts. `binMi` computes the mutual information
##' for each bin using the observed and expected counts. Finally,
##' `binAbsDif` computes the absolute difference between observed
##' and expected counts. Each function first computes a value on
##' every bin independently and stores all these values in memory
##' before using the function provided in the optional argument `agg`
##' to aggregate these values.
##' @param bins a list of bins, each a list with elements `x`, `y`,
##' `depth`, `bnds` (list with elements `x` and `y`), `expn`, `n`
##' @param agg function which is aggregates the individual statistics
##' computed over each bin
##' @return A list with elements `residuals` and `stat` reporting the
##' individual statistic values (possibly transformed) and the
##' aggegrated statistic value.
##' @examples
##' binList1 <- list(list(x = c(1,2), y = c(3,1), depth = 1, n = 2,
##'                       expn = 2),
##'                 list(x = c(3,4), y = c(2,4), depth = 1, n = 2,
##'                      expn = 2))
##' binList2 <- list(list(x = c(1,2), y = c(3,1), depth = 6, n = 2,
##'                       expn = 4),
##'                 list(x = c(), y = c(), depth = 1, n = 0, expn = 1))
##' binChi(binList1)
##' binChi(binList2)
##' binMI(binList1)
##' binMI(binList2)
##' binAbsDif(binList2)
##' @author Chris Salahub
##' @describeIn binstatistics Chi-squared statistic
binChi <- function(bins, agg = sum) {
    obs <- sapply(bins, function(bn) bn$n)
    ex <- sapply(bins, function(bn) bn$expn)
    resids <- (obs - ex)^2/ex
    signs <- sign(obs - ex) # signs of residuals
    list(residuals = signs*sqrt(resids), stat = agg(resids))
}
##' @describeIn binstatistics Mutual information
binMI <- function(bins, agg = sum) {
    obs <- sapply(bins, function(bin) bin$n)
    ex <- sapply(bins, function(bin) bin$expn)
    n <- sum(obs)
    resids <- log(obs/ex)
    resids[obs == 0] <- 0
    probs <- obs/n
    list(residuals = resids, stat = agg(resids*probs))
}
##' @describeIn binstatistics Absolute difference between observed
##' and expected
binAbsDif <- function(bins, agg = sum) {
    obs <- sapply(bins, function(bin) bin$n)
    ex <- sapply(bins, function(bin) bin$expn)
    resids <- abs(obs - ex)
    signs <- sign(obs - ex)
    list(residuals = signs*resids, stat = agg(resids))
}
