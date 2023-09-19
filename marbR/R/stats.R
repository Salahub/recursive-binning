
## wrapper functions to compute different statistics over the bins
## a chi^2 measure
binChi <- function(bins, agg = sum) {
    obs <- sapply(bins, function(bn) bn$n)
    ex <- sapply(bins, function(bn) bn$expn)
    resids <- (obs - ex)^2/ex
    signs <- sign(obs - ex) # signs of residuals
    list(residuals = signs*sqrt(resids), stat = agg(resids))
}

## a mutual information measure
binMI <- function(bins, agg = sum) {
    obs <- sapply(bins, function(bin) bin$n)
    ex <- sapply(bins, function(bin) bin$expn)
    n <- sum(obs)
    resids <- log(obs/ex)
    resids[obs == 0] <- 0
    probs <- obs/n
    list(residuals = resids, stat = agg(resids*probs))
}

## absolute values of the counts
binAbsDif <- function(bins, agg = sum) {
    obs <- sapply(bins, function(bin) bin$n)
    ex <- sapply(bins, function(bin) bin$expn)
    resids <- abs(obs - ex)
    signs <- sign(obs - ex)
    list(residuals = signs*resids, stat = agg(resids))
}
