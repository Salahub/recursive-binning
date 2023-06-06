## wrapper functions to compute different statistics over the bins
## a chi^2 measure
binChi <- function(bins, agg = sum) {
    feats <- getObsExp(bins) # observed and expected
    resids <- with(feats, (ob-ex)^2/ex)
    signs <- with(feats, sign(ob-ex)) # signs of points
    list(residuals = signs*sqrt(resids), stat = agg(resids))
}

## a mutual information measure
binMI <- function(bins, agg = sum) {
    feats <- getObsExp(bins)
    resids <- with(feats, log(ob/ex))
    resids[feats$ob == 0] <- 0 # remove NAs
    probs <- with(feats, ob/n)
    list(residuals = resids, stat = agg(resids*probs))
}

## simply sum the absolute values of the counts
binCount <- function(bins, agg = sum) {
    feats <- getObsExp(bins)
    resids <- with(feats, abs(ob - ex))
    signs <- with(feats, sign(ob-ex))
    list(residuals = signs*resids, stat = agg(resids))
}
