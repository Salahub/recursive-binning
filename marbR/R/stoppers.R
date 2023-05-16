## bin size stopper
sizeStop <- function(binList, lim = 100) {
  lens <- sapply(binList, function(bin) bin$n) # bin sizes
  lens < lim # check against limit
}
