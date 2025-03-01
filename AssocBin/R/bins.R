##' @title Make a bin
##' @description Creating a new bin object
##' @details `makeBin` creates a bin list based on the arguments
##' provided to it. Should some be missing, basic defaults ensure
##' that the complete set of bin characteristics are created in
##' the resulting list representing the bin object.
##' @param x numeric vector of observations on the first variable
##' @param y numeric vector of observations on the second variable
##' @param bnds list of length two with named elements `x` and `y`
##' each a vector of length two giving respective bin boundaries
##' @param expn expected number of points in the bin, can be
##' non-integer
##' @param n observed count of points in the bin
##' @param depth number of splits from the initial bin to the bin
##' @param stopped logical; should the bin be split further?
##' @return A list with named elements matching these arguments
##' @examples
##' makeBin(x = 1:10, y = sample(1:10),
##' bnds = list(x = c(0,10), y = c(0, 10)), expn = 10, n = 10,
##' depth = 0, stopped = FALSE)
##' @author Chris Salahub
makeBin <- function(x, y,
                    bnds = list(x = range(x) - c(1,0),
                                y = range(y) - c(1,0)),
                    expn = length(x), n = length(x), depth = 0,
                    stopped = FALSE) {
    list(x = x, y = y, bnds = bnds, expn = expn, n = n,
         depth = depth, stopped = stopped)
}
