##' stoppers: accept a list of bins and produce TRUE if the bin should
##' be kept for the final bin list and FALSE if the bin should be
##' split further

## potentially useful...
##' make criteria by capturing expressions
makeCriteria <- function(...) {
    cl <- match.call() # capturing inputs
    crits <- as.list(cl) # change to a list
    ## remove self reference, collapse into single OR
    paste(sapply(crits[-1], deparse), collapse = " | ")
}

##' take critera, apply to all bins in list
stopper <- function(binList, criteria) {
    sapply(binList,
           function(b) eval(parse(test = critera), envir = b))
}

## simple helper to check bin against stop criteria
#checkStop <- function(bin) {
#    eval(parse(text = bin$criteria), envir = bin)
#}

## helper to remove boundary conditions
#removeBoundaries <- function(expr) { # helper to check expressions
#    parse(text = gsub("([\\<\\>])\\=", "\\1", expr)) # remove =
#}

## one to check splits
#checkSplits <- function(blw, abv, criteria) {
#    eval(removeBoundaries(criteria), envir = blw) |
#        eval(removeBoundaries(criteria), envir = abv)
#}
