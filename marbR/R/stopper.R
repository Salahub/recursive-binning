##' @title Make stop crteria
##' @description Capture a sequence of logical statements and append
##' them into a single expression.
##' @details This function, along with `stopper` dictates the stop
##' behaviour of recursive binning. It accepts an arbitrary number
##' of arguments, each a logical statement, and appends them all into
##' a string separated by the pipe character.
##' @param ... an arbitrary number of expressions which evaluate to
##' logicals
##' @return A string which appends all expressions together.
##' @examples
##' @author Chris Salahub
makeCriteria <- function(...) {
    cl <- match.call() # capturing inputs
    crits <- as.list(cl) # change to a list
    ## remove self reference, collapse into single OR
    paste(sapply(crits[-1], deparse), collapse = " | ")
}

##' @title Check bins against stop criteria
##' @description Evaluate the stop `criteria` for each bin in
##' `binList`
##' @details This function makes use of R's lexical scoping to
##' evaluate `criteria` (a string), within each bin of `binList`.
##' @param binList a list of bins, each a list which can be cast as an
##' environment for evaluation
##' @param criteria string of logical expressions separated by pipes
##' to be evaluated within each bin of `binList`
##' @return A logical vector of the same length as `binList`.
##' @examples
##' @author Chris Salahub
stopper <- function(binList, criteria) {
    sapply(binList,
           function(b) eval(parse(text = criteria), envir = b))
}
