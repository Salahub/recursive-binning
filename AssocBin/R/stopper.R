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
##' makeCriteria(depth >= 5, n < 1)
##' @author Chris Salahub
makeCriteria <- function(...) {
    cl <- match.call() # capturing inputs
    crits <- as.list(cl) # change to a list
    ## remove self reference, collapse into single OR
    paste(c(sapply(crits[-1], deparse), "stopped"), collapse = " | ")
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
##' crits <- makeCriteria(depth >= 5, n < 1)
##' binList1 <- list(makeBin(x = c(1,2), y = c(3,1), depth = 1, n = 2),
##'                 makeBin(x = c(3,4), y = c(2,4), depth = 1, n = 2))
##' binList2 <- list(makeBin(x = c(1,2), y = c(3,1), depth = 6, n = 2),
##'                 makeBin(x = c(), y = c(), depth = 1, n = 0))
##' stopper(binList1, crits)
##' stopper(binList2, crits)
##' @author Chris Salahub
stopper <- function(binList, criteria) {
    sapply(binList,
           function(b) eval(parse(text = criteria), envir = b))
}
