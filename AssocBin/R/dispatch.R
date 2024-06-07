##' @title Pairwise binning of complex data
##' @description Measure association for every pair in a data set
##' by providing functions for each combination of data types.
##' @details This wrapper coordinates the application of recursive
##' binning to complex data containing potentially many types.
##' @param data `data.frame` or object coercible to a `data.frame`
##' @param stopCriteria output of `makeCriteria` providing criteria
##' used to stop binning to be passed to binning functions
##' @param dualCategorical function which bins pairs of categorical
##' variables `x` and `y`
##' @param categoricalContinuous function which bins mixed variable
##' combinations and always assumes `x` is the categorical variable
##' @param dualContinuous function which bins pairs of continuous
##' variables `x` and `y`
##' @param measure function which measures association on the binnings
##' of all pairs
##' @param returnBinnings logical; should binnings be returned?
##' @return A list of lists each with elements `x`, `y`, `bnds`,
##' `expn`, `n`, and `stopped`.
##' @author Chris Salahub
pairwiseAssociation <- function(data, stopCriteria,
                                dualCategorical,
                                categoricalContinuous,
                                dualContinuous,
                                measure = binChi,
                                returnBinnings = FALSE) {
    ## argument checking
    if (!is.data.frame(data)) stop("`data` must be a data frame")

    ## necessary function definitions
    stopFn <- function(bns) stopper(bns, stopCriteria)
    if (missing(dualCategorical)) {
        dualCategorical <- catBinner
    }
    if (missing(categoricalContinuous)) {
        categoricalContinuous <-
            function(x, y) {
                uniBinner(x, y, stopper = stopFn,
                          splitter = uniRIntSplit)
            }
    }
    if (missing(dualContinuous)) {
        dualContinuous <-
            function(x, y) {
                binner(x, y, stopper = stopFn,
                       splitter = rIntSplit)
            }
    }
    
    ## pre-processing
    vars <- names(data) # get all variable names
    types <- sapply(data, class) # get all classes
    chars <- types == "character"
    ints <- types == "integer"
    logs <- types == "logical"
    if (any(ints)) { # regularize type names
        data[ints] <- lapply(data[ints], as.numeric)
        types[ints] <- "numeric"
    }
    if (any(chars)) {
        data[chars] <- lapply(data[chars], as.factor)
        types[chars] <- "factor"
    }
    if (any(logs)) {
        data[logs] <- lapply(data[logs],
                             function(x) as.factor(as.numeric(x)))
        types[logs] <- "factor"
    }
    combs <- combn(ncol(data), 2) # get all pairs
    scndRwInds <- types[combs[2, ]] == "factor"
    scndRwFs <- combs[2, scndRwInds]
    combs[2, scndRwInds] <- combs[1, scndRwInds]
    combs[1, scndRwInds] <- scndRwFs # factors always come first
    typecomb <- apply(combs, 2,
                      function(x) paste(types[x], collapse = ":"))
    nlev <- sapply(data, function(var) length(levels(var)))

    ## get ranks
    ranks <- data
    #ranks[, types=="factor"] <- apply(data[, types=="factor"], 2, rank,
    #                                  ties.method = "average")
    ranks[, types=="numeric"] <- apply(data[, types=="numeric"], 2, rank,
                                       ties.method = "random")

    ## bin everything
    binnings <- vector(mode = "list", length(typecomb))
    names(binnings) <- apply(combs, 2,
                             function(x) paste(vars[x],
                                               collapse = ":"))
    facFac <- which(typecomb == "factor:factor")
    facCon <- which(typecomb == "factor:numeric")
    conCon <- which(typecomb == "numeric:numeric")
    for (ii in facFac) { # categorical-categorical combos
        binnings[[ii]] <- dualCategorical(x = ranks[, combs[1, ii]],
                                          y = ranks[, combs[2, ii]])
    }
    for (ii in facCon) { # categorical-numeric combos
        binnings[[ii]] <- categoricalContinuous(x = ranks[, combs[1, ii]],
                                                y = ranks[, combs[2, ii]])
    }
    for (ii in conCon) { # numeric-numeric combos
        binnings[[ii]] <- dualContinuous(x = ranks[, combs[1, ii]],
                                         y = ranks[, combs[2, ii]])
    }

    ## compute objects to return
    binStats <- lapply(binnings, measure)
    K <- sapply(binStats, function(x) length(x$residuals))
    adj <- apply(matrix(nlev[combs], nrow = 2), 2, # constraints
                 function(cl) {
                     max(cl[1] + cl[2] - 1, cl[1], cl[2], 1)
                 })
    dfs <- K - adj
    if (returnBinnings) {
        return(list(df = dfs, statistic = binStats,
                    binning = binnings))
    } else {
        return(list(df = dfs, statistic = binStats))
    }
}
