#
# Block creation functions
# 
#
# block will be a list of two vectors rows and cols
# containing the row and column indices of the block lattice.
#   - a collection of blocks will just be a list.
#   - full N x N lattice case is assumed even in cases where one
#     or both variates are categorical.
#
# helper functions
# 
get_index_lims <- function(N, nsplit, nkeep, minExpected){  
    # Should return only good splits.
    #   ... works but probably needs a little checking
    #   I think it occasionally returns something that will 
    #   be rejected by having a split that does not meet the
    #   minExpected threshold.  I'm sure it's a simple fix ... rwo
    from <- ceiling(N * minExpected / nkeep)
    to <- nsplit - from
    list(from = from, to = to, 
         new_min_e = from * nkeep/N, success = (to >= from))
}
#
#  Functions for splitting blocks
#  these might not be exported if you the user is to
#  be prevented from using them.
#  ... I'd let the user have access.
#  
#  
#  This is the cts vs cts 
#  
split_block <- function(block, N = NULL, 
                        split_loc = NULL, 
                        minExpected = 5){
    # cts vs cts ... full N by N lattice
    rows <- block$rows
    nrows <- length(rows) 
    cols <- block$cols
    ncols <- length(cols) 
    if (is.null(N)) {  N <- max(nrows, ncols)  }
    if (nrows*ncols/N < 2 * minExpected){
        result <- list(block) 
    } else {
        goodSplit <- FALSE
        
        if (nrows > ncols){
            
            lims <- get_index_lims(N = N, nsplit = nrows, 
                                   nkeep = ncols, 
                                   minExpected = minExpected)
            if (!lims$success) return(list(block))
            possible_splits <- lims$from:lims$to
            
            
            while (!goodSplit) {
                if(is.null(split_loc)) split_loc <- sample(possible_splits, 1)
                rows1 <- rows[1:max(1, split_loc)]
                rows2 <- rows[min(nrows, split_loc+1):nrows]
                e1 <- length(rows1)*ncols /N
                e2 <- length(rows2)*ncols /N
                if((e1 >= minExpected) & (e2 >= minExpected)) { 
                    goodSplit <- TRUE
                } else {
                    split_loc <- NULL 
                }
            }
            result <- list(list(rows = rows1, cols = cols),
                           list(rows = rows2, cols = cols))
        } else {
            if (nrows < ncols){
                lims <- get_index_lims(N = N, nsplit = ncols, nkeep = nrows, minExpected = minExpected)
                if (!lims$success) return(list(block))
                possible_splits <- lims$from:lims$to
                while (!goodSplit) {
                    if(is.null(split_loc)) split_loc <- sample(sample(possible_splits, 1), 1)
                    cols1 <- cols[1:max(1, split_loc)]
                    cols2 <- cols[min(ncols, split_loc+1):ncols]
                    e1 <- nrows * length(cols1) /N
                    e2 <- nrows * length(cols2) /N
                    if((e1 >= minExpected) & (e2 >= minExpected))  { 
                        goodSplit <- TRUE 
                    } else {
                        split_loc <- NULL
                    }
                }
                result <- list(list(rows = rows, cols = cols1),
                               list(rows = rows, cols = cols2))
                
            } else {
                lims <- get_index_lims(N = N, nsplit = ncols, nkeep = nrows, minExpected = minExpected)
                if (!lims$success) return(list(block))
                possible_splits <- lims$from:lims$to
                if (runif(1) <= 0.5){
                    while (!goodSplit) {
                        if(is.null(split_loc)) split_loc <- sample(sample(possible_splits, 1), 1)
                        rows1 <- rows[1:max(1, split_loc)]
                        rows2 <- rows[min(nrows, split_loc+1):nrows]
                        e1 <- length(rows1) * ncols /N
                        e2 <- length(rows2) * ncols /N
                        if((e1 >= minExpected) & (e2 >= minExpected)) { 
                            goodSplit <- TRUE 
                        } else {
                            split_loc <- NULL
                        }
                        
                    }
                    result <- list(list(rows = rows1, cols = cols),
                                   list(rows = rows2, cols = cols))
                    
                } else {
                    while (!goodSplit) {
                        if(is.null(split_loc)) split_loc <- sample(possible_splits, 1)
                        cols1 <- cols[1:max(1, split_loc)]
                        cols2 <- cols[min(ncols, split_loc+1):ncols]
                        e1 <- nrows * length(cols1) /N
                        e2 <- nrows * length(cols2) /N
                        if((e1 >= minExpected) & (e2 >= minExpected))  {  
                            goodSplit <- TRUE
                        } else {
                            split_loc <- NULL
                        }
                        
                    }
                    result <- list(list(rows = rows, cols = cols1),
                                   list(rows = rows, cols = cols2))
                }
            }
        }
    }
    result
}

#
# As the name suggests. 
#   equal_bins = TRUE means all of the same size (up to rounding)
#   ... special case of binning.  Primarily to be used when 
#       we wanted to fake a cat vs cat situation ... then splitting
#       was done with same_bins = TRUE.  Just to divide the lattice
#       in rows across the columns defining each category.
# 
split_cat_into_nbins <- function(cat, nbins, 
                                 equal_bins = TRUE,
                                 minExpected = 0) {
    rows <- cat$rows
    cols <- cat$cols
    nrows <- length(rows)
    bins <- list()
    if (equal_bins){
        nperbin <- floor(nrows/nbins)
        ends <- cumsum(rep(nperbin, nbins))
    } else {
        adjust <- nbins * minExpected
        if (adjust >= nrows){stop("N not large enough")}
        nperbin <- minExpected + 
            rmultinom(1, size = (nrows - adjust), prob = rep(1/nbins, nbins))
        ends <- cumsum(nperbin)
    }
    starts <- ends - (nperbin -1)
    ends[nbins] <- nrows
    for (i in 1:nbins){
        bins[[i]] <- list(rows = rows[starts[i]:ends[i]], cols = cols)
    }
    bins 
}

#
# Splitting a block within a category.
# 
split_cat_block <- function(block, N, minExpected = 0, split_loc = NULL){
    rows <- block$rows
    nrows <- length(rows) 
    cols <- block$cols
    ncols <- length(cols)
    if (nrows*ncols/N < 2 * minExpected){
        result <- list(block)  # Not sure if this is what should be returned in this case.
    } else {
        goodSplit <- FALSE
        lims <- get_index_lims(N = N, 
                               nsplit = nrows,  # category shares columns always; split rows
                               nkeep = ncols, 
                               minExpected = minExpected)
        if (!lims$success) return(list(block))
        # stop(paste("N =", N, "  nrows =", nrows, "  ncols =", ncols,
        #                                   " from =", lims$from, "  to =", lims$to))
        possible_splits <- seq(lims$from, lims$to, 1)
        
        while (!goodSplit) {
            if(is.null(split_loc)) split_loc <- sample(possible_splits, 1)
            rows1 <- rows[1:max(1, split_loc)]
            rows2 <- rows[min(nrows, split_loc+1):nrows]
            e1 <- length(rows1)*ncols /N
            e2 <- length(rows2)*ncols /N
            if((e1 >= minExpected) & (e2 >= minExpected)) { 
                goodSplit <- TRUE
            } else {
                
                print(paste("fail:",
                            "split_loc =", split_loc,  " from =", lims$from, "  to =", lims$to))
                split_loc <- NULL
            }
        }
        result <- list(list(rows = rows1, cols = cols),
                       list(rows = rows2, cols = cols))
    }
    result
}

#
#  Creating blocks by random recursive binning (squarified)
#  ... i.e., from top level
#  
#  - handy if user wants to just reuse the same blocks on differet
#    data sets.

#  Cts versus cts
#  
create_random_blocks <- function(N, depth = 2, minExpected = 0){
    rows <- 1:N
    cols <- 1:N
    blocks <- list(list(rows = rows, cols = cols))
    while(depth > 0){
        blocks <- unlist(Map(function(block) {
            split_block(block, N = N, minExpected = minExpected)}, 
            blocks), 
            recursive = FALSE)
        depth <- depth - 1
    }
    blocks
}

# Case where one variable is a category
# 
# Note that this is generating random categories too.
# That is, it is really to generate simulated categorical
# data AND bin it.
# 
# For the package, we might want to rewrite and have the categorical
# variate counts as input.  A much simpler function to write.
# 
# 
# 
# Lots of ways to specify.
# - ncat and prob  ... generates category number and (probable) size
#                  ... really just for simulation purposes.
#                  
# - equal_bins ... used with nbins, 
#                  same number (subject to  rounding) in each bin  
#                  within every ccategory
# - same_bins ... used with nbins, 
#                 use the same bins limits across categories
#                 can be used to create contingency tables on N x N lattice
# 

create_cat_blocks  <- function(N, 
                               depth = 2, 
                               minExpected = 5, 
                               cat_counts = NULL,
                               nbins = NULL, 
                               equal_bins = FALSE,
                               same_bins = FALSE
                               ){
    if (is.null(cat_counts)) stop("Argument cat_counts (number in each category) is required.")
    if (N != sum(cat_counts)) warning("counts in categories do not sum to N")
    cat_ends <- cumsum(cat_counts)
    cat_starts <- c(1, 1 + cat_ends[-ncat])
    cat_indices <- matrix(c(cat_starts, cat_ends), 
                          nrow = ncat, ncol = 2,
                          byrow = FALSE)
    
    rows <- 1:N
    cols <- 1:N
    cats <- list()
    for (i in 1:ncat){
        cats[[i]] <- list(rows = rows, 
                          cols = cat_indices[i,1]:cat_indices[i,2])
    }
    if (is.numeric(nbins)){
        if (same_bins){
            
            cats[[1]] <- split_cat_into_nbins(cat = cats[[1]], 
                                              nbins = nbins, 
                                              equal_bins = equal_bins,
                                              minExpected = minExpected)
            for (i in 2:ncat){
                cats[[i]] <- Map(function(cat) {
                    list(rows = cat$rows, cols = cats[[i]]$cols)}, 
                    cats[[1]])
            }
            blocks <- unlist(cats, recursive = FALSE)
        } else {
            blocks <- unlist(Map(function(cat) {
                split_cat_into_nbins(cat, 
                                     nbins = nbins, 
                                     equal_bins = equal_bins,
                                     minExpected = minExpected)}, 
                cats), 
                recursive = FALSE)
        }
    } else { # random split
        remaining_depth <- depth
        while(remaining_depth > 0){
            for (i in 1:ncat){
                cat <- cats[[i]]
                cats[[i]] <-  split_cat_block(cat, 
                                              N = N, 
                                              minExpected = minExpected)
            }
            remaining_depth <- remaining_depth - 1
            cats <- unlist(cats, recursive = FALSE)
            ncat <- length(cats)
        }
        blocks <- cats
    }
    
    blocks 
}
    
generate_cat_blocks  <- function(N, depth = 2, minExpected = 5, 
                               ncat = 3, prob = rep(1/ncat, ncat),
                               nbins = NULL, equal_bins = FALSE, 
                               same_bins = FALSE){
   
    if (is.null(nbins)) {
        # make sure there are enough in each category
        adjust <- (2^depth -1) * minExpected 
    } else {
        adjust <- nbins * minExpected
    }
    if (adjust*ncat >= N){stop("N not large enough")}
    
    cat_counts <- adjust + 
        rmultinom(1, size = (N - adjust * ncat), prob = prob)  
    
    blocks <- create_cat_blocks(N, depth = depth, 
                                minExpected = minExpected, 
                                cat_counts = cat_counts,
                                nbins = nbins, 
                                equal_bins = equal_bins,
                                same_bins = same_bins)
    blocks
}

#
#  Number of blocks by category (in columns)
#  
#  
block_count_by_cat <- function(blocks){
    # categories are defined to be groups of identical columns
    # (assumed for all N ... no check)
    cats <- unique(Map(function(block){block$cols}, blocks))
    row_counts_per_category <- 
        Map(function(cat){
            unlist(
                Map(function(block){
                    length(block$rows)},
                    Filter(function(block){identical(block$cols, cat)}, 
                           blocks)
                ), 
                recursive = FALSE)
        }, 
        cats)
    row_counts_per_category     
}

#
# As the name suggests
# 
nbins_per_category <- function(blocks){
    row_counts_per_category <- block_count_by_cat(blocks)
    unlist(Map(length, row_counts_per_category))
}


#
#  Getting expected values
#  

get_expected <- function(N, blocks) {
    unlist(Map(function(block){
        length(block$rows) * length(block$cols) /N
    }, blocks))
}

#
#  Calculating X2 values
#  

get_X2 <- function(x, y, blocks, 
                   expected = NULL, 
                   pseudo_obs = FALSE){
    # x and y are ranks if pseudo_obs = FALSE
    # and numbers between 0 and 1 if pseudo_obs = TRUE
    N <- length(x)
    nBlocks <- length(blocks)
    if(is.null(expected)) { expected <- get_expected(N, blocks)}
    X2 <- 0
    for (i in 1:nBlocks){
        block <- blocks[[i]]
        rangeX <- range(block$cols)
        rangeY <- range(block$rows)
        if (pseudo_obs){
            rangeX <- rangeX/(N+1)
            rangeY <- rangeY/(N+1)
        }
        casesInBlock <- (x <= max(rangeX)) & (x >= min(rangeX)) & 
            (y <= max(rangeY)) & (y >= min(rangeY))
        Oi <- sum(casesInBlock)
        Ei <- expected[i]
        X2 <- X2 + (Oi - Ei)^2/Ei
    }
    X2
}

#
# calculating p-value via prob. integral
# 

pit_pvalues <- function(rank_x, rank_y, blocks, Nreps = 100) {
    p_values <- numeric(length = Nreps)
    N <- length(rank_x)
    b <- length(blocks)
    expected_values <- get_expected(N, blocks)
    for (i in 1:Nreps){
        unif_x <- sort(runif(N))[rank_x]
        unif_y <- sort(runif(N))[rank_y]
        X2 <- get_X2(x = unif_x, y = unif_y, 
                     blocks = blocks,
                     expected = expected_values,
                     pseudo_obs = TRUE)
        p_values[i] <- pchisq(X2, df = (b - 1), lower.tail = FALSE)
    }
    p_values
}
