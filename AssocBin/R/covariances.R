#
# Covariance matrix for a set of blocks of rank pairs
# 
#
# Covariance matrix for block counts
# 
covMat_blocks <- function(N, blocks){
    nblocks <- length(blocks)
    result <- matrix(0, nrow = nblocks, ncol = nblocks)
    # ii <- 1
    for (ii in 1:nblocks){
        # jj <- ii
        block1 <- blocks[[ii]]
        for (jj in ii:nblocks){
            block2 <- blocks[[jj]]
            result[ii, jj] <- blockCov(N, 
                                       block1$rows, block1$cols, 
                                       block2$rows, block2$cols)
            result[jj, ii] <- result[ii,jj]
            
            # jj <- jj+1
        }
        # ii <- ii + 1
    }
    result
}

#
# Actual calculations as determined by paper Appendix
# 
# covariance between counts of two blocks (Note could even over lap, or be identical)
# 
# Handle larger N via logs.
blockCov <- function(N, 
                     block_1_row_indices, block_1_col_indices,
                     block_2_row_indices, block_2_col_indices,
                     counts =  FALSE){
    
    if ((min(c(block_1_row_indices, block_1_col_indices,
               block_2_row_indices, block_2_col_indices)) < 1) |
        (max(c(block_1_row_indices, block_1_col_indices,
               block_2_row_indices, block_2_col_indices)) > N)) stop("Some indices out of range")
    r_A <- length(block_1_row_indices)
    c_A <- length(block_1_col_indices)
    r_B <- length(block_2_row_indices)
    c_B <- length(block_2_col_indices)
    r_I <- length(intersect(block_1_row_indices, block_2_row_indices))
    c_I <- length(intersect(block_1_col_indices, block_2_col_indices))
    
    blockCov_calc_logs(N, r_A, c_A, r_B, c_B, r_I, c_I)
}


blockCov_calc_general <- function(N, r_A, c_A, r_B, c_B, r_I, c_I){
    denominator <- (N^2 * (N - 1))
    r_A_r_I <- r_A - r_I
    c_A_c_I <- c_A - c_I
    r_B_1 <- r_B - 1
    c_B_1 <- c_B - 1
    
    
    denominator <-  2 * log(N) + log(N-1)
    
    part1 <- log(r_A_r_I) + log(c_A_c_I) + log(r_B) + log(c_B)
    part1 <- exp(part1 - denominator)
    
    part2 <- log(r_A_r_I) + log(c_I) + log(r_B) + log(c_B_1)
    part2 <- exp(part2 - denominator)
    
    part3 <- log(r_I) + log(c_A_c_I) + log(r_B_1) + log(c_B)
    part3 <- exp(part3 - denominator)
    
    part4 <- log(r_I) + log(c_I) + log(r_B_1) + log(c_B_1)
    part4 <- exp(part4 - denominator)
    
    partA <- part1 + part2 + part3 + part4
    
    denominator <- 2 * log(N)
    
    part1 <- log(r_A_r_I) + log(c_I) + log(r_B)
    part1 <- exp(part1 - denominator)
    
    part2 <- log(r_I) + log(c_A_c_I) + log(c_B)
    part2 <- exp(part2 - denominator)
    
    part3 <- log(r_I) + log(c_I) + log(c_B_1)
    part3 <- exp(part3 - denominator)
    
    part4 <- log(r_I) + log(c_I) + log(r_B_1)
    part4 <- exp(part4 - denominator)
    
    partB <- part1 + part2 + part3 + part4
    
    partC <- exp(log(N-1) + log(r_I) + log(c_I) - denominator)
    
    result <- partA - partB + partC
    
    
    
    result
}

blockCov_calc_logs <- function(N, r_A, c_A, r_B, c_B, r_I, c_I){
    if((r_I == 0) & (c_I == 0)){
        log_numerator <- sum(log(c(r_A, c_A, r_B, c_B))) # r_A * c_A * r_B * c_B)
        log_denominator <- 2 * log(N) + log(N-1)  # (N^2 * (N - 1))
        result <- exp(log_numerator - log_denominator)
    } else {
        if (r_I == 0) {
            log_part1 <- log(r_A) + log(c_A) + log(r_B) - log(N) - log(N-1)
            part2 <- (c_B/N) - (c_I/c_A)
            result <- exp(log_part1) * part2
        } else {
            if (c_I == 0) {
                log_part1 <- log(r_A) + log(c_A) + log(c_B) - log(N) - log(N-1)
                part2 <- (r_B/N) - (r_I/r_A)
                result <- exp(log_part1) * part2
            }  else {
                
                result <- blockCov_calc_general(N, r_A, c_A, r_B, c_B, r_I, c_I)
            }
        }
    }
    result
}  

# Same calculation without logarithms ... probably not needed.
# 
blockCov_calc_no_logs <- function(N, r_A, c_A, r_B, c_B, r_I, c_I){
    denominator <- (N^2 * (N - 1))
    if((r_I == 0) & (c_I == 0)){
        numerator <- (r_A * c_A * r_B * c_B)
        result <- numerator/denominator
    } else {
        if (r_I == 0) {
            numerator <- (r_A * r_B) * (c_A * c_B - N * c_I)
            result <- numerator/denominator
        } else {
            if (c_I == 0) {
                numerator <- (c_A * c_B) * (r_A * r_B - N * r_I)
                result <- numerator/denominator
            }  else {
                # Now always uses logs.
                result <- blockCov_calc_general(N, 
                                                r_A, c_A,
                                                r_B, c_B, 
                                                r_I, c_I)
            }
        }
    }
    result
}
  
#
# N by N Variance covariance of the Z_ijs , indicators on the rank pair matrix.
# 
VarCovZ <- function(N){
    V <- matrix(0, nrow = N^2, ncol = N^2)
    for(j1 in 1:N){
        for (i1 in 1:N){
            index1 <- i1 + (j1 - 1) * N
            for(j2 in 1:N){
                for (i2 in 1:N){
                    index2 <- i2 + (j2 - 1) * N
                    if((i1 == i2) & (j1 == j2)){
                        V[index1, index2] <- (N -1)/N^2
                    } else{
                        if((i1 == i2) | (j1 == j2)) {
                            V[index1, index2] <- -1/N^2
                        } else {
                            V[index1, index2] <- 1/(N^2 * (N-1))
                        }
                    }
                }
            }
        }
    }
    V
}

#
#  b by b variance-covariance matrix for counts within a set of b blocks
#  on the N by N lattice 
#  - need not be all of the blocks,
#  - need not be a tesselation, could be overlapping (or even nested!) blocks.
#
# Covariance matrix for counts in blocks
covMat_blocks <- function(N, blocks){
    nblocks <- length(blocks)
    result <- matrix(0, nrow = nblocks, ncol = nblocks)
    # ii <- 1
    for (ii in 1:nblocks){
        # jj <- ii
        block1 <- blocks[[ii]]
        for (jj in ii:nblocks){
            block2 <- blocks[[jj]]
            result[ii, jj] <- blockCov(N, block1$rows, block1$cols, block2$rows, block2$cols)
            result[jj, ii] <- result[ii,jj]
            
            # jj <- jj+1
        }
        # ii <- ii + 1
    }
    result
}
