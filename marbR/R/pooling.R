## some pool functions
poolChi <- function(p, k) {
    M <- length(p) # dimension
    pchisq(sum(qchisq(p, df = k, lower.tail = FALSE)), df = M*k,
           lower.tail = FALSE)
}

poolNorm <- function(p, mu = 0, sd = 1) {
    M <- length(p) # dimension
    pnorm(sum(qnorm(p, mean = mu, sd = sd, lower.tail = FALSE)),
          mean = M*mu, sd = sqrt(M)*sd, lower.tail = FALSE)
}

poolGamma <- function(p, shape = 1, rate = 1/2) {
    M <- length(p) # dimension
    pgamma(sum(qgamma(p, shape = shape, rate = rate,
                      lower.tail = FALSE)),
           shape = M*shape, rate = rate,
           lower.tail = FALSE) # sum shapes
}

poolTip <- function(p) {
    1 - (1 - min(p))^(length(p))
}

## chi-specific central, marginal tendencies
chiPc <- function(k, M, alpha = 0.05) {
    pchisq(qchisq(alpha, df = M*k, lower.tail = FALSE)/M, df = k,
           lower.tail = FALSE)
}

chiPr <- function(k, M, alpha = 0.05) {
    pchisq(qchisq(alpha, df = M*k, lower.tail = FALSE), df = k,
           lower.tail = FALSE)
}

## two ways of computing the centrality quotient:
## conditional distribution
centQuot <- function(k, M, alpha = 0.05) {
    chiCrit <- qchisq(alpha, df = M*k, lower.tail = FALSE)
    1 - pchisq(chiCrit, df = k, lower.tail = FALSE)/
        pchisq(chiCrit/M, df = k, lower.tail = FALSE)
}
## directly:
centQuotRat <- function(k, M, alpha = 0.05) {
    1 - chiPr(k, M, alpha)/chiPc(k, M, alpha)
}

## the inverse: compute kappa for a given centrality quotient
chiKappa <- function(cq, M, alpha = 0.05, interval = c(0,100),
                     tol = .Machine$double.eps^0.5) {
    uniroot(function(x) centQuotRat(x, M, alpha) - cq,
            interval = interval, tol = tol)$root
}

## the heard-rubin-delanchy case
## pooled p value functions
hrBeta <- function(p, w = 1) {
    w*sum(log(p)) - (1 - w)*sum(log(1 - p))
}

## simulated rejection bounds hrb
hrBnd <- function(w, alpha = 0.05, M = 2, nsim = 1e4) {
    dat <- matrix(runif(M*nsim), ncol = M) # simulate data
    vals <- apply(dat, 1, hrBeta, w = w)
    quantile(vals, 1 - alpha) # get simulated quantiles
}

## correctly calibrated pooled hrb p-value based on simulation
poolHR <- function(w = 1, M = 2, nsim = 1e5) { # simulated pooled p
    dat <- matrix(runif(M*nsim), ncol = M)
    pools <- apply(dat, 1, hrBeta, w = w)
    function(p) { # closure to limit repeated computation
        mean(hrBeta(p, w) >= pools) # obs quant
    }
}

## test based on the hrb and simulation
testHR <- function(alpha = 0.05, M = 2, nsim = 1e4, neval = 50) {
    dat <- matrix(runif(M*nsim), ncol = M)
    ws <- seq(0, 1, length.out = neval)
    qs <- sapply(ws, # simulated bounds across multiple simulations
                 function(w.) quantile(apply(dat, 1, hrBeta, w = w.),
                                       1 - alpha))
    function(p, w) { # closure to reduce computation
        pool <- hrBeta(p, w)
        wind <- sum(ws <= w)
        bnd <- (ws[wind]*(ws[wind + 1] - w) + ws[wind+1]*(w - ws[wind]))/
            (ws[wind + 1] - ws[wind])
        pool < bnd
    }
}
