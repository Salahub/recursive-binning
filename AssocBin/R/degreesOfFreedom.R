##` Degrees of freedom
##' @title Computing a binning's degrees of freedom
##' @description Functions which compute the degrees of freedom of a
##' binning for a chi-squared approximation or parameters for
##' a gamma approximation based on empirical results.
##' @details These exported functions are used to compute parameters
##' needed to approximate the distribution of the chi-squared
##' statistic computed over bins. A full discussion can be found in
##' the accompanying paper.
##' @param nbins the number of bins resulting from recursive random
##' binning
##' @param ncat if one variable is categorical, the number of
##' values the variable can take
##' @return A numeric estimate of the paramter. In the case of
##' degrees of freedom, this is generally not an integer.
##' @author Chris Salahub
##' @describeIn degreesoffreedom Dual continuous fitted df
numNumFittedDf <- function(nbins){
    #(1.0015633 * sqrt(nbins)  -0.8576945)^2
    (sqrt(nbins) - 0.85769)^2
}
##' @describeIn degreesoffreedom Dual continuous simple df
numNumSimpleDf <- function(nbins){
    (sqrt(nbins) - 1)^2
}
##' @describeIn degreesoffreedom Dual continuous gamma shape
numNumGammaShape <- function(nbins){
    (0.1199774 + 0.7214124 * sqrt(ctsSimpleDf(nbins)))^2
}
##' @describeIn degreesoffreedom Dual continuous gamma scale
numNumGammaScale <- function(nbins){
    log_df <- log(ctsSimpleDf(nbins))
    exp(0.4329157 + (1 - 0.9571741) * log_df)
}
##' @describeIn degreesoffreedom Mixed type simple df
facNumSimpleDf <- function(nbins, ncat){
    (ncat -1) * (nbins/ncat -1)
}
##' @describeIn degreesoffreedom Mixed type fitted df
facNumFittedDf <- function(nbins, ncat){
    0.201221 + 0.992706 * catSimpleDf(nbins, ncat)
}
##' @describeIn degreesoffreedom Mixed type gamma shape
facNumGammaShape <- function(nbins, ncat){
    df_simple <- catSimpleDf(nbins, ncat)
    1.102814 * ctsGammaShape(df_simple)
}
##' @describeIn degreesoffreedom Mixed type gamma scale
facNumGammaScale <- function(nbins, ncat){
    df_simple <- catSimpleDf(nbins, ncat)
    log_df <- log(df_simple)
    exp(0.3742961 + (1 - 0.9674642) * log_df)
}
