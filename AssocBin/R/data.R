##' Null replicates of different binning methods
##'
##' These data sets are statistic values for maximized chi, maximized
##' mutual information, and random binning applied to many replicates
##' of independent uniform pairs for a range of different maximum
##' depths. They can be used to determine empirical p-values of the
##' maximized methods. The number following 'null' corresponds to the
##' sample size of the independent uniform pairs, so 'null100' is the
##' data set for samples of 100 independent uniform pairs.
##'
##' @format
##' A list with four elements, the first gives a vector of the maximum
##' depths used, the follow three are arrays corresponding to
##' maximized chi, maximized mutual information, and random splitting.
##' The first dimension of each matrix gives a different statistic on
##' the final bins: the chi squared statistic, the mutual
##' information, the number of bins, and the maximum depth reached
##' before stopping. The second dimension corresponds to the maximum
##' depth from the first element of the list. The final dimension
##' provides the value of each statistic for each maximum depth over
##' 10,000 replicates.
##'
##' @usage data(null100)
"null100"

##' Null replicates of different binning methods
##'
##' These data sets are statistic values for maximized chi, maximized
##' mutual information, and random binning applied to many replicates
##' of independent uniform pairs for a range of different maximum
##' depths. They can be used to determine empirical p-values of the
##' maximized methods. The number following 'null' corresponds to the
##' sample size of the independent uniform pairs, so 'null100' is the
##' data set for samples of 100 independent uniform pairs.
##'
##' @format
##' A list with four elements, the first gives a vector of the maximum
##' depths used, the follow three are arrays corresponding to
##' maximized chi, maximized mutual information, and random splitting.
##' The first dimension of each matrix gives a different statistic on
##' the final bins: the chi squared statistic, the mutual
##' information, the number of bins, and the maximum depth reached
##' before stopping. The second dimension corresponds to the maximum
##' depth from the first element of the list. The final dimension
##' provides the value of each statistic for each maximum depth over
##' 10,000 replicates.
##'
##' @usage data(null1000)
"null1000"

##' Null replicates of different binning methods
##'
##' These data sets are statistic values for maximized chi, maximized
##' mutual information, and random binning applied to many replicates
##' of independent uniform pairs for a range of different maximum
##' depths. They can be used to determine empirical p-values of the
##' maximized methods. The number following 'null' corresponds to the
##' sample size of the independent uniform pairs, so 'null100' is the
##' data set for samples of 100 independent uniform pairs.
##'
##' @format
##' A list with four elements, the first gives a vector of the maximum
##' depths used, the follow three are arrays corresponding to
##' maximized chi, maximized mutual information, and random splitting.
##' The first dimension of each matrix gives a different statistic on
##' the final bins: the chi squared statistic, the mutual
##' information, the number of bins, and the maximum depth reached
##' before stopping. The second dimension corresponds to the maximum
##' depth from the first element of the list. The final dimension
##' provides the value of each statistic for each maximum depth over
##' 10,000 replicates.
##'
##' @usage data(null10000)
"null10000"

##' De-Garched S&P 500 returns
##'
##' This data uses code from the 'zenplots' package to process S&P 500
##' consituent stock returns into uniform pseudo-observations for
##' measuring association.
##'
##' @format
##' A matrix with 755 rows and 461 columns, the rows correspond to
##' dates between 2007 and 2009 and the columns correspond to the
##' different S&P 500 constituent stocks.
##'
##' @usage data(sp500pseudo)
"sp500pseudo"
