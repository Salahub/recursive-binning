# Measuring statistical dependence using recursive binary splits

This repository contains all of the code required to run a recursive
binning algorithm to measure dependence.

## AssocBin

The AssocBin directory contains an implementation of the algorithm to
measure association by recursive binning in R. The package includes
two demos and one real data set summarizing S&P 500 return
data. Additionally, data providing null realizations of the splitting
algorithm under random, chi-maximizing, and mutual
information-maximizing splits is provided in three data sets for 
convenience. The package can be downloaded via the command
`install_github("Salahub/recursive-binning", subdir="AssocBin")` using
`install_github` from the `devtools` package in R.

## Scripts

This directory contains some extraneous data and a script used to
undertake the explorations in my thesis. Two scripts are present along
with some examples of null data. **experiments.R** contains sections
which mirror the demos in AssocBin and the thesis, as these
experiments were excellent motivation for demo
content. **sp500process.R** has code modified from Hofert and Oldford 
(2018) to extract and process S&P 500 data.
