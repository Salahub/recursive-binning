# Measuring statistical dependence using recursive binary splits

This repository contains all of the code required to run a recursive
binning algorithm to measure dependence.

## marbR

The marbR directory contains an implementation of the algorithm to
**m**easure **a**ssociation by **r**ecursive **b**inning in **R**. The
package includes two demos and one real data set summarizing S&P 500
return data. Additionally, data providing null realizations of the
splitting algorithm under random, chi-maximizing, and mutual
information-maximizing splits is provided in three data sets for
convenience.

## R

This directory contains some extraneous data and a script used to
undertake the explorations in my thesis. The script, **experiments.R**,
contains sections which mirror the demos in marbR, as these
experiments served as excellent motivation for demo content.
