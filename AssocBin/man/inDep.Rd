% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dispatch.R
\name{inDep}
\alias{inDep}
\title{Test pairwise variable independence}
\usage{
inDep(
  data,
  stopCriteria,
  catCon = uniRIntSplit,
  conCon = rIntSplit,
  dropPoints = FALSE
)
}
\arguments{
\item{data}{`data.frame` or object coercible to a `data.frame`}

\item{stopCriteria}{output of `makeCriteria` providing criteria
used to stop binning to be passed to binning functions}

\item{catCon}{splitting function to apply to pairs of one
cateogorical and one continuous variable}

\item{conCon}{splitting function to apply to pairs of continuous
variables}

\item{dropPoints}{logical; should returned bins contain points?}
}
\value{
An `inDep` object, with slots `data`, `types`, `pairs`,
`binnings`, `residuals`, `statistics`, `dfs`, `logps`, and
`pvalues` that stores the results of using recursive binning with
the specified splitting logic to test independence on a data set.
`data` gives the name of the data object in the global environment
which was split, `types` is a character vector giving the data
types of each pair, `pairs` is a character vector of the variable
names of each pair, `binnings` is a list of lists where each list
is the binning fir to the corresponding pair by the recursive
binning algorithm, `residuals` is list of numeric vectors giving
the residual for each bin of each pairwise binning, `statistics`
is a numeric vector giving the chi-squared statistic for each
binning, `dfs` is a numeric vector giving the degrees of freedom
of each binning based on the variable type combination and the
final number of bins, `logps` gives the natural logarithm of
the statistic's p-value, and finally `pvalues` is a numeric
vector of p-values for `statistics` assuming a chi-squared null
distribution with `dfs` degrees of freedom.  Internally, the
p-values are computed on the log scale to better distinguish
between strongly dependent pairs and the `pvalues` returned are
computed by calling `exp(logps)`. The order of all returned values
is by increasing `logps`.
}
\description{
This is a high-level function which accepts a data
set, stop criteria, and split functions for continuous variables
and then applies a chi-square test for independence to bins
generated by recursively binning the ranks of continuous variables
or implied by the combinations of levels of categorical variables.
}
\details{
The output of `inDep` is a list, the first element of
which is a list of lists, each of which records the details of the
binning of a particular pair of variables
}
\author{
Chris Salahub
}
