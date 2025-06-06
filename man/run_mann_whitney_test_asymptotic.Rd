% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mann_whitney_test.R
\name{run_mann_whitney_test_asymptotic}
\alias{run_mann_whitney_test_asymptotic}
\title{Run Mann-Whitney test (asymptotic)}
\usage{
run_mann_whitney_test_asymptotic(
  Y_list,
  x,
  Z,
  side = "two_tailed",
  implementation = "custom",
  alpha = 0.1
)
}
\arguments{
\item{Y_list}{a list of response vectors}

\item{x}{the treatment vector}

\item{Z}{the covariate matrix}

\item{side}{the side of the test, one of "left", "right", or "two_tailed"}

\item{implementation}{which implementation to use for the MW test, either "custom" (for the in-house implementation) or "r" (for the standard R implementation)}

\item{alpha}{nominal FDR of the BH procedure}
}
\value{
a data frame containing the results
}
\description{
Run Mann-Whitney test (asymptotic)
}
