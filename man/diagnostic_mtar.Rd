\name{diagnostic_mtar}
\alias{diagnostic_mtar}

\title{Residual diagnosis for model MTAR}
\description{
Tests to help evaluate some assumptions about the MTAR model. calculating some tests and graphs.
}
\usage{
diagnostic_mtar(regim_model,lagmax = NULL)
}
\arguments{
 \item{regim_model}{Object of class \dQuote{\code{regim_model}}.}
 \item{lag.max}{maximum lag at which to calculate the acf and pacf.}
}

\details{
  For the graphical tests it returns: \dQuote{\code{Residuals plot}} and \dQuote{\code{Residuals density plot}} (overlaps a standard normal density),\dQuote{\code{Residuals plot}} and \dQuote{\code{Residuals plot}}, \dQuote{\code{CUSUM}} and \dQuote{\code{CUSUMQ}} statistic for residuals (alpha=0.05), \dQuote{\code{ACF}} and \dQuote{\code{PACF}} plots for residuals series.
}
\value{
Returns a list of ggplot objects with the graphics mentioned before.
}
\references{
  Calderon V, S. A., & Nieto, F. H. (2017).  \emph{Bayesian analysis of multivariate       threshold autoregressive models with missing data}. Communications in Statistics-Theory   and Methods, 46(1), 296-318.
}
\author{Valeria Bejarano \email{vbejaranos@unal.edu.co} & Andrey Rincon \email{adrincont@unal.edu.co}}
\examples{
library(ggplot2)
data('data_ns')
diagnostic.mtar(estim1)
}
