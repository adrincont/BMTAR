\name{diagnostic_mtar}
\alias{diagnostic_mtar}

\title{Residual diagnosis for model MTAR}
\description{
Tests to help evaluate some assumptions about the MTAR model. calculating some tests and graphs.
}
\usage{
diagnostic_mtar(regime_model,lagmax = NULL,CusumSQ = NULL)
}
\arguments{
 \item{regime_model}{Object of class \dQuote{\code{regime_model}}. Not NULL}
 \item{lagmax}{maximum lag at which to calculate the acf and pacf. Default NULL}
 \item{CusumSQ}{Critical Values of the One-Side Cusumsq Test.}
}

\details{
  For the graphical tests it returns: \dQuote{\code{Residuals plot}} and \dQuote{\code{Residuals density plot}} (overlaps a standard normal density),\dQuote{\code{Residuals plot}} and \dQuote{\code{Residuals plot}}, \dQuote{\code{CUSUM}} statistic for residuals, \dQuote{\code{ACF}} and \dQuote{\code{PACF}} plots for residuals series. \dQuote{\code{CUSUMQ}} statistic for residuals is not returned if the critical value CusumSQ is not given (see Edgerton 1994 https://www.researchgate.net/publication/4994469_Critical_Values
  _for_the_CUSUMSQ_Statistic_in_Medium_and_Large_Sized_Samples).
}
\value{
Returns a list of ggplot objects with the graphics mentioned before.
}
\references{
  Calderon V, S. A., & Nieto, F. H. (2017).  \emph{Bayesian analysis of multivariate       threshold autoregressive models with missing data}. Communications in Statistics-Theory   and Methods, 46(1), 296-318.

  Edgerton, D., & Wells, C. (1994). Critical values for the CUSUMSQ statistic in medium and large sized samples. Oxford Bulletin of Economics and Statistics, 56(3), 355-365.
}
\author{Valeria Bejarano \email{vbejaranos@unal.edu.co}, Sergio Calderón \email{sacalderonv@unal.edu.co} & Andrey Rincon \email{adrincont@unal.edu.co}}
\examples{
library(ggplot2)
data("datasim")
data = datasim$Sim$Zt
parameters = list(l = 1,orders = list(pj = 1))
initial = mtarinipars(tsregime_obj = tsregime(data),
                      list_model = list(pars = parameters))
estim1 = mtarns(ini_obj = initial,niter = 500,chain = TRUE,burn = 500)
diagnostic_mtar(estim1)
}

