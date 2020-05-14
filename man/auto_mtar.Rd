\name{auto_mtar}
\alias{auto_mtar}

\title{Estimation of a MTAR model for some data}
\description{
Compute by Bayesian methodology a MTAR model for some data
}
\usage{
auto_mtar(Yt, Zt = NULL, Xt = NULL, l0 = 3, maxorders = list(pj = 2,qj = 0,dj = 0),niter = 5000, chain = FALSE, method = 'KUO',parallel = FALSE)
}

\arguments{
  \item{Yt}{matrix type object, observed process. Not NULL}
  \item{Zt}{matrix type object, threshold process. Default NULL}
  \item{Xt}{matrix type object, covariate process. Default NULL}
  \item{l0}{numeric type between 2 and 4, number of regimes to consider. Default 3}
  \item{maxorders}{list type object with names (pj,qj,dj), maximum lags consider for the processes in each regime. Default pj = 2, qj = 0,dj = 0}
  \item{niter}{numeric type, number of runs for every estimation. Default 5000}
  \item{chain}{logical type, if return chains of estimations parameters and values (if missing)}
  \item{method}{character type, must be one ``KUO'' or ``SSVS''}
  \item{parallel}{logical, if \code{parallel} package should be used. Default FALSE}
}

\details{
  The default arguments are designed for rapid estimation of models for any data (Yt, Zt and Xt). Returns the fit of MTAR model.
  The function conducts Bayesian estimation with ``niter'' chains of the number of regimes with maximum ``l0'' and within the maximum lags orders provided ``maxorders''.
  It can be a little be faster when used ``parallel'' for \code{parallel} package.
}
\value{
  Return list type object

  \item{tsregim}{class ``tsregime'' object, if missing values completed with estimations}
  \item{numreg}{numeric type, number of regimes estimated}
  \item{pars}{class ``regime_model'' object with final estimations of parameters}
}
\references{
  Calderon V, S. A., & Nieto, F. H. (2017).  \emph{Bayesian analysis of multivariate threshold autoregressive models with missing data}. Communications in Statistics-Theory and Methods, 46(1), 296-318.
}
\author{Valeria Bejarano \email{vbejaranos@unal.edu.co} & Andrey Rincon \email{adrincont@unal.edu.co}}

\examples{
  data('datasim')
  data = datasim$Sim
  auto_mtar(Yt = data$Yt, Zt = data$Zt)
}

\keyword{Bayesian estimation}
