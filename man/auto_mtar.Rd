\name{auto_mtar}
\alias{auto_mtar}

\title{Estimation of a MTAR model for some data}
\description{
Compute by Bayesian methodoly a MTAR model for some data
}
\usage{
auto_mtar(Yt, Zt = NULL, Xt = NULL, l0 = 3, maxorders = list(pj = 2,qj = 0,dj = 0),niter = 5000, chain = FALSE, method = 'KUO')
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
}

\details{
  Detalles de los metodos usados en la funcion corespondiente
}
\value{
  Return list type object

  $tsregim class ``tsregime'' object, if missing values completed with estimations
  $numreg numeric type, number of regimes estimated
  $pars class ``regime_model'' object with final estimations of parameters
}
\references{
  Fox, J. (2016)
  \emph{Applied Regression Analysis and Generalized Linear Models},
  Third Edition. Sage.
}
\author{Valeria Bejarano \email{vbejaranos@unal.edu.co}}
\author{Andrey Rincon \email{adrincont@unal.edu.co}}

\examples{
  poner ejemplo de uso de funcion con datos reales
}

\keyword{Bayesian estimation}
