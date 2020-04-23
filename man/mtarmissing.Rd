\name{mtarmissing}
\alias{mtarmissing}

\title{Estimation of missing values of observed, covariate and threshold processes}
\description{
Estimation using Bayesian methology of missing values in observed, covariate and threshold processes
}
\usage{
mtarmissing(ini_obj,niter = 1000, chain = FALSE, level = 0.95, burn = NULL, cU = 0.5, b = NULL)
autoplot(regime_missing)
print(regime_missing)
}

\arguments{
  \item{ini_obj}{class ``regime_inipars'' object, here specificate in pars: l, orders and r known. Not NULL}
  \item{niter}{numeric type, number of runs of MCMC. Default 1000}
  \item{chain}{logical type, if return chains of parameters. Default FALSE}
  \item{level}{numeric type, confident interval for estimations. Default 0.95}
  \item{burn}{numeric type, number of initial runs. Default NULL (10\% of niter)}
  \item{cU}{numeric type, coeficient of the diagonal covariance matrix of process Ut = (Zt,Xt). Default 0.5}
  \item{b}{numeric type greater or equal 1, autoregressive order of Ut = (Zt,Xt). Default NULL meaning 1}
}

\details{
  Detalles de los metodos usados en la funcion corespondiente
}
\value{
  Return list type object of class ``regime_missing''
  $tsregim = ini_obj$tsregim_obj
  $estimates confident interval and mean of estimated missing values
  $Chain if chain TRUE, chains of the estimated missing values
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
\keyword{State Space Form}
