\name{mtarstr}
\alias{mtarstr}

\title{Estimation of structural parameters of MTAR model}
\description{
Estimate structural and non-structural parameters of a MTAR model
}
\usage{
mtarstr(ini_obj, level = 0.95, niter = 1000, burn = NULL, chain = FALSE, r_init = NULL)
autoplot(regime_model)
print(regime_model)
}

\arguments{
  \item{ini_obj}{class ``regime_inipars'' object, here specificate in pars: l known, orders not known. Not NULL}
  \item{level}{numeric type, confident interval for estimations. Default 0.95}
  \item{burn}{numeric type, number of initial runs. Default NULL (10\% of niter)}
  \item{niter}{numeric type, number of runs of MCMC. Default 1000}
  \item{chain}{logical type, if return chains of parameters. Default FALSE}
  \item{r_init}{numeric type of length l - 1. If r not known, starting value of the chain. Default NULL}
}

\details{
  Detalles de los metodos usados en la funcion corespondiente
}
\value{
  Return a list type object of class ``regime_model''
  $Nj number of observations in each regime
  $estimates list for each regime with confident interval and mean value of the parameters
  $regime ``regime'' class objects with final estimations
  $Chain  if chain TRUE list type object with parameters chains
  $logLikj log-likelihood of each regime with final estimations
  $data list type object $Yt and $Ut = (Zt,Xt)
  $r final threshold value estimation
  $orders list type object with names (pj,qj,dj) final estimations
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
\keyword{MCMC}
