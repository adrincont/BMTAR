\name{mtarsim}
\alias{mtarsim}

\title{Multivariate threshold autoregressive process simulation}
\description{
Simulate a MTAR process
}
\usage{
mtarsim(N, Rg, r = NULL, Xt = NULL, Zt = NULL, seed = NULL, ...)
}

\arguments{
  \item{N}{numeric type greater than 0. Number of observation to simulate. Not NULL}
  \item{Rg}{list type object of length l number of regimes of the process with names (R1, ..., Rl), each a class ``regime'' object. Not NULL}
  \item{r}{numeric type of length l - 1, threshold value. Default NULL}
  \item{Xt}{matrix type object, covariate process. Default NULL}
  \item{Zt}{matrix type object, threshold process. Default NULL}
  \item{seed}{numeric type, set a seed for simulation}
}

\details{
  Detalles de los metodos usados en la funcion corespondiente
}
\value{
  Return a list type object of class ``mtarsim''
  $Sim object class ``tsregim''
  $Reg list type object with names (R1, ..., Rl) each one class ``regime''
  $pj vector of autoregressive orders in each regime
  $qj vector of covariate lags orders in each regime
  $dj vector of lags orders of threshold process in each regime
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

\keyword{MTAR}
\keyword{Covariate process}
\keyword{Threshold process}

