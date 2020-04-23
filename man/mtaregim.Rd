\name{mtaregim}
\alias{mtaregim}

\title{Object class ``regime'' creation}
\description{
Organization of parameters among regimes
}
\usage{
mtaregim(orders = list(p = 1,q = 1,d = 1), cs = NULL, Phi, Beta = NULL, Delta = NULL, Sigma, ...)
}

\arguments{
  \item{orders}{list type object with names (p,q,d), number of lags of Yt, Xt and Zt, respectively. Default p = 1, q = 0, d = 0}
  \item{cs}{matrix type object, the constant term of the regime specification. Default NULL}
  \item{Phi}{list type object with names (phi1, ..., phip), each one a matrix type object autoregressive specification. Not NULL}
  \item{Beta}{list type object with names (beta1, ..., betaq), each one a matrix type object covariate parameters specification Default NULL}
  \item{Delta}{list type object with names (delta1, ..., deltad), each one a matrix type object parameter specification of Threshold process. Default NULL}
  \item{Sigma}{matrix type object, specfication of errors covariate matrix. Not NULL}
}

\details{
  Detalles de los metodos usados en la funcion corespondiente
}
\value{
  Return list type object of class ``regime'' with the values of the arguments
}
\references{
  Fox, J. (2016)
  \emph{Applied Regression Analysis and Generalized Linear Models},
  Third Edition. Sage.
}
\author{Valeria Bejarano \email{vbejaranos@unal.edu.co} & Andrey Rincon \email{adrincont@unal.edu.co}}

\seealso{\code{\link[MTAR]{mtaregim}}, \code{\link[MTAR]{mtarns}},
  \code{\link[MTAR]{mtaregim}}, \code{\link{showLabels}}}

\examples{
  poner ejemplo de uso de funcion con datos reales
}

\keyword{Threshold process}
\keyword{Multivariate threshold autoregressive model}
\keyword{Regime}

