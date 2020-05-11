\name{mtaregim}
\alias{mtaregim}

\title{Object class ``regime'' creation}
\description{
Create an object of class ``regime'' given nonstructural and structural parameters for regime i.
}
\usage{
mtaregim(orders = list(p = 1,q = 1,d = 1), cs = NULL, Phi, Beta = NULL, Delta = NULL, Sigma, ...)
}
\arguments{
  \item{orders}{list type object with names (p,q,d), number of lags of Yt, Xt and Zt, respectively. Default p = 1, q = 0, d = 0}
  \item{cs}{matrix type object, the constant term of the regime specification. Default NULL}
  \item{Phi}{list type object with names (phi1, ..., phip), each one a matrix \eqn{(kxk)} type object autoregressive specification. Not NULL}
  \item{Beta}{list type object with names (beta1, ..., betaq), each one a matrix \eqn{(kx\nu)} type object covariate parameters specification Default NULL}
  \item{Delta}{list type object with names (delta1, ..., deltad), each one a matrix \eqn{(kx1)} type object parameter specification of Threshold process. Default NULL}
  \item{Sigma}{a positive-definite symmetric matrix \eqn{(kxk)}, specfication of errors covariate matrix. Not NULL}
}
\details{
  Causes creation of the object class ``regime''. Sigma matrix corresponds to \eqn{\Sigma} (root of the covariance matrix). When cs is not specified or only matrices are delivered for some lags, the function assumes unspecified 0 (matrix). Rows number in the Phi, Beta and Delta matrix should be the same (k dimension of variables in Yt).
}
\value{
  Return list type object of class ``regime'' with the values of the arguments
}
\references{
  Calderon V, S. A., & Nieto, F. H. (2017).  \emph{Bayesian analysis of multivariate       threshold autoregressive models with missing data}. Communications in Statistics-Theory   and Methods, 46(1), 296-318.
}
\author{Valeria Bejarano \email{vbejaranos@unal.edu.co} & Andrey Rincon \email{adrincont@unal.edu.co}}

\seealso{\code{\link[MTAR]{mtaregim}}, \code{\link[MTAR]{mtarns}},
  \code{\link[MTAR]{mtaregim}}}

\examples{
# Creation of parameters for regimen with orders = c(2,1,1), nu = 1 y k = 2.
## previous objects
orders = list(p = 2,q = 1,d = 1)
Phi = list(phi2 = matrix(c(0.1,0.6,-0.4,0.5),2,2, byrow = T))
Beta = list(beta1 = matrix(c(0.3,-0.4),2, 1))
Delta = list(delta1 = matrix(c(0.6,1),2,1))
Sigma = matrix(c(1,0.6,0.6,1.5),2,2,byrow = T)
cs = matrix(c(1,-1),nrow = 2)
## creacion de la clase regime
Ri = mtaregim(orders = orders,Phi = Phi,Beta = Beta,Delta = Delta,Sigma = Sigma,cs = cs)
class(Ri)
}
\keyword{Threshold process}
\keyword{Multivariate threshold autoregressive model}
\keyword{Regime}

