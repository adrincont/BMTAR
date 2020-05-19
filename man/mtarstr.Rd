\name{mtarstr}
\alias{mtarstr}

\title{Estimation of structural parameters of MTAR model}
\description{
Estimate structural and non-structural parameters of a MTAR model when the number of regimes is fixed.
}
\usage{
mtarstr(ini_obj, level = 0.95, niter = 1000, burn = NULL, chain = FALSE,
r_init = NULL, parallel = FALSE)
}

\arguments{
  \item{ini_obj}{class ``regime_inipars'' object, here specificate in pars: l known, orders not known. Not NULL}
  \item{level}{numeric type, confident interval for estimations. Default 0.95}
  \item{burn}{numeric type, number of initial runs. Default NULL (10\% of niter)}
  \item{niter}{numeric type, number of runs of MCMC. Default 1000}
  \item{chain}{logical type, if return chains of parameters. Default FALSE}
  \item{r_init}{numeric type of length l - 1. If r not known, starting value of the chain. Default NULL}
  \item{parallel}{logical type, if package \code{parallel} should be used. Default FALSE}
}

\details{
  If the number of regimes \eqn{l} is known or fixed, we can estimate other structural parameters of the MTAR model: Thresholds(\eqn{r_1,\cdots,r_{l-1}}), and autoregressive orders(\eqn{p_j,q_j,d_j}). Of course, the non-structural parameters are also estimated. The problem of estimation the autoregressive orders is addressed to the problem of Bayesian variable selection in regression using Gibss Variable selection(GVS) or Kuo and Mallick Methodologies. Samples of the full conditional distribution for Threshold values are extracted using Random Walk Metropolis-Hastings Algorithm.
}
\value{
  Return a list type object of class ``regime_model''
  \item{Nj}{number of observations in each regime}
  \item{estimates}{list for each regime with confident interval and mean value of the parameters}
  \item{regime}{``regime'' class objects with final estimations}
  \item{Chain}{if chain TRUE list type object with parameters chains}
  \item{fitted.values}{matrix type object with fitted.values of the estimated model}
  \item{residuals}{matrix type object with residuals of the estimated model}
  \item{logLikj}{log-likelihood of each regime with final estimations}
  \item{data}{list type object $Yt and $Ut = (Zt,Xt)}
  \item{r}{final threshold value estimation with acceptance percentage}
  \item{orders}{list type object with names (pj,qj,dj) final estimations}
}
\references{
  Calderon V, S. A., & Nieto, F. H. (2017).  \emph{Bayesian analysis of multivariate       threshold autoregressive models with missing data}. Communications in Statistics-Theory   and Methods, 46(1), 296-318.
}
\author{Valeria Bejarano \email{vbejaranos@unal.edu.co}, Sergio Calderón \email{sacalderonv@unal.edu.co} & Andrey Rincon \email{adrincont@unal.edu.co}}

\examples{
data("datasim")
data = datasim
# KUO method
initial = mtarinipars(tsregim_obj = data$Sim,method = 'KUO',
list_model = list(pars = list(l = 2),orders = list(pj = c(2,2))))
estruc = mtarstr(ini_obj = initial,niter = 500,chain = TRUE)
autoplot.regim_model(estruc,1)
autoplot.regim_model(estruc,2)
autoplot.regim_model(estruc,3)
autoplot.regim_model(estruc,4)
autoplot.regim_model(estruc,5)
# SSVS method
initial = mtarinipars(tsregim_obj = data$Sim,method = 'SSVS',
list_model = list(pars = list(l = 2),orders = list(pj = c(2,2))))
estruc = mtarstr(ini_obj = initial,niter = 500,chain = TRUE)
autoplot.regim_model(estruc,1)
autoplot.regim_model(estruc,2)
autoplot.regim_model(estruc,3)
autoplot.regim_model(estruc,4)
autoplot.regim_model(estruc,5)
}

\keyword{Bayesian estimation}
\keyword{MCMC}
