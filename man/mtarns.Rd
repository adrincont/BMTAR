\name{mtarns}
\alias{mtarns}

\title{Estimation of non-structural parameters for MTAR model}
\description{
Bayesian method for estimating non-structural parameters of a MTAR model with prior conjugate.
}
\usage{
mtarns(ini_obj, level = 0.95, burn = NULL, niter = 1000,
chain = FALSE, r_init = NULL)
}

\arguments{
  \item{ini_obj}{class ``regime_inipars'' object, here specificate l and orders known, might know r or Sigma. Not NULL}
  \item{level}{numeric type, confident interval for estimations. Default 0.95}
  \item{burn}{numeric type, number of initial runs. Default NULL (10\% of niter)}
  \item{niter}{numeric type, number of runs of MCMC. Default 1000}
  \item{chain}{logical type, if return chains of parameters. Default FALSE}
  \item{r_init}{numeric type of length l - 1. If r not known, starting value of the chain. Default NULL}
}

\details{
  Based on the equation of the Multivariate Thresold Autoregressive(MTAR) Model

 \eqn{Y_t= \phi^{(j)}_{0}+ \sum_{i=1}^{p_j} \phi_{i}^{(j)}Y_{t-i}+  \sum_{i=1}^{q_j}\beta_{i}^{(j)}X_{t-i} + \sum{i=1}^{d_j}\delta_{i}^{(j)}Z_{t-i} +\Sigma_{(j)}^{1/2} \epsilon_{t}  if  r_{j-1}< Z_t \le r_{j},}

where process \eqn{\{\epsilon_{t}\}} is a k-variate independent Gaussian process, \eqn{\{Y_t\}} is k-variate process, \eqn{\{X_t\} is a \nu - variate process.}
The function implements Bayesian estimation of non-structural parameters of each regime j(\eqn{\phi^{(j)}_{0} \phi_{i}^{(j)}, \beta_{i}^{(j)}, \delta_{i}^{(j)} and \Sigma_{(j)}^{1/2}}) is carried out. The structural parameters: Number of Regimes(l), Thresholds(\eqn{r_1,\cdots,r_{l-1}}), and autoregressive orders(\eqn{p_j,q_j,d_j}) must be known. Prior distributions where selected in order to get conjugate distributions.
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
  \item{r}{final threshold value with acceptance percentage or r if known}
  \item{orders}{list type object with names (pj,qj,dj) known}
}
\references{
  Calderon V, S. A., & Nieto, F. H. (2017).  \emph{Bayesian analysis of multivariate       threshold autoregressive models with missing data}. Communications in Statistics-Theory   and Methods, 46(1), 296-318.
}
\author{Valeria Bejarano \email{vbejaranos@unal.edu.co}, Sergio Calderón \email{sacalderonv@unal.edu.co} & Andrey Rincon \email{adrincont@unal.edu.co}}

\examples{
data("datasim")
data = datasim
#r known
parameters = list(l = 2,
                  orders = list(pj = c(1,1),dj = c(0,0),qj = c(0,0)),
                  r = data$Sim$r)
initial = mtarinipars(tsregim_obj = data$Sim,
                      list_model = list(pars = parameters))
estim1 = mtarns(ini_obj = initial,niter = 1000,chain = TRUE,burn = 1000)
print.regim_model(estim1)
autoplot.regim_model(estim1,2)
autoplot.regim_model(estim1,3)
autoplot.regim_model(estim1,5)
diagnostic.mtar(estim1)
#r unknown
parameters = list(l = 2,orders = list(pj = c(1,1),dj = c(0,0),qj = c(0,0)))
initial = mtarinipars(tsregim_obj = data$Sim,
list_model = list(pars = parameters))
estim2 = mtarns(ini_obj = initial,niter = 500,chain = TRUE)
print.regim_model(estim2)
autoplot.regim_model(estim2,1)
autoplot.regim_model(estim2,2)
autoplot.regim_model(estim2,3)
autoplot.regim_model(estim2,5)
}

\keyword{Bayesian estimation}
\keyword{MCMC}

