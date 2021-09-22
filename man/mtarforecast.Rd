\name{mtarforecast}
\alias{mtarforecast}

\title{Forecast for MTAR model}
\description{Bayesian method to compute forecasts of a MTAR model}
\usage{
mtarforecast(regimemodel,h,level = 0.95, chain = FALSE, modelU)
}
\arguments{
  \item{regimemodel}{object of class \dQuote{\code{regime_model}} of the process
  to forecast. Not NULL}
  \item{h}{numeric type, steps for forecast. Not NULL}
  \item{level}{numeric type, confident interval for estimations. Default 0.95}
  \item{chain}{logical type, if return chains of forecasts. Default TRUE}
  \item{modelU}{object of class \dQuote{\code{regime_model}} of the covariate
  and threshold processes. Not NULL}
}
\details{
Based on the equation of the Multivariate Threshold Autoregressive(MTAR) Model

 \eqn{Y_t= \phi^{(j)}_{0}+ \sum_{i=1}^{p_j} \phi_{i}^{(j)}Y_{t-i}+  \sum_{i=1}^{q_j}\beta_{i}^{(j)}X_{t-i} + \sum{i=1}^{d_j}\delta_{i}^{(j)}Z_{t-i} +\Sigma_{(j)}^{1/2} \epsilon_{t}  if  r_{j-1}< Z_t \le r_{j},}

where process \eqn{\{\epsilon_{t}\}} is a k-variate independent Gaussian process, \eqn{\{Y_t\}} is k-variate process, \eqn{\{X_t\} is a \nu - variate process.}
The function implements Bayesian estimation that gives us the forecasts
for the response and exogenous variable vectors. The coverage percentages
of the forecast intervals and the variability of the predictive distributions.
}
\value{
  Return a list type object
    \item{forecast}{list type object with confident interval, mean value of the forecasts and if the chains if TRUE}
    \item{tsregime}{object of class tsregime of the final Output process with forecast}
    \item{FNDP}{Frobenius norm of the covariance matrix of
the predictive distribution for Yt and Ut}
}
\references{
Calderon, S., & Nieto, F. H. (2021). \emph{Forecasting with Multivariate Threshold Autoregressive Models.} Revista Colombiana de Estadística, 44(2), 369-383.
}
\author{Valeria Bejarano \email{vbejaranos@unal.edu.co} & Andrey Rincon \email{adrincont@unal.edu.co}}
\seealso{
\code{\link[BMTAR]{mtarns}}
\code{\link[BMTAR]{mtarinipars}}
\code{\link[BMTAR]{mtarstr}}
}
\examples{
  data = datasim
  parameters = list(l = 2,
                    orders = list(pj = c(1,1)),
                    r = data$Sim$r)
  initial = mtarinipars(tsregime_obj = data$Sim,
                        list_model = list(pars = parameters))

  estimyt = mtarns(ini_obj = initial,niter = 1000,chain = TRUE)

  Ut = tsregime(datasim$Sim$Zt)
  parameters = list(l = 1,
                    orders = list(pj = 1))
  initial = mtarinipars(tsregime_obj = Ut,
                        list_model = list(pars = parameters))
  estimut = mtarns(ini_obj = initial,niter = 1000,chain = TRUE)

  pred1 = mtarforecast(regimemodel = estimyt,h = 100,modelU = estimut)

}
\keyword{Threshold process}
\keyword{Multivariate threshold autoregressive model}
\keyword{Forecast}
