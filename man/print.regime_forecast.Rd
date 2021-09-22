\name{print.regime_forecast}
\alias{print.regime_forecast}

\title{print regime_forecast object for the function outputs mtarforecast}

\description{
Print estimates for the results of the mtarforecast functions.
}
\usage{
\method{print}{regime_forecast}(object, ...)
}

\arguments{
\item{object}{Object of class \dQuote{\code{regime_forecast}}. Not NULL}
\item{...}{Other print parameters that affect.}
}

\details{
Print estimates outputs corresponding to function \dQuote{\code{mtarforecast}} which return an object of class \dQuote{\code{regime_forecast}}.
}
\value{
  Return to console.
}
\references{
Calderon, S. and Nieto, F. (2017) \emph{Bayesian analysis of multivariate threshold autoregress models with missing data.} Communications in Statistics - Theory and Methods 46 (1):296–318. doi:10.1080/03610926.2014.990758.
}
\seealso{\code{\link[BMTAR]{mtarforecast}}}
\author{Valeria Bejarano \email{vbejaranos@unal.edu.co} & Andrey Rincon \email{adrincont@unal.edu.co}}

\examples{
data("datasim")
data('datasim')
estim = mtarns(mtarinipars(tsregime(datasim$Sim$Yt,datasim$Sim$Zt),
list_model = list(pars = list(l = datasim$Sim$l,r = datasim$Sim$r,
orders = list(pj = datasim$pj)))))
estimU = mtarns(mtarinipars(tsregime(datasim$Sim$Zt),
list_model = list(pars = list(l = 1,orders = list(pj = 1)))))
forecast = mtarforecast(estim,h = 10,modelU = estimU)
print.regime_model(forecast)
}


