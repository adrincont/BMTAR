\name{autoplot.regime_forecast}
\alias{autoplot.regime_forecast}

\title{regime_forecast object ggplot for the outputs on the function outputs mtarforecast}

\description{
Produces a ggplot object for the results of the mtarforecast function.
}
\usage{
\method{autoplot}{regime_forecast}(object, type = 1,...)
}

\arguments{
\item{object}{Object of class \dQuote{\code{regim_forecast}}. Not NULL}
\item{type}{numeric giving the type of plot to be computed. Allowed
values are 1 for list of \dQuote{\code{Output, Threshold and Covariate process forecast chains}(the default) or 2 for a list "\code{Output,Threshold and Covariate process}" with the estimated forecast and its credibilit bands.
}
\item{...}{other arguments passed to specific methods}
}
\details{
Graph the strings for the outputs corresponding to the functions \dQuote{\code{mtarforecast}} which return an object of class \dQuote{\code{regim_forecast}}. The chains corresponding to the samplings in each case do not contain the burning period.
}
\value{
  Return a ggplot object.
}
\references{
Calderon, S. and Nieto, F. (2017) \emph{Bayesian analysis of multivariate threshold autoregress models with missing data.} Communications in Statistics - Theory and Methods 46 (1):296–318. doi:10.1080/03610926.2014.990758.
}
\seealso{\code{\link[BMTAR]{mtarforecast}}}
\author{Valeria Bejarano \email{vbejaranos@unal.edu.co} & Andrey Rincon \email{adrincont@unal.edu.co}}

\examples{
library(ggplot2)
data('datasim')
estim = mtarns(mtarinipars(tsregime(datasim$Sim$Yt,datasim$Sim$Zt),
list_model = list(pars = list(l = datasim$Sim$l,r = datasim$Sim$r,
orders = list(pj = datasim$pj)))))
estimU = mtarns(mtarinipars(tsregime(datasim$Sim$Zt),
list_model = list(pars = list(l = 1,orders = list(pj = 1)))))
forecast = mtarforecast(estim,h = 10,modelU = estimU)
autoplot.regime_forecast(forecast,1)
}


