\name{autoplot.regime_number}
\alias{autoplot.regime_number}

\title{regime_number object ggplot for the outputs on the function outputs mtarnumreg}

\description{
Produces a ggplot object for the results of the mtarnumreg function.
}
\usage{
\method{autoplot}{regime_number}(object, type = 1,...)
}

\arguments{
\item{object}{Object of class \dQuote{\code{regim_number}}. Not NULL}
\item{type}{numeric giving the type of plot to be computed. Allowed
values are 1 for list of \dQuote{\code{Threshold value chains}} for each number of regimes consider (the default) or 2 for "\code{Threshold process}" with the estimated threshold values for each number of regimes consider.
}
\item{...}{other arguments passed to specific methods}
}
\details{
Graph the strings for the outputs corresponding to the functions \dQuote{\code{mtarnumreg}} which return an object of class \dQuote{\code{regim_numreg}}. The chains corresponding to the samplings in each case do not contain the burning period.
}
\value{
  Return a ggplot object.
}
\references{
Calderon, S. and Nieto, F. (2017) \emph{Bayesian analysis of multivariate threshold autoregress models with missing data.} Communications in Statistics - Theory and Methods 46 (1):296–318. doi:10.1080/03610926.2014.990758.
}
\seealso{\code{\link[BMTAR]{mtarnumreg}}}
\author{Valeria Bejarano \email{vbejaranos@unal.edu.co} & Andrey Rincon \email{adrincont@unal.edu.co}}

\examples{
library(ggplot2)
data('datasim_numreg')
autoplot.regime_numreg(datasim_numreg,1)
}

