\name{mtarmissing}
\alias{mtarmissing}

\title{Estimation of missing values of observed, covariate and threshold processes}
\description{
Estimation using Bayesian methodology of missing values in observed, covariate and threshold processes
}
\usage{
mtarmissing(ini_obj,niter = 1000, chain = FALSE, level = 0.95,
burn = NULL, cU = 0.5, b = NULL)
}

\arguments{
  \item{ini_obj}{class ``regime_inipars'' object, here specificate in pars: l, orders and r known. Not NULL}
  \item{niter}{numeric type, number of runs of MCMC. Default 1000}
  \item{chain}{logical type, if return chains of parameters. Default FALSE}
  \item{level}{numeric type, confident interval for estimations. Default 0.95}
  \item{burn}{numeric type, number of initial runs. Default NULL (10\% of niter)}
  \item{cU}{numeric type, coeficient of the diagonal covariance matrix of process Ut = (Zt,Xt). Default 0.5}
  \item{b}{numeric type greater or equal 1, autoregressive order of Ut = (Zt,Xt). Default NULL meaning 1}
}

\details{
  Detalles de los metodos usados en la funcion corespondiente
}
\value{
  Return list type object of class ``regime_missing''
  \item{tsregim}{ini_obj$tsregim_obj with estimated observations}
  \item{estimates}{confident interval and mean of estimated missing values}
  \item{Chain}{if chain TRUE, chains of the estimated missing values}
}
\references{
  Calderon V, S. A., & Nieto, F. H. (2017).  \emph{Bayesian analysis of multivariate       threshold autoregressive models with missing data}. Communications in Statistics-Theory   and Methods, 46(1), 296-318.
}
\author{Valeria Bejarano \email{vbejaranos@unal.edu.co} & Andrey Rincon \email{adrincont@unal.edu.co}}

\examples{
data("datasim")
yt = datasim$Sim
# some missing data
data_yt = yt$Yt
data_zt = yt$Zt
posNA = sample(c(1:500),8)
data_yt[c(posNA),] = c(NA,NA)
posNA = sample(c(1:500),8)
data_zt[c(posNA)] = NA
data_final = tsregim(data_yt,data_zt,r = yt$r)
autoplot.tsregim(data_final,1)
autoplot.tsregim(data_final,2)

initial = mtarinipars(tsregim_obj = data_final,
list_model = list(pars = list(l = 2,r = datasim$Sim$r,
orders = list(pj = c(1,1), qj = c(0,0),dj = c(0,0)))))
missingest = mtarmissing(ini_obj = initial,chain = TRUE,
niter = 500,burn = 500)
print(missingest)
autoplot.regime_missing(missingest,1)
datasim$Sim$Yt[is.na(data_yt[,1]),]
missingest$tsregim$Yt[is.na(data_yt[,1]),]
}

\keyword{Bayesian estimation}
\keyword{State Space Form}
