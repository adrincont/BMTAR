\name{print.regim_missing}
\alias{print.regim_missing}

\title{Print estimates of a regim_missing object of the function output mtarmissing}

\description{
Print estimates output of mtarmissing function.
}
\usage{
\method{print}{regim_missing}(object, ...)
regim_missing(x)
}

\arguments{
\item{object}{Object of class \dQuote{\code{regim_model}}.}
\item{...}{Other plotting parameters to affect the plot.}
}

\details{
Print the estimates for the outputs corresponding to the functions \dQuote{\code{mtarmissing}} which return an object of class \dQuote{\code{regim_missing}}.
}
\value{
  Return to console.
}
\references{
  Calderon V, S. A., & Nieto, F. H. (2017).  \emph{Bayesian analysis of multivariate       threshold autoregressive models with missing data}. Communications in Statistics-Theory   and Methods, 46(1), 296-318.
}
\seealso{\code{\link[MTAR]{mtarmising}}}
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
missingest = mtarmissing(ini_obj = initial,chain = TRUE,niter = 500,burn = 500)
print.regim_missing(missingest)
}
