\name{print.regime_number}
\alias{print.regime_number}

\title{print regime_number object for the function outputs mtarnumreg}

\description{
Print estimates for the results of the mtarnumreg function.
}
\usage{
\method{print}{regime_number}(object, ...)
}

\arguments{
\item{object}{Object of class \dQuote{\code{regime_number}}. Not NULL}
\item{...}{Other print parameters that affect.}
}

\details{
Print estimates outputs corresponding to function \dQuote{\code{mtarnumreg}} which return an object of class \dQuote{\code{regime_number}}.
}
\value{
  Return to console.
}
\references{
  Calderon V, S. A., & Nieto, F. H. (2017).  \emph{Bayesian analysis of multivariate       threshold autoregressive models with missing data}. Communications in Statistics-Theory   and Methods, 46(1), 296-318.
}
\seealso{\code{\link[MTAR]{mtarnumreg}}}
\author{Valeria Bejarano \email{vbejaranos@unal.edu.co} & Andrey Rincon \email{adrincont@unal.edu.co}}

\examples{
data("datasim")
data = datasim$Sim$Zt
parameters = list(l = 1,orders = list(pj = 1))
initial = mtarinipars(tsregim_obj = tsregime(data),
                      list_model = list(pars = parameters))
estim1 = mtarns(ini_obj = initial,niter = 500,chain = TRUE,burn = 500)
print.regime_model(estim1)
}