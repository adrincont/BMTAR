\name{print.regim_model}
\alias{print.regim_model}

\title{print regim_model object for the function outputs mtarns and mtastr}

\description{
Print estimates for the results of the mtarns and mtarstr functions.
}
\usage{
\method{print}{regim_model}(object, ...)
print.regim_model(object)
}

\arguments{
\item{object}{Object of class \dQuote{\code{regim_model}}.}
\item{...}{Other print parameters that affect.}
}

\details{
Print estimates outputs corresponding to the functions \dQuote{\code{mtarns}} and \dQuote{\code{mtarstr}} which return an object of class \dQuote{\code{regim_model}}.
}
\value{
  Return to console.
}
\references{
  Calderon V, S. A., & Nieto, F. H. (2017).  \emph{Bayesian analysis of multivariate       threshold autoregressive models with missing data}. Communications in Statistics-Theory   and Methods, 46(1), 296-318.
}
\seealso{\code{\link[MTAR]{mtarns}}, \code{\link[MTAR]{mtarstr}}}
\author{Valeria Bejarano \email{vbejaranos@unal.edu.co} & Andrey Rincon \email{adrincont@unal.edu.co}}

\examples{
data("datasim")
data = datasim
#r known
parameters = list(l = 2,
                  orders = list(pj = c(1,1),dj = c(0,0),qj = c(0,0)),
                  r = data$Sim$r)
initial = mtarinipars(tsregim_obj = data$Sim,
                      list_model = list(pars = parameters))
estim1 = mtarns(ini_obj = initial,niter = 500,chain = TRUE,burn = 500)
print.regim_model(estim1)
}
