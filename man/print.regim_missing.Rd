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
Print the estimates for the outputs corresponding to the function \dQuote{\code{mtarmissing}} which return an object of class \dQuote{\code{regim_missing}}.
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
data('missingest')
print.regim_missing(missingest)
}
