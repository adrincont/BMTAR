\name{print.tsregime}
\alias{print.tsregime}

\title{Print tsregime object}

\description{
Print the structure of a object class \dQuote{\code{tsregime}}.
}
\usage{
\method{print}{tsregime}(object, ...)
}

\arguments{
\item{object}{Object of class \dQuote{\code{tsregime}}. Not NULL}
\item{...}{Other parameters that affect print.}
}

\details{
Print the structure of a object class \dQuote{\code{tsregime}}.
}
\value{
  Return to console.
}
\references{
  Calderon V, S. A., & Nieto, F. H. (2017).  \emph{Bayesian analysis of multivariate       threshold autoregressive models with missing data}. Communications in Statistics-Theory   and Methods, 46(1), 296-318.
}
\seealso{\code{\link[MTAR]{tsregime}}}
\author{Valeria Bejarano \email{vbejaranos@unal.edu.co} & Andrey Rincon \email{adrincont@unal.edu.co}}

\examples{
data("datasim")
data = datasim
print.tsregime(data$Sim)
}
