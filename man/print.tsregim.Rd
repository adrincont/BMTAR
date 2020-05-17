\name{print.tsregim}
\alias{print.tsregim}

\title{summary tsregim object}

\description{
Produces a ggplot object for the results of the mtarns and mtarstr functions.
}
\usage{
\method{print}{tsregim}(x, ...)
print.tsregim(x)
}

\arguments{
\item{x}{Object of class \dQuote{\code{tsregim}}. Not NULL}
\item{...}{Other parameters that affect print.}
}

\details{
Print the structure of a object class \dQuote{\code{tsregim}}.
}
\value{
  Return to console.
}
\references{
  Calderon V, S. A., & Nieto, F. H. (2017).  \emph{Bayesian analysis of multivariate       threshold autoregressive models with missing data}. Communications in Statistics-Theory   and Methods, 46(1), 296-318.
}
\seealso{\code{\link[MTAR]{tsregim}}}
\author{Valeria Bejarano \email{vbejaranos@unal.edu.co} & Andrey Rincon \email{adrincont@unal.edu.co}}

\examples{
data("datasim")
data = datasim
print.tsregim(data$Sim)
}
