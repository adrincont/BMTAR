\name{autoplot.tsregim}
\alias{autoplot.tsregim}

\title{tsregim object ggplot for the outputs on the function tsregim}

\description{
Produces a ggplot object for the results of the tsregim function.
}
\usage{
\method{autoplot}{tsregim}(object, type = 1)
}

\arguments{
\item{object}{Object of class \dQuote{\code{tsregim}}.}
\item{type}{character string giving the type of plot to be computed. Allowed
values are 1 for "\code{Output processes}" (the default), 2 for \dQuote{\code{Threshold processes}}, 3 for \dQuote{\code{Covariates processes}}
}.
}
\details{
Graph for the stochastic processes of the object \ dQuote {\ code {tsregim}} for: Output process and if they exist for: Threshold process and Covariates process. In the case that there are missing data, a red line is drawn at the indicated time.
}
\value{
  Return a ggplot object.
}
\references{
  Calderon V, S. A., & Nieto, F. H. (2017).  \emph{Bayesian analysis of multivariate       threshold autoregressive models with missing data}. Communications in Statistics-Theory   and Methods, 46(1), 296-318.
}
\seealso{\code{\link[MTAR]{tsregim}}}
\author{Valeria Bejarano \email{vbejaranos@unal.edu.co} & Andrey Rincon \email{adrincont@unal.edu.co}}

\examples{
data("datasim")
yt = datasim$Sim
Yt = yt$Yt
Zt = yt$Zt
data = tsregim(Yt,Zt)
autoplot.tsregim(data)
}
