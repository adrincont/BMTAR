\name{tsregim}
\alias{tsregim}

\title{Creation of class ``tsregim'' for some data}
\description{
The function tsregim is used to create time-series-regimen objects.
}
\usage{
tsregim(Yt, Zt = NULL, Xt = NULL, r = NULL, ...)
autoplot(tsregim)
print(tsregim)
}

\arguments{
  \item{Yt}{matrix (Nxk) type object, observed process (admit NA values). Not NULL}
  \item{Zt}{matrix (Nx1) type object, threshold process (admit NA values). Default NULL}
  \item{Xt}{matrix (Nxnu) type object, covariate process (admit NA values). Default NULL}
  \item{r}{numeric type, threshold value (within the range of Zt) if known. Default NULL}
}

\details{
 Create a class ``tsregim'' object composed of: \eqn{\left\{Y}_{t}\right\}} and \eqn{\left\{\mathrm{X}_{t}\right\}} stochastics processes such that \eqn{\mathrm{Y}_{t}=\left(\mathrm{Y}_{1 t}, \ldots, \mathrm{Y}_{k t}\right)^{\prime}}, \eqn{\mathrm{X}_{t}=\left(\mathrm{X}_{1 t}, \ldots, \mathrm{X}_{v t}\right)^{\prime}} and \eqn{\left\{\mathrm{Z}_{t}\right\}} is a univariate process. Where \eqn{\left\{\mathrm{Y}_{t}\right\}} follows a MTAR model with threshold variable \eqn{\mathrm{Z}_{t}}
  \eq{
  \begin{equation}\mathrm{Y}_{t}=\phi_{0}^{(j)}+\sum_{i=1}^{p_{j}} \phi_{i}^{(j)} \mathrm{Y}_{t-i}+\sum_{i=1}^{q_{j}} \boldsymbol{\beta}_{i}^{(j)} \mathrm{X}_{t-i}+\sum_{i=1}^{d_{j}} \boldsymbol{\delta}_{i}^{(j)} \mathrm{Z}_{t-i}+\mathbf{\Sigma}_{(j)}^{1 / 2} \varepsilon_{t} \text { if } r_{j-1}<\mathrm{Z}_{t} \leq r_{j}\end{equation}
  }
Missing data is allowed for processes \eqn{\left\{Y}_{t}\right\}}, \eqn{\left\{\mathrm{X}_{t}\right\}} and \eqn{\mathrm{Z}_{t}} (can then be estimated with the `` mtarmissing '' function). In the case of known r, the output returns the percentages of observations found in each regimen.
}
\value{
  Return a list type object of class ``tsregim'': $Yt $Xt $Zt stochastics processes (if enter), $N number of observations, $k number of variables. If r known: $r threshold value,$Ind numeric type, number of the regime each observation belongs,$Summary_r data.frame type, number and proportion of observations in each regime.
}
\references{
  Calderon V, S. A., & Nieto, F. H. (2017).  \emph{Bayesian analysis of multivariate       threshold autoregressive models with missing data}. Communications in Statistics-Theory   and Methods, 46(1), 296-318.
}
\author{Valeria Bejarano \email{vbejaranos@unal.edu.co} & Andrey Rincon \email{adrincont@unal.edu.co}}

\seealso{\code{\link[MTAR]{mtaregim}}, \code{\link[MTAR]{mtarns}},
  \code{\link[MTAR]{mtarstr}}, \code{\link{mtarmissin}},code{\link{mtarnumreg}}}

\examples{
data("datasim")
yt = datasim$Sim
Yt = yt$Yt
Zt = yt$Zt
(datos = tsregim(Yt,Zt))
autoplot.tsregim(datos,1)
autoplot.tsregim(datos,2)
}

\keyword{Observed process}

