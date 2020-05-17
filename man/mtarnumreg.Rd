\name{mtarnumreg}
\alias{mtarnumreg}

\title{Estimation of the number of regimes in a MTAR model}
\description{
Compute estimation of number of regimes by NAIC or Carlin and Chip methodology for a MTAR model
}
\usage{
mtarnumreg(ini_obj, r_init = NULL, level = 0.95, burn_m = NULL,
niter_m = 1000,iterprev = 500, chain_m = FALSE, list_m = FALSE,
NAIC = FALSE,ordersprev = list(maxpj = 2,maxqj = 0,maxdj = 0),
parallel = FALSE)
}

\arguments{
  \item{ini_obj}{class ``regime_inipars'' object, here specificate only l0. Not NULL}
  \item{r_init}{numeric type of length l - 1. If r not known, starting value of the chain. Default NULL}
  \item{level}{numeric type, confident interval for estimations. Default 0.95}
  \item{burn_m}{numeric type, number of initial runs. Default NULL (10\% of niter)}
  \item{niter_m}{numeric type, number of runs of MCMC. Default 1000}
  \item{iterprev}{numeric type, number of runs for pseudo values. Default 500}
  \item{chain_m}{logical type, if return chains of parameters. Default FALSE}
  \item{list_m}{logical type, if return list of regimes considered. Default FALSE}
  \item{NAIC}{logical type, if return estimation of number of regimes by NAIC (not run Carlin and Chip for l). Default FALSE}
  \item{ordersprev}{list type object with names (maxpj,maxqj,maxdj), maximum number of lags of each process consider in the pseudo values for each number of regimes considered . Default maxpj = 2,maxqj = 0, maxdj = 0}
  \item{parallel}{logical type, if package \code{parallel} should be used. Default FALSE}
}

\details{
  Detalles de los metodos usados en la funcion corespondiente
}
\value{
  Return a list type object of class ``regime_number''
  \item{tsregim}{ini_obj$tsregim_obj}
  \item{list_m}{if list_m TRUE list of models considered}
  \item{m_chain}{if chain_m TRUE chains of m}
  \item{estimates}{table of the proportions of m estimated}
  \item{final_m}{numeric type, final number of regimes estimated}

  if NAIC TRUE
  \item{tsregim}{ini_obj$tsregim_obj}
  \item{list_m}{list of consider models}
  \item{NAIC}{list type of NAIC for each considered model}
  \item{NAIC_final_m}{numeric type, final number of regimes by this criteria}
}
\references{
  Calderon V, S. A., & Nieto, F. H. (2017).  \emph{Bayesian analysis of multivariate       threshold autoregressive models with missing data}. Communications in Statistics-Theory   and Methods, 46(1), 296-318.
}
\author{Valeria Bejarano \email{vbejaranos@unal.edu.co} & Andrey Rincon \email{adrincont@unal.edu.co}}

\examples{
data("datasim")
data = datasim
initial = mtarinipars(tsregim_obj = data$Sim,
list_model = list(l0 = 3),method = 'KUO')
estim = mtarnumreg(ini_obj = initial,iterprev = 500,niter_m = 500,
burn_m = 500, list_m = T,ordersprev = list(maxpj = 2,maxqj = 0,maxdj = 0))
estim$final_m
}

\keyword{Bayesian estimation}
\keyword{Carling and Chip}

