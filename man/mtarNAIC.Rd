\name{mtarNAIC}
\alias{mtarNAIC}

\title{Compute NAIC of a MTAR model}
\description{
Compute the Non-linear Akaike information criterion (NAIC) of a ``regime_model'' class object.
}
\usage{
mtarNAIC(regimemodel)
}

\arguments{
  \item{regimemodel}{object of class ``regime_model''}
}

\details{
Estimation of thresholds was made before starting the Bayesian procedure via the Non-linear Akaike information criterion (NAIC) (Tong, 1990), in MTAR model. The NAIC for a MTAR model with l regimes is:
 \deqn{
 NAIC = {\sum_{j = 1}^l AICj(r)}/{\sum_{j = 1}^l N_j}
 }
  \deqn{
  AICj(r) = Nj ln(|Sj/Nj|) + 2k \eta_j
 }
 Nj: number of observations in each regime.
 \deqn{
 S_{j}=\sum_{t:j_t = j}(y_t-YP)'(y_t-YP)
 }
  \deqn{
  YP =  \Phi_{0}^(j)+\sum_{i=1}^{p_j}\Phi_{i}^{(j)} Y_{t-i}+\sum_{i=1}^{q_j} \beta_{i}^{(j)} X_{t-i} + \sum_{i=1}^{d_j} \delta_{i}^{(j)} Z_{t-i}
 }
}
\value{
  Return a list type object
  $AICj numeric type, AIC for each regime
  $NAIC numeric type, NAIC value
}
\references{
  Calderon V, S. A., & Nieto, F. H. (2017).  \emph{Bayesian analysis of multivariate       threshold autoregressive models with missing data}. Communications in Statistics-Theory   and Methods, 46(1), 296-318.
}
\author{Valeria Bejarano \email{vbejaranos@unal.edu.co}, Sergio Calderón \email{sacalderonv@unal.edu.co} & Andrey Rincon \email{adrincont@unal.edu.co}}

\examples{
data("datasim")
data = datasim
# Estimate number of regimes with NAIC
initial2 = mtarinipars(tsregim_obj = data$Sim,method = 'KUO',
list_model = list(pars = list(l = 2),
orders = list(pj = c(1,1),qj = c(1,1),dj = c(1,1))))
estruc2 = mtarstr(ini_obj = initial2,niter = 100,chain = TRUE,burn = 100)
initial3 = mtarinipars(tsregim_obj = data$Sim,method = 'KUO',
list_model = list(pars = list(l = 3),
orders = list(pj = c(1,1,1),qj = c(1,1,1),dj = c(1,1,1))))
estruc3 = mtarstr(ini_obj = initial3,niter = 100,chain = TRUE,burn = 100)
mtarNAIC(estruc2)
mtarNAIC(estruc3)
}
\keyword{AIC}
\keyword{NAIC}
