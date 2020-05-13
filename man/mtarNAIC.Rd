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
  AICj(r) = Nj*ln(|Sj/Nj|) + 2k\eta_j
 }
}
\value{
  Return a list type object
  $AICj numeric type, AIC for each regime
  $NAIC numeric type, NAIC value
}
\references{
  Fox, J. (2016)
  \emph{Applied Regression Analysis and Generalized Linear Models},
  Third Edition. Sage.
}
\author{Valeria Bejarano \email{vbejaranos@unal.edu.co} & Andrey Rincon \email{adrincont@unal.edu.co}}

\examples{
  poner ejemplo de uso de funcion con datos reales
}

\keyword{AIC}
\keyword{NAIC}
