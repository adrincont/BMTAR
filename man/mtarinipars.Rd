\name{mtarinipars}
\alias{mtarinipars}

\title{Organization and check model specification}
\description{
Model specification of data, known or unknown parameters and its prior distributions
}
\usage{
mtarinipars(tsregim_obj, list_model = list(pars = list(l = 2, orders = list(pj = c(1,1), qj = c(0,0), dj = c(0,0)), r = NULL, Sigma = NULL),orders = NULL,l0 = NULL), method = NULL, theta_prior = NULL, sigma_prior = NULL, gamma_prior = NULL, r_prior = NULL)
}

\arguments{
  \item{tsregim_obj}{class ``tsregim'' object. Not NULL}
  \item{list_model}{list type object with at least one of the names (pars, orders, l0).
  pars list type object with one of the names (l, orders, r, Sigma) for parameters of the model known
    l number of regimes of the model (not known use l0). Default 2
    orders list type object with names (pj,qj,dj) each a vector of length l, specificate lags orders in each regime. Default list(pj = c(1,1), qj = c(0,0), dj = c(0,0))
    r threshold value. Default NULL
    Sigma list type object with names (R1, ..., Rl) each a matrix type object, specification of error covariate matrix. Default NULL
  orders list type object with names (pj,qj,dj) each a vector of length l,specificate maximum lags orders in each regime if not known. Default NULL
  l0 numeric type, number maximum of regimes to consider for the model if not known. Default NULL
  }
  \item{method}{character type object, if orders not known or enter l0 it must be one ``KUO'' or ``SSVS'', when known it is not necesary. Default NULL}
  \item{theta_prior}{list type object with names (R1, ..., Rl), each one a list type object with at least one of the names (theta0j,cov0j)(if method ``SSVS'' (theta0j,cov0j,Cij,Tauij,R)),specification of mean and covariate matrix of the prior distribution for $\theta$ parameter. Default NULL}
  \item{sigma_prior}{list type object with names (R1, ..., Rl), each one a list type object with at least one of the names (S0j,nu0j) specification of matrix and degrees of freedom of the prior distribution for $\Sigma$ parameter. Default NULL}
  \item{gamma_prior}{list type object with names (R1, ..., Rl), each one a vector of prior probabilities for $\gamma$ parameter. Default NULL}
  \item{r_prior}{list type object with at least one name (za, zb, val_rmh), each one a numeric type object with the minimum, maximum value for r and its parameter for Metropolis-Hasting respectively. Default NULL}
}

\details{
  Detalles de los metodos usados en la funcion corespondiente
}
\value{
  Return a list type object of class ``regim_inipars''
  $tsregim_obj = tsregim_obj
  $pars = list_model$pars
  $orders = list_model$orders
  $method = method
  $init$r = r_prior
  init$Theta = theta_prior
  $init$Sigma = sigma_prior
  $init$Gamma = gamma_prior
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

\keyword{prior distribution}
\keyword{Metropolis - Hasting}


