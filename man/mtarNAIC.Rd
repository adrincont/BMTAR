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
 Nj Corresponds number observations in each regime.
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
\author{Valeria Bejarano \email{vbejaranos@unal.edu.co} & Andrey Rincon \email{adrincont@unal.edu.co}}

\examples{
  data("datasim")
yt = datasim$Sim
# Variable threshold estimation using NAIC (Using parallel package)

## Number of cores we have available
micluster = parallel::makeCluster(4)
doParallel::registerDoParallel(micluster)
parallel::clusterEvalQ(micluster, library(MTAR))
parallel::clusterExport(cl = micluster,varlist = list('yt'))

rprop = function(r){
  parameters = list(l = 2,orders = list(pj = c(1,1),dj = c(0,0),qj = c(0,0)),r = r)
  initial = mtarinipars(tsregim_obj = yt,list_model = list(pars = parameters))
  estim = mtarns(ini_obj = initial,niter = 1000,chain = TRUE,burn = 500)
  naic = mtarNAIC(estim)
  return(naic$NAIC)
}

rprop = Vectorize(rprop)
r = seq(-0.5,0.5,length.out = 40)

## divide the iterations in each cluster
clus_r = list(r[1:10],r[11:20],r[21:30],r[31:40])
rNAIC_parallel = parallel::parLapply(micluster,clus_r,rprop)

## close the parallel process
parallel::stopCluster(micluster)
rnaic = unlist(rNAIC_parallel)

cat('r = ',r[rnaic == min(rnaic)])
ggplot2::ggplot(data = NULL, ggplot2::aes(x = r,y = rnaic)) + ggplot2::geom_point() +
  ggplot2::geom_line() + ggplot2::theme_bw() +
  ggplot2::geom_vline(xintercept = r[rnaic == min(rnaic)], color = 2)

#===============================================================================================#
# Estimate number of regimes with NAIC
    initial2 = mtarinipars(tsregim_obj = data$Sim,method = 'KUO',
                           list_model = list(pars = list(l = 2), orders = list(pj = c(1,1),qj = c(1,1),dj = c(1,1))))
    estruc2 = mtarstr(ini_obj = initial2,niter = 100,chain = T,burn = 100)
    initial3 = mtarinipars(tsregim_obj = data$Sim,method = 'KUO',
                           list_model = list(pars = list(l = 3), orders = list(pj = c(1,1,1),qj = c(1,1,1),dj = c(1,1,1))))
    estruc3 = mtarstr(ini_obj = initial3,niter = 100,chain = T,burn = 100)
    initial4 = mtarinipars(tsregim_obj = data$Sim,method = 'KUO',
                           list_model = list(pars = list(l = 4), orders = list(pj = c(1,1,1,1),qj = c(1,1,1,1),dj = c(1,1,1,1))))
    estruc4 = mtarstr(ini_obj = initial4,niter = 100,chain = T,burn = 100)
mtarNAIC(estruc2)
mtarNAIC(estruc3)
mtarNAIC(estruc4)
}

\keyword{AIC}
\keyword{NAIC}
