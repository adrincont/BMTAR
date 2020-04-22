# Date: 14/04/2020
# Description:
# Function:
mtarNAIC = function(regimemodel){
  if (class(regimemodel) != 'regim_model') {
    stop('regimemodel must be an object of type (regim_model)')
  }
  l = length(regimemodel$regime)
  k = nrow(regimemodel$regime[[1]]$sigma)
  nuaux = NULL
  for (lj in 1:l) {
    nuaux[lj] = length(regimemodel$regime[[lj]]$beta[[1]][1,])
  }
  nu = max(nuaux)
  pj = qj = dj = vector(length = l)
  for (lj in 1:l) {
    pj[lj] = length(regimemodel$regime[[lj]]$phi)
    qj[lj] = length(regimemodel$regime[[lj]]$beta)
    dj[lj] = length(regimemodel$regime[[lj]]$delta)
  }
  etaj = 1 + pj*k + qj*nu + dj
  Nj = c(regimemodel$Nj)
  logLikj = as.numeric(regimemodel$logLikj)
  AICj = as.numeric(Nj*logLikj + 2*k*etaj,nrow = 1,row.names = NULL)
  NAIC = sum(AICj)/sum(Nj)
  cat('NAIC=',round(NAIC,4),'\n')
  return(list(AICj = AICj,NAIC = NAIC))
}
# Example:
data("datasim")
yt = datasim$Sim
parameters = list(l = 2,orders = list(pj = c(1,1),dj = c(0,0),qj = c(0,0)),r = yt$r)
initial = mtarinipars(tsregim_obj = yt,list_model = list(pars = parameters))
estim1 = mtarns(ini_obj = initial,niter = 500,chain = TRUE,burn = 200)
parameters = list(l = 2,orders = list(pj = c(1,1),dj = c(0,0),qj = c(0,0)),r = qnorm(0.1))
initial = mtarinipars(tsregim_obj = yt,list_model = list(pars = parameters))
estim2 = mtarns(ini_obj = initial,niter = 500,chain = TRUE,burn = 200)
mtarNAIC(estim1)
mtarNAIC(estim2)
rprop = function(r){
  parameters = list(l = 2,orders = list(pj = c(1,1),dj = c(0,0),qj = c(0,0)),r = r)
  initial = mtarinipars(tsregim_obj = yt,list_model = list(pars = parameters))
  estim = mtarns(ini_obj = initial,niter = 500,chain = TRUE,burn = 200)
  naic = mtarNAIC(estim)
  return(naic$NAIC)
}
rprop = Vectorize(rprop)
r = seq(-1,1,length.out = 50)
rnaic = rprop(r)
ggplot2::ggplot(data = NULL, ggplot2::aes(x = r,y = rnaic)) + ggplot2::geom_point() +
  ggplot2::geom_line() + ggplot2::theme_bw() + ggplot2::geom_vline(xintercept = yt$r, color = 2)


initial1 = mtarinipars(tsregim_obj = data$Sim,method = 'KUO',
                      list_model = list(pars = list(l = 2), orders = list(pj = c(2,2),qj = c(1,1),dj = c(1,1))))
estruc1 = mtarstr(ini_obj = initial1,niter = 100,chain = T,burn = 100)
initial2 = mtarinipars(tsregim_obj = data$Sim,method = 'KUO',
                       list_model = list(pars = list(l = 2), orders = list(pj = c(1,1),qj = c(1,1),dj = c(1,1))))
estruc2 = mtarstr(ini_obj = initial2,niter = 100,chain = T,burn = 100)
initial3 = mtarinipars(tsregim_obj = data$Sim,method = 'KUO',
                       list_model = list(pars = list(l = 3), orders = list(pj = c(2,2,1),qj = c(1,1,1),dj = c(1,1,1))))
estruc3 = mtarstr(ini_obj = initial3,niter = 100,chain = T,burn = 100)
initial4 = mtarinipars(tsregim_obj = data$Sim,method = 'KUO',
                       list_model = list(pars = list(l = 4), orders = list(pj = c(2,2,1,1),qj = c(1,1,1,1),dj = c(1,1,1,1))))
estruc4 = mtarstr(ini_obj = initial4,niter = 100,chain = T,burn = 100)
mtarNAIC(estruc1)
mtarNAIC(estruc2)
mtarNAIC(estruc3)
mtarNAIC(estruc4)
