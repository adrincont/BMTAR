# Date: 29/07/2019
# Description:
#-> Calculo para el criterio NAIC para un abjeti tipo "regime-model".
# Function:
# Date: 29/07/2019
# Description:
#-> Calculo para el criterio NAIC para un abjeti tipo "regime-model".
# Function:
mtarNAIC = function(regimemodel){
  if (class(regimemodel) != "regime-model") {
    stop("regimemodel must be an object of type (regime-model)")
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
  cat("NAIC=",round(NAIC,4),"\n")
  return(list(AICj = AICj,NAIC = NAIC))
}
# Example:
yt = datasim #(Buscar datos)
estim1 = mtarns(Yt = yt$Yt,Ut = Ut,l = 2,orders = list(pj = yt$pj,dj = yt$dj,qj = yt$qj)
                ,r = qnorm(0.4),niter = 300)
estim2 = mtarns(Yt = yt$Yt,Ut = Ut,l = 2,orders = list(pj = yt$pj,dj = yt$dj,qj = yt$qj)
                ,r = qnorm(0.1),niter = 300)
mtarNAIC(estim1)
mtarNAIC(estim2)
rprop = function(r){
  estim = mtarns(Yt = yt$Yt,Ut = Ut,l = 2,orders = list(pj = yt$pj,dj = yt$dj,qj = yt$qj)
                  ,r = r,niter = 300)
  naic = mtarNAIC(estim)
  return(naic$NAIC)
}
rprop = Vectorize(rprop)
r = seq(-1,1,length.out = 50)
rnaic = rprop(r)
ggplot2::ggplot(data = NULL, ggplot2::aes(x = r,y = rnaic)) + ggplot2::geom_point() + 
  ggplot2::geom_line() + ggplot2::theme_bw() + ggplot2::geom_vline(xintercept = qnorm(0.4), color = 2)

estruc1 = mtarstr(Yt = yt$Yt,Ut = Ut,l = 2,orders = list(pjmax = c(2,2),qjmax = c(1,1),djmax = c(1,1)),niter = 100,method = 'KUO',chain = T)
estruc2 = mtarstr(Yt = yt$Yt,Ut = Ut,l = 2,orders = list(pjmax = c(1,1),qjmax = c(1,1),djmax = c(1,1)),niter = 100,method = 'KUO',chain = T)
estruc3 = mtarstr(Yt = yt$Yt,Ut = Ut,l = 3,orders = list(pjmax = c(2,2,1),qjmax = c(1,1,1),djmax = c(1,1,1)),niter = 100,method = 'KUO',chain = T)
estruc4 = mtarstr(Yt = yt$Yt,Ut = Ut,l = 4,orders = list(pjmax = c(2,2,1,1),qjmax = c(1,1,1,1),djmax = c(1,1,1,1)),niter = 100,method = 'KUO',chain = T)
mtarNAIC(estruc1)
mtarNAIC(estruc2)
mtarNAIC(estruc3)
mtarNAIC(estruc4)
