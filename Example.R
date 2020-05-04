#=======================================================================================#
# Example mtaregim:
orders = list(p = 2,q = 1,d = 1)
Phi = list(phi2 = matrix(c(0.1,0.6,-0.4,0.5),2,2, byrow = T))
Beta = list(beta1 = matrix(c(0.3,-0.4),2, 1))
Delta = list(delta1 = matrix(c(0.6,1),2,1))
Sigma = matrix(c(1,0.6,0.6,1.5),2,2,byrow = T)
cs = matrix(c(1,-1),nrow = 2)
Ri = mtaregim(orders = orders,Phi = Phi,Beta = Beta,Delta = Delta,Sigma = Sigma,cs = cs)
#=======================================================================================#
# Example mtarsim
## get Ut data process
Tlen = 500
Sigma_ut = 2
Phi_ut = list(phi1 = 0.3)
R_ut = list(R1 = mtaregim(orders = list(p = 1,q = 0,d = 0),Phi = Phi_ut,Sigma = Sigma_ut))
Ut = mtarsim(N = Tlen,Rg = R_ut,seed = 124)
Zt = Ut$Sim$Yt
# Yt process
k = 2
## R1 regime
Phi_R1 = list(phi1 = matrix(c(0.1,0.6,-0.4,0.5),k,k,byrow = T))
Sigma_R1 = matrix(c(1,0,0,1),k,k,byrow = T)
R1 = mtaregim(orders = list(p = 1,q = 0,d = 0),Phi = Phi_R1,Sigma = Sigma_R1)
## R2 regime
Phi_R2 = list(phi1 = matrix(c(0.3,0.5,0.2,0.7),2,2,byrow = T))
Sigma_R2 = matrix(c(2.5,0.5,0.5,1),2,2,byrow = T)
R2 = mtaregim(orders = list(p = 1,q = 0,d = 0),Phi = Phi_R2,Sigma = Sigma_R2)
## create list of regime-type objects
Rg = list(R1 = R1,R2 = R2)
r = 0.3
## get the simulation
simul = mtarsim(N = Tlen,Rg = Rg,r = r,Zt = Zt,seed = 124)
autoplot.tsregim(simul$Sim,1)
autoplot.tsregim(simul$Sim,2)
#=======================================================================================#
# Example tsregim
data("datasim")
yt = datasim$Sim
Yt = yt$Yt
Zt = yt$Zt
datos = tsregim(Yt,Zt)
#=======================================================================================#
# Example mtarinipars
data("datasim")
tsregim_obj = datasim$Sim
# ns: l siempre conocido, Sigma = NULL = list(R1,R2) puede ser
# conocido r = NULL puede ser conocido
# Sigma conocido y r conocido
parameters = list(l = length(datasim$Reg),
                  Sigma = list(R1 = Sigma_R1,R2 = Sigma_R2),
                  r = tsregim_obj$r,
                  orders = list(pj = datasim$pj, qj = datasim$qj, dj = datasim$dj))
initpars_Sr = mtarinipars(tsregim_obj,list_model = list(pars = parameters))
#r conocido
parameters = list(l = length(datasim$Reg),Sigma = NULL, r = tsregim_obj$r,
                  orders = list(pj = datasim$pj, qj = datasim$qj, dj = datasim$dj))
initpars_r = mtarinipars(tsregim_obj,list_model = list(pars = parameters))
#r desconocido
parameters = list(l = length(datasim$Reg),Sigma = NULL, r = NULL,
                  orders = list(pj = datasim$pj, qj = datasim$qj, dj = datasim$dj))
initpars = mtarinipars(tsregim_obj,list_model = list(pars = parameters))
#str: l siempre conocido
parameters = list(l = length(datasim$Reg))
orders = list(pj = c(2,2),qj = c(1,1),dj = c(1,1))
initpars_KUO = mtarinipars(tsregim_obj,
                           list_model = list(
                             pars = parameters,orders = orders),method = 'KUO')
initpars_SSVS = mtarinipars(tsregim_obj,
                            list_model = list(
                              pars = parameters,orders = orders),method = 'SSVS')
#========================================================================================#
# Example mtarns
data("datasim")
data = datasim
#r known
parameters = list(l = 2,
                  orders = list(pj = c(1,1),dj = c(0,0),qj = c(0,0)),
                  r = data$Sim$r)
initial = mtarinipars(tsregim_obj = data$Sim,
                      list_model = list(pars = parameters))
estim1 = mtarns(ini_obj = initial,niter = 500,chain = TRUE)
print.regim_model(estim1)
autoplot.regim_model(estim1,2)
autoplot.regim_model(estim1,3)
#r unknown
parameters = list(l = 2,orders = list(pj = c(1,1),dj = c(0,0),qj = c(0,0)))
initial = mtarinipars(tsregim_obj = data$Sim,list_model = list(pars = parameters))
estim2 = mtarns(ini_obj = initial,niter = 500,chain = TRUE)
print.regim_model(estim2)
autoplot.regim_model(estim2,1)
autoplot.regim_model(estim2,2)
autoplot.regim_model(estim2,3)
#========================================================================================#
# Example mtarstr
data("datasim")
data = datasim
# Metodo KUO
library(compiler)
initial = mtarinipars(tsregim_obj = data$Sim,method = 'KUO',
                      list_model = list(pars = list(l = 2), orders = list(pj = c(2,2),qj = c(0,0),dj = c(0,0))))
estruc = mtarstr(ini_obj = initial,niter = 500,chain = T)
autoplot.regim_model(estruc,1)
autoplot.regim_model(estruc,2)
autoplot.regim_model(estruc,3)
autoplot.regim_model(estruc,4)
# Metodo SSVS
initial = mtarinipars(tsregim_obj = data$Sim,method = 'SSVS',
                      list_model = list(pars = list(l = 2), orders = list(pj = c(2,2),qj = c(0,0),dj = c(0,0))))
estruc = mtarstr(ini_obj = initial,niter = 500,chain = T)

autoplot.regim_model(estruc,1)
autoplot.regim_model(estruc,2)
autoplot.regim_model(estruc,3)
autoplot.regim_model(estruc,4)

# Example mtarNAIC
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

# Example mtarmissing
data("datasim")
yt = datasim$Sim
# Simulacion de datos faltantes
data_yt = yt$Yt
data_zt = yt$Zt
posNA = sample(c(1:500),8)
data_yt[c(posNA),] = c(NA,NA)
posNA = sample(c(1:500),8)
data_zt[c(posNA)] = NA
data_final = tsregim(data_yt,data_zt,r = yt$r)
autoplot.tsregim(data_final,1)
autoplot.tsregim(data_final,2)

initial = mtarinipars(tsregim_obj = data_final,
                      list_model = list(pars = list(l = 2,r = datasim$Sim$r,
                                                    orders = list(pj = c(1,1), qj = c(0,0),dj = c(0,0)))))
missingest = mtarmissing(ini_obj = initial,chain = TRUE,niter = 500,burn = 500)
print(missingest)
autoplot.regime_missing(missingest,1)
datasim$Sim$Yt[is.na(data_yt[,1]),]

# Example mtarnumreg
data("datasim")
data = datasim
initial = mtarinipars(tsregim_obj = data$Sim,list_model = list(l0 = 3),method = 'KUO')
estim = mtarnumreg(ini_obj = initial,iterprev = 500,niter_m = 500,burn_m = 500, list_m = T,
                   ordersprev = list(maxpj = 2,maxqj = 0,maxdj = 0),parallel = TRUE)
estim$final_m

# Example auto_mtar
data("datasim")
final = auto_mtar(Yt = datasim$Sim$Yt, Zt = datasim$Sim$Zt, l0 = 3, maxorders = list(pj = 2,qj = 0,dj = 0),
                  niter = 5000, chain = FALSE, method = 'KUO')
