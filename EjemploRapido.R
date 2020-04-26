rm(list = ls())
library(MTAR)
library(ggplot2)
Tlen = 500

sigmarZ = matrix(c(0.5),1,1)
AZ = matrix(c(0.4),1,1,byrow = T)
R_Zt = list(R1 = mtaregim(orders = list(p = 1,q = 0,d = 0),Phi = list(phi1 = AZ),Sigma = sigmarZ))
Zt = mtarsim(N = Tlen,Rg = R_Zt)
print(Zt$Sim)

R1 = mtaregim(orders = list(p = 1,q = 0,d = 0),
              Phi = list(phi1 = matrix(c(0.5,-0.2,-0.2,0.8),2,2,byrow = T)),
              Sigma = matrix(c(1,0.6,0.6,1.5),2,2,byrow = T))
R2 = mtaregim(orders = list(p = 1,q = 0,d = 0),
              Phi = list(phi1 = matrix(c(0.3,0.5,0.2,0.7),2,2,byrow = T)),
              Sigma = matrix(c(2.5,0.5,0.5,1),2,2,byrow = T),
              cs = matrix(c(5,2),nrow = 2))
datasim = mtarsim(N = Tlen,Rg = list(R1 = R1,R2 = R2),r = qnorm(0.4), Xt = NULL, Zt = Zt$Sim$Yt)
data_cm = datasim$Sim
autoplot(data_cm,1)
autoplot(data_cm,2)
autoplot(data_cm,3)

initial = mtarinipars(tsregim_obj = data_cm,list_model = list(l0 = 3),method = 'KUO')
numregest_KUO = mtarnumreg(ini_obj = initial,niter_m = 500,chain_m = T,list_m = T,iterprev = 500,
                           ordersprev = list(maxpj = 2, maxqj = 0,maxdj = 0))
numregest_KUO$estimates
initial = mtarinipars(tsregim_obj = data_cm,method = 'KUO',
                      list_model = list(pars = list(l = numregest_KUO$final_m),
                                        orders = list(pj = c(2,2),qj = c(0,0),dj = c(0,0))))
est_final_KUO = mtarstr(ini_obj = initial,niter = 2000,chain = T)
print(est_final_KUO)
est_final_KUO$regime
autoplot(est_final_KUO,1)
autoplot(est_final_KUO,2)
autoplot(est_final_KUO,3)
autoplot(est_final_KUO,4)

initial = mtarinipars(tsregim_obj = data_cm,list_model = list(l0 = 3),method = 'SSVS')
numregest_SSVS = mtarnumreg(ini_obj = initial,NAIC = T, iterprev = 500,
                            ordersprev = list(maxpj = 2, maxqj = 0,maxdj = 0))
numregest_SSVS$NAIC_final_m
initial = mtarinipars(tsregim_obj = data_cm,method = 'SSVS',
                      list_model = list(pars = list(l = numregest_SSVS$NAIC_final_m),
                                        orders = list(pj = c(2,2),qj = c(0,0),dj = c(0,0))))
est_final_SSVS = mtarstr(ini_obj = initial,niter = 2000,chain = T)
print(est_final_SSVS)
est_final_SSVS$regime
autoplot(est_final_SSVS,1)
autoplot(est_final_SSVS,2)
autoplot(est_final_SSVS,3)
autoplot(est_final_SSVS,4)

save.image('EjemploRapido.RData')
