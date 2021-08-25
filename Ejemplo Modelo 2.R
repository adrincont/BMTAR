# Ejemplos de Textos
## Ut = (Zt,Xt)
Tlen = 2000
Sigma_ut = expm::sqrtm(matrix(c(1,0.4,0.2,
                                0.4,1,0.4,
                                0.2,0.4,1),3,3,byrow = TRUE))
Phi_ut = list(phi1 = matrix(c(0.5,0,0,
                              0.1,0.1,0.3,
                              0,0.2,0.3),3,3))
R_ut = list(R1 = mtaregime(orders = list(p = 1,q = 0,d = 0),Phi = Phi_ut,Sigma = Sigma_ut))
Ut = mtarsim(N = Tlen,Rg = R_ut,seed = 1997)
autoplot(Ut$Sim)
## Modelo 2 ===============================================================.####
k = 3
nu = 2
## R1 regime
Phi_R1 = list(phi1 = matrix(c(0.1,0.2,0,
                              0.1,0.1,0.3,
                              0,0.2,0.5),k,k),
              phi2 = matrix(c(0.3,0.1,0,
                              0,0.1,0.2,
                              0,0.3,0.3),k,k))
Beta_R1 = list(beta1 = matrix(c(0.3,-0.1,0,
                                0.1,0,0.2),k,nu))
Delta_R1 = list(delta1 = matrix(c(0.6,0,0.1),k,1))
Sigma_R1 = matrix(c(1,0.1,0.2,
                    0.1,1,0.4,
                    0.2,0.4,1),k,k,byrow = TRUE)
R1 = mtaregime(orders = list(p = 2,q = 1,d = 1),Phi = Phi_R1,Beta = Beta_R1,Delta = Delta_R1,Sigma = Sigma_R1,cs = c(1,0,-1))
## R2 regime
Phi_R2 = list(phi1 = matrix(c(0.1,0.2,0,
                              0,0.1,0.3,
                              0.2,0.2,0.5),k,k))
Sigma_R2 = matrix(c(1,0,0.2,
                    0,1,0.3,
                    0.2,0.3,1),k,k,byrow = TRUE)
R2 = mtaregime(orders = list(p = 1,q = 0,d = 0),Phi = Phi_R2,Sigma = Sigma_R2,cs = c(5,2,0))
## R3 regime
Phi_R3 = list(phi1 = matrix(c(0.1,0.2,0,
                              0,0.1,0.3,
                              0.2,0,0.4),k,k))
Beta_R3 = list(beta1 = matrix(c(0.4,0.05,0,
                                -0.1,0,0.2),k,nu))
Sigma_R3 = matrix(c(1,0,0.1,
                    0,1,0.4,
                    0.1,0.4,1),k,k,byrow = TRUE)
R3 = mtaregime(orders = list(p = 1,q = 1,d = 0),Phi = Phi_R3,Beta = Beta_R3,Sigma = Sigma_R3,cs = c(5,2,0))
## create list of regime-type objects
Rg = list(R1 = R1,R2 = R2,R3 = R3)
r = quantile(Ut$Sim$Yt[,1],c(0.4,0.7))
# Simulacion
datasim = mtarsim(N = Tlen,Rg = Rg,r = r,Zt = Ut$Sim$Yt[,1],Xt = Ut$Sim$Yt[,-1],seed = 1997)
autoplot.tsregime(datasim$Sim,1)
autoplot.tsregime(datasim$Sim,2)
autoplot.tsregime(datasim$Sim,3)
### Simular datos Faltantes -----------------------------------------------.####
pos_Nas_Y_1 = sample(1:2000,20)
pos_Nas_Y_3 = sample(1:2000,20)
pos_Nas_X_1 = sample(1:2000,10)
pos_Nas_Z = sample(1:2000,5)

Yt_sim = datasim$Sim$Yt
Yt_sim[pos_Nas_Y_1,1] = NA
Yt_sim[pos_Nas_Y_3,3] = NA
Xt_sim = datasim$Sim$Xt
Xt_sim[pos_Nas_X_1,1] = NA
Zt_sim = datasim$Sim$Zt
Zt_sim[pos_Nas_Z,1] = NA

data_miss = tsregime(Yt = Yt_sim,Zt = Zt_sim,Xt = Xt_sim)
print.tsregime(data_miss)
autoplot.tsregime(data_miss,1)
autoplot.tsregime(data_miss,2)
autoplot.tsregime(data_miss,3)
### Estimacion de parametros y datos faltantes ----------------------------.####
#### 1) Completamos con promedio
Y_temp = data_miss$Yt
Y_temp[,1][is.na(Y_temp[,1])] = mean(Y_temp[,1],na.rm = TRUE)
Y_temp[,3][is.na(Y_temp[,3])] = mean(Y_temp[,3],na.rm = TRUE)
X_temp = data_miss$Xt
X_temp[,1][is.na(X_temp[,1])] = mean(X_temp[,1],na.rm = TRUE)
Z_temp = data_miss$Zt
Z_temp[,1][is.na(Z_temp[,1])] = mean(Z_temp[,1],na.rm = TRUE)
### 2) Estimar numero de regimenes (Con datos completados por media)
data_temp = tsregime(Y_temp,Z_temp,X_temp)
initial = mtarinipars(tsregime_obj = data_temp,list_model = list(l0_max = 3),
                      method = 'KUO')
estim_nr = mtarnumreg(ini_obj = initial,iterprev = 1000,niter_m = 2000,burn_m = 1000,
                      list_m = TRUE,#
                      ordersprev = list(maxpj = 2,maxqj = 2,maxdj = 2),
                      parallel = FALSE)
print(estim_nr)

autoplot(estim_nr$list_m$m2$par,1)
autoplot(estim_nr$list_m$m2$par,2)
autoplot(estim_nr$list_m$m2$par,3)
autoplot(estim_nr$list_m$m2$par,4)
autoplot(estim_nr$list_m$m2$par,5)
diagnostic_mtar(estim_nr$list_m$m2$par)

autoplot(estim_nr$list_m$m3$par,1)
autoplot(estim_nr$list_m$m3$par,2)
autoplot(estim_nr$list_m$m3$par,3)
autoplot(estim_nr$list_m$m3$par,4)
autoplot(estim_nr$list_m$m3$par,5)
diagnostic_mtar(estim_nr$list_m$m3$par)
