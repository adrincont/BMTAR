# Ejemplos de Textos
## Ut = (Zt,Xt)
Tlen = 1000
Sigma_ut = expm::sqrtm(Sigma_ut)
Phi_ut = list(phi1 = matrix(c(0.5,0.4,0.1,0.5),2,2))
R_ut = list(R1 = mtaregime(orders = list(p = 1,q = 0,d = 0),Phi = Phi_ut,Sigma = Sigma_ut))
Ut = mtarsim(N = Tlen,Rg = R_ut,seed = 1997)
## Modelo 1 ===============================================================.####
### Yt process ------------------------------------------------------------.####
k = 2
## R1 regime
Phi_R1 = list(phi1 = matrix(c(0.5,-0.2,-0.2,0.8),k,k),phi2 = matrix(c(0.1,-0.4,0.6,0.5),k,k))
Beta_R1 = list(beta1 = matrix(c(0.3,-0.4),k,1))
Delta_R1 = list(delta1 = matrix(c(0.6,1),k,1))
Sigma_R1 = matrix(c(1,0.6,0.6,1.5),k,k,byrow = TRUE)
R1 = mtaregime(orders = list(p = 2,q = 1,d = 1),Phi = Phi_R1,Beta = Beta_R1,Delta = Delta_R1,Sigma = Sigma_R1,cs = c(1,-1))
## R2 regime
Phi_R2 = list(phi1 = matrix(c(0.3,0.2,0.5,0.7),k,k))
Sigma_R2 = matrix(c(2.5,0.5,0.5,1),k,k,byrow = TRUE)
R2 = mtaregime(orders = list(p = 1,q = 0,d = 0),Phi = Phi_R2,Sigma = Sigma_R2,cs = c(5,2))
## create list of regime-type objects
Rg = list(R1 = R1,R2 = R2)
r = quantile(Ut$Sim$Yt[,1],0.4)
# Simulacion
datasim = mtarsim(N = Tlen,Rg = Rg,r = r,Zt = Ut$Sim$Yt[,1],Xt = Ut$Sim$Yt[,2],seed = 1997)
autoplot.tsregime(datasim$Sim,1)
autoplot.tsregime(datasim$Sim,2)
### Simular datos Faltantes -----------------------------------------------.####
pos_Nas_Y_1 = sample(1:1000,20)
pos_Nas_Y_2 = sample(1:1000,20)
pos_Nas_X = sample(1:1000,10)
pos_Nas_Z = sample(1:1000,5)

Yt_sim = datasim$Sim$Yt
Yt_sim[pos_Nas_Y_1,1] = NA
Yt_sim[pos_Nas_Y_2,2] = NA
Xt_sim = datasim$Sim$Xt
Xt_sim[pos_Nas_X,1] = NA
Zt_sim = datasim$Sim$Zt
Zt_sim[pos_Nas_Z,1] = NA

data_miss = tsregime(Yt = Yt_sim,Zt = Zt_sim,Xt = Xt_sim)
print(data_miss)                                                                #------------------------> Revisar Print de objetos
autoplot.tsregime(data_miss,1)
autoplot.tsregime(data_miss,2)                                                  #------------------------> Nombre del grafico
autoplot.tsregime(data_miss,3)
### Estimacion de parametros y datos faltantes ----------------------------.####
#### 1) Completamos con promedio
Y_temp = data_miss$Yt
Y_temp[,1][is.na(Y_temp[,1])] = mean(Y_temp[,1],na.rm = TRUE)
Y_temp[,2][is.na(Y_temp[,2])] = mean(Y_temp[,2],na.rm = TRUE)
X_temp = data_miss$Xt
X_temp[,1][is.na(X_temp[,1])] = mean(X_temp[,1],na.rm = TRUE)
Z_temp = data_miss$Zt
Z_temp[,1][is.na(Z_temp[,1])] = mean(Z_temp[,1],na.rm = TRUE)
### 2) Estimar numero de regimenes (Con datos completados por media)
data_temp = tsregime(Y_temp,Z_temp,X_temp)
initial = mtarinipars(tsregime_obj = data_temp,list_model = list(l0_max = 3),method = 'KUO')
estim_nr = mtarnumreg(ini_obj = initial,iterprev = 500,niter_m = 1000,burn_m = 500, list_m = TRUE,ordersprev = list(maxpj = 2,maxqj = 2,maxdj = 2),parallel = FALSE)
print(estim_nr)
### 3) Estimacion de parametros estructurales y no estructurales
initial = mtarinipars(tsregime_obj = data_temp,method = 'KUO',
                      list_model = list(pars = list(l = estim_nr$final_m),
                                        orders = list(pj = c(2,2))))
estruc = mtarstr(ini_obj = initial,niter = 500,chain = TRUE,parallel = FALSE)
autoplot.regime_model(estruc,1)
autoplot.regime_model(estruc,2)
autoplot.regime_model(estruc,3)
autoplot.regime_model(estruc,4)
autoplot.regime_model(estruc,5)
diagnostic_mtar(estruc)
### 4) Estimar los datos faltantes
list_model = list(pars = list(l = estim_nr$final_m,r = estruc$estimates$r[,2],orders = estruc$orders))
initial = mtarinipars(tsregime_obj = datasim_miss,list_model = list_model)
missingest = mtarmissing(ini_obj = initial,chain = TRUE, niter = 500,burn = 500)
print(missingest)
autoplot.regime_missing(missingest,1)
data_c = missingest$tsregim
### 5) Estimar nuevamente numero de regimenes y parametros
nitial = mtarinipars(tsregime_obj = data_c,list_model = list(l0_max = 3),method = 'KUO')
estim_nr = mtarnumreg(ini_obj = initial,iterprev = 500,niter_m = 500,burn_m = 500,
                      list_m = TRUE,ordersprev = list(maxpj = 2,maxqj = 2,maxdj = 2))
print(estim_nr)

initial = mtarinipars(tsregime_obj = data_c,method = 'KUO',
                      list_model = list(pars = list(l = estim_nr$final_m),orders = list(pj = c(2,2))))
estruc = mtarstr(ini_obj = initial,niter = 500,chain = TRUE)
autoplot.regime_model(estruc,1)
autoplot.regime_model(estruc,2)
autoplot.regime_model(estruc,3)
autoplot.regime_model(estruc,4)
autoplot.regime_model(estruc,5)
diagnostic_mtar(estruc)
### 6) Predicciones

