##### Ejemplos uso libreria dlm para Modelos espacio y estado filtro Kalman ######
library(dlm)
library(ks)
#Yt = F_t*Theta_t+v_t vt \sim N(0,V_t)
#Theta_t = G_t*Theta_{t-1}+w_t w_t \sim N(0,Wt)
# Inicial \Theta_0 \sim N(m0,C0)
#--Definir el modelo--#
##-- Supongamos que tenemos un modelo VAR(2)--#
Tlen = 1000
k = 2
nu = 0
Sigma_ut = expm::sqrtm(matrix(c(1,0.4,0.4,2),k,k))
Phi_ut = list(phi1 = matrix(c(0.2,0.3,0.3,0.2),k,k,byrow = T),
              phi2 = matrix(c(0.1,0.2,0.2,0.1),k,k,byrow = T))
R_ut = list(R1 = mtaregim(orders = list(p = 2,q = 0,d = 0),Phi = Phi_ut,Sigma = Sigma_ut))
Ut = mtarsim(N = Tlen,Rg = R_ut,seed = 123)$Yt
forecast::autoplot(ts(Ut),facets = TRUE)
#---Escribimos los coeficientes en modelo espacio estado--#
FF = cbind(diag(k),matrix(0,k,k))
V = matrix(0,k,k)
GG = rbind(cbind(Phi_ut$phi1,Phi_ut$phi2),cbind(diag(k),matrix(0,k,k)))
W = (rbind(Sigma_ut,matrix(0,k,k))) %*% diag(2) %*% t(rbind(Sigma_ut,matrix(0,k,k)))
m0 = matrix(0,2*k,1)
C0 = invvec(solve(diag(4*k*k) - GG %x% GG) %*% vec(W))
model = dlm(FF = FF, V = V, GG = GG, W = W, m0 = m0, C0 = C0)
#---Procesos con el filtro de Kalman---#
Filt = dlmFilter(ts(Ut),model)
plot(ts(Ut[,1]), type = 'l', col = "seagreen")
lines(dropFirst(Filt$m)[,2], type = 'o',pch = 20, col = "brown")

Smooth = dlmSmooth(Filt)
plot(ts(Ut[,1]), type = 'l', col = "seagreen")
lines(dropFirst(Smooth$s)[,2], type = 'o',pch = 20, col = "brown")

Forcast = dlmForecast(Filt,nAhead = 50)
plot(ts(Ut[,1]), type = 'l', col = "seagreen",xlim = c(0,1010))
lines(Forcast$f[,1])

#----MARSS----#
#---Simulamos un proceso 1 solo regimen---#
Tlen = 1000
k = 2
nu = 0
Sigma_ut = expm::sqrtm(matrix(c(1,0.4,0.4,2),k,k))
Phi_ut = list(phi1 = matrix(c(0.5,0.1,0.4,0.5),k,k,byrow = T))
R_ut = list(R1 = mtaregim(orders = list(p = 1,q = 0,d = 0),Phi = Phi_ut,Sigma = Sigma_ut))
Ut = mtarsim(N = Tlen,Rg = R_ut)$Yt
forecast::autoplot(ts(Ut),facets = TRUE)
## R1 regime
Phi_R1 = list(phi2 = matrix(c(0.1,0.6,-0.4,0.5),k,k,byrow = T))
Beta_R1 = list(beta1 = matrix(c(0.6,1),k,1))
Sigma_R1 = matrix(c(1,0,0,1),k,k,byrow = T)
cs_R1 = matrix(c(1,-1),nrow = 2)
R1 = mtaregim(orders = list(p = 2,q = 1,d = 0),Phi = Phi_R1,
              Beta = Beta_R1,Sigma = Sigma_R1,cs = cs_R1)
## get the simulation
Xt = Ut[,1]
datasim = mtarsim(N = Tlen,Rg = list(R1 = R1),Ut = Xt)
forecast::autoplot(ts(datasim$Yt),facets = T)

m = list()
C = list()

FF = cbind(diag(k),matrix(0,k,k + 1))
V = matrix(0,k,k)
GG = rbind(cbind(R1$phi$phi1,R1$phi$phi2,R1$beta$beta1),
           rbind(cbind(diag(k),matrix(0,k,k + 1)),matrix(0,1,2*k + 1)))
W = rbind(R1$sigma,matrix(0,k + 1,k)) %*% diag(k) %*% t(rbind(R1$sigma,matrix(0,k + 1,k)))
m0 = {rbind(R1$cs,matrix(0,3,1)) + rbind(matrix(0,4,1),1) %*% t(Xt[1])}
C0 = invvec(solve(diag((2*k + 1)*(2*k + 1)) - GG %x% GG) %*% vec(W))
model = dlm(FF = FF, V = V, GG = GG, W = W, m0 = m0, C0 = C0)
Filt = dlmFilter(ts(datasim$Yt[2:3,]),model,simplify = TRUE)
# Este numero es 1 es C0,2 es C1 y asi
Cn = Filt$U.C[[1]] %*% diag(Filt$D.C[1,]^2) %*% t(Filt$U.C[[1]])
mn = Filt$m[1,]
Smooth = dlmSmooth(Filt)
Cn2 = Smooth$s[1,]
mn2 = Smooth$U.S[[1]] %*% diag(Smooth$D.S[1,]^2) %*% t(Smooth$U.S[[1]])


#----Programacion de el filtro Kalman para optener los momentos---#
Tlen = 1000
k = 2
nu = 0
Sigma_ut = expm::sqrtm(matrix(c(1,0.4,0.4,2),k,k))
Phi_ut = list(phi1 = matrix(c(0.2,0.3,0.3,0.2),k,k,byrow = T),
              phi2 = matrix(c(0.1,0.2,0.2,0.1),k,k,byrow = T))
R_ut = list(R1 = mtaregim(orders = list(p = 2,q = 0,d = 0),Phi = Phi_ut,Sigma = Sigma_ut))
Ut = mtarsim(N = Tlen,Rg = R_ut,seed = 123)$Yt
forecast::autoplot(ts(Ut),facets = TRUE)
##---Crear las matrices---##
yt = t(Ut)
H = rbind(cbind(Phi_ut$phi1,Phi_ut$phi2),cbind(diag(k),matrix(0,k,k)))
K = cbind(diag(k),matrix(0,k,k))
L = matrix(0,4,1)
#M = matrix()
R = rbind(R_ut$R1$sigma,matrix(0,2,2))
##---Propagacion---##
PtC = AlphatC = vector("list",Tlen + 1)
# iniciales
ytnas = yt 
ytnas[1,2] = 0
AlphatC[[1]] = solve(diag(4) - H) %*% L
PtC[[1]] = 500*diag(4)
t = 2
# iteraciones
for (i in 2:{Tlen + 1}) {
  AlphatC[[i]] = H %*% AlphatC[[i - 1]] + L # %*% ut[i,]
  PtC[[i]] = H %*% PtC[[i - 1]] %*% t(H) + R %*% diag(k) %*% t(R)
}
# prediccion
QtC = ytC = vector("list",Tlen)
for (i in 1:{Tlen}) {
  ytC[[i]] = K %*% AlphatC[[i]]
  QtC[[i]] = K %*% PtC[[i]] %*% t(K)
}
# correccion
Pt = Alphat = vector("list",Tlen)
for (i in 1:{Tlen}) {
  St = PtC[[i]] %*% t(K) %*% solve(QtC[[i]])
  Alphat[[i]] = AlphatC[[i]] + St %*% {ytnas[,t] - ytC[[1]]}
  Pt[[i]] = (diag(2*k) - St %*% K) %*% PtC[[i]]
}
# Smoother Pt|T
PT = AlphaT = vector("list",Tlen)
AlphaT[[Tlen]] = Alphat[[Tlen]]
PT[[Tlen]] = Pt[[Tlen]]
for (i in rev(1:{Tlen})[-1]) {
  Ptx = Pt[[i]] %*% t(H) %*% solve(PtC[[i + 1]])
  AlphaT[[i]] = Alphat[[i]] + Ptx %*% {AlphaT[[i + 1]] - H %*% Alphat[[i]]} 
  PT[[i]] = Pt[[i]] + Ptx %*% {PT[[i + 1]] - PtC[[i + 1]]} %*% t(Ptx)
}
