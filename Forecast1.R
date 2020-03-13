# Estimacion de datos h pasaos adelante dado que conocemos Los parametros del modelo
# y la dencidad kernel para Ut.
yt = datasim;Ut = Ut;b = 1
## Kernel de Ut por ahora para modelo AR
modelU = mtarns(Yt = Ut,l = 1,orders = list(pj = 1,qj = 0,dj = 0),
                niter = 1000,chain = FALSE)$regime$R1
Ugen = function(t,U){
  cs = modelU$cs
  At = as.matrix(as.data.frame(modelU$phi))
  Sig = as.matrix(modelU$sigma)
  uti = c()
  for (w in 1:b) {uti = c(uti,U[,t - w])}
  val = MASS::mvrnorm(1, cs + At %*% uti,Sig %*% Sig)
  return(c(val))
}
## Corrida para Thetay
estimyt = mtarns(Yt = yt$Yt,Ut = Ut,l = 2,orders = list(pj = yt$pj,dj = yt$dj,qj = yt$qj)
                ,r = qnorm(0.4),niter = 1000,chain = TRUE)
#Creando objetos de entrada para generar Yt
Yt = t(estimyt$data$yt)
Ut = t(estimyt$data$Ut)
nu = nrow(Ut) - 1
l = length(estimyt$regime)
r = estimyt$r
k = nrow(Yt)
N = ncol(Yt)
pj = sapply(estimyt$regime,function(x){length(x$phi)})
qj = sapply(estimyt$regime,function(x){length(x$beta)})
dj = sapply(estimyt$regime,function(x){length(x$delta)})
eta = 1 + pj*k + qj*nu + dj
Zt = Ut[1,]
Xt = matrix(Ut[-1,],nrow = nu,ncol = N,byrow = TRUE)
# Vector de parametros
theta_iter = estimyt$Chain$Theta
sigma_iter = estimyt$Chain$Sigma
i2 = 1
# Generacion de Yt
Ygen = function(t,Y,U,...){
  rj = matrix(nrow = 2,ncol = l)
  if (l == 1) {
    rj[,1] = c(-Inf,Inf)
    }else{
      rj[,1] = c(-Inf,r[1])
      rj[,l] = c(rev(r)[1],Inf)
    }
  if (l > 2) {for (i2 in 2:{l - 1}) {rj[,i2] = c(r[i2 - 1],r[i2])}}
  # indicator variable for the regime
  for (j in 1:l) {
    if (U[1,t] > rj[1,j] & U[1,t] <= rj[2,j]) {Ind = j}
    }
  p = pj[Ind]
  q = qj[Ind]
  d = dj[Ind]
  maxj = max(p,q,d)
  # matrix wj =(1,lagY,lagX,lagZ)
  yti = c()
  for (w in 1:p) {yti = c(yti,Y[,t - w])}
  xti = c()
  for (w in 1:q) {xti = c(xti,U[-1,t - w])}
  zti = c()
  for (w in 1:d) {zti = c(zti,U[1,t - w])}
  if (q == 0 & d != 0) {
    wtj = c(1,yti,zti)
    }else if (d == 0 & q != 0) {
      wtj = c(1,yti,xti)
      }else if (d == 0 & q == 0) {
        wtj = c(1,yti)
      }else{
          wtj = c(1,yti,xti,zti)}
  mean = {t(wtj) %x% diag(k)} %*% theta_iter[[Ind]][,i2]
  cov = ks::invvec(estimyt$Chain$Sigma[[Ind]][,i2])
  val = MASS::mvrnorm(1,mean,cov %*% cov)
  return(val)
}
#Funcion solo con ns----------------------------------------------------------------------#
## Entradas de la funcion
yt = datasim
estim = mtarns(Yt = yt$Yt,Ut = Ut,l = 2,orders = list(pj = yt$pj,dj = yt$dj,qj = yt$qj)
               ,r = qnorm(0.4),niter = 1500,chain = TRUE)
niter = 1000
b = 1
h = 50
level = 0.95
## Estimacion de Ut kernel y funcion para generar
modelU = mtarns(Yt = Ut,l = 1,orders = list(pj = 1,qj = 0,dj = 0),
                niter = niter,chain = FALSE)$regime$R1
Ugen = function(t,U){
  cs = modelU$cs
  At = as.matrix(as.data.frame(modelU$phi))
  Sig = as.matrix(modelU$sigma)
  uti = c()
  for (w in 1:b) {uti = c(uti,U[,t - w])}
  val = MASS::mvrnorm(1, cs + At %*% uti,Sig %*% Sig)
  return(c(val))
}
## Funcion para generar Yt
Yt = t(estim$data$yt)
Ut = t(estim$data$Ut)
nu = nrow(Ut) - 1
l = length(estim$regime)
r = estim$r
k = nrow(Yt)
N = ncol(Yt)
pj = sapply(estim$regime,function(x){length(x$phi)})
qj = sapply(estim$regime,function(x){length(x$beta)})
dj = sapply(estim$regime,function(x){length(x$delta)})
eta = 1 + pj*k + qj*nu + dj
### Vector de parametros
theta_iter = estim$Chain$Theta
sigma_iter = estim$Chain$Sigma
### generar Yt
Ygen = function(t,Y,U,i2...){
  rj = matrix(nrow = 2,ncol = l)
  if (l == 1) {
    rj[,1] = c(-Inf,Inf)
  }else{
    rj[,1] = c(-Inf,r[1])
    rj[,l] = c(rev(r)[1],Inf)
  }
  if (l > 2) {for (i2 in 2:{l - 1}) {rj[,i2] = c(r[i2 - 1],r[i2])}}
  # indicator variable for the regime
  for (j in 1:l) {
    if (U[1,t] > rj[1,j] & U[1,t] <= rj[2,j]) {Ind = j}
  }
  p = pj[Ind]
  q = qj[Ind]
  d = dj[Ind]
  maxj = max(p,q,d)
  # matrix wj =(1,lagY,lagX,lagZ)
  yti = c()
  for (w in 1:p) {yti = c(yti,Y[,t - w])}
  xti = c()
  for (w in 1:q) {xti = c(xti,U[-1,t - w])}
  zti = c()
  for (w in 1:d) {zti = c(zti,U[1,t - w])}
  if (q == 0 & d != 0) {
    wtj = c(1,yti,zti)
  }else if (d == 0 & q != 0) {
    wtj = c(1,yti,xti)
  }else if (d == 0 & q == 0) {
    wtj = c(1,yti)
  }else{
    wtj = c(1,yti,xti,zti)}
  mean = {t(wtj) %x% diag(k)} %*% theta_iter[[Ind]][,i2]
  cov = ks::invvec(sigma_iter[[Ind]][,i2])
  val = MASS::mvrnorm(1,mean,cov %*% cov)
  return(val)
}
### Cadenas para las estimaciones
ChainYt = matrix(ncol = k*h,nrow = niter)
ChainUt = matrix(ncol = (nu + 1)*h,nrow = niter)
Yp = Yt;Up = Ut
Ti = 5
Upi = matrix(nrow = nrow(Ut),ncol = h)
Ypi = matrix(nrow = nrow(Yt),ncol = h)
### Proceso iterativo
pb = txtProgressBar(min = 1, max = niter, style = 3)
for (i2 in 1:niter) {
  for (i3 in 1:h) {
    Upi[,i3] = Ugen(Ti + i3,Up)
    Up = cbind(Up,Upi[,i3])
    Ypi[,i3] = Ygen(Ti + i3,Yp,Up,i2)
    Yp = cbind(Yp,Ypi[,i3])
  }
  ChainYt[i2,] = ks::vec(Ypi)
  ChainUt[i2,] = ks::vec(Upi)
  setTxtProgressBar(pb,i2)
}
### Calculo de estimaciones y intervalos
estimYt = matrix(nrow = k*h,ncol = 4)
estimUt = matrix(nrow = (nu + 1)*h,ncol = 4)
colnames(estimUt) = colnames(estimYt) = 
  c(paste('lower limit ',(1 - level)/2*100,'%',sep = ""),'mean',paste('upper limit ',(1 + level)/2*100,'%',sep = ""),"RVPD")
#rchain = matrix(r_iter[,-c(1:burn)],ncol = niter,nrow = l - 1)
estimUt[,1] = apply(ChainUt,2,quantile,probs = (1 - level)/2)
estimUt[,2] = apply(ChainUt,2,mean)
estimUt[,3] = apply(ChainUt,2,quantile,probs = (1 + level)/2)
estimUt[,4] = apply(ChainUt,2,sd)
estimYt[,1] = apply(ChainYt,2,quantile,probs = (1 - level)/2)
estimYt[,2] = apply(ChainYt,2,mean)
estimYt[,3] = apply(ChainYt,2,quantile,probs = (1 + level)/2)
estimUt[,4] = apply(ChainYt,2,sd)
row.names(estimYt) = paste(paste("t",(1:h) %x% rep(1,k) ,sep = "+"),rep(1,h) %x% (1:k),sep = ".")
row.names(estimUt) = paste(paste("t",(1:h) %x% rep(1,nu + 1) ,sep = "+"),rep(1,h) %x% (1:{nu + 1}),sep = ".")
Yth = ks::invvec(estimYt[,2],ncol = h,nrow = k)
Uth = ks::invvec(estimUt[,2],ncol = h,nrow = nu + 1)

forecast::autoplot(ts(ChainUt),facets = TRUE)
forecast::autoplot(ts(ChainYt),facets = TRUE)

forecast::autoplot(ts(t(Yth)),facets = TRUE)
forecast::autoplot(ts(t(Uth)),facets = TRUE)

plot(f2[1,],type = "l")
abline(f3[1,])

plot(f3[1,],type = "l")
lines(f2[1,],type = "l",col = "blue")

serielab = factor(rep(1,h) %x% 1:2)
time = ks::vec(matrix(1:h,k,h,T))
dataplot = data.frame(estim = estimYt[,2],Lw = estimYt[,3],Ul = estimYt[,1],time = time,lab = serielab)
dataplot = dataplot[1:995,]
dataplot = cbind(dataplot,t2 = ks::vec(Yt[,-c(1:5)]))
ggplot(data = dataplot,aes(x = time,y = estim)) + geom_line() +
  facet_grid(lab~.) + geom_ribbon(aes(ymin = Lw,ymax = Ul),fill = "gray",alpha = 0.5)
  
#Funcion con str----------------------------------------------------------------------------------####
## Entradas de la funcion
yt = datasim
estim = mtarstr(Yt = yt$Yt,Ut = Ut,l = 2,orders = list(pjmax = c(2,2),
        qjmax = c(1,1),djmax = c(1,1)),niter = 1000,method = 'KUO',chain = T)

niter = 1000
b = 1
h = 1000
level = 0.95
## Estimacion de Ut kernel y funcion para generar
modelU = mtarns(Yt = Ut,l = 1,orders = list(pj = 1,qj = 0,dj = 0),
                niter = niter,chain = FALSE)$regime$R1
Ugen = function(t,U){
  cs = modelU$cs
  At = as.matrix(as.data.frame(modelU$phi))
  Sig = as.matrix(modelU$sigma)
  uti = c()
  for (w in 1:b) {uti = c(uti,U[,t - w])}
  val = MASS::mvrnorm(1, cs + At %*% uti,Sig %*% Sig)
  return(c(val))
}
## Funcion para generar Yt
Yt = t(estim$data$yt)
Ut = t(estim$data$Ut)
nu = nrow(Ut) - 1
l = length(estim$regime)
r = estim$r
k = nrow(Yt)
N = ncol(Yt)
pj = c(2,2)
qj = c(1,1)
dj = c(1,1)
eta = 1 + pj*k + qj*nu + dj
### Vector de parametros
theta_iter = estim$Chain$Theta
sigma_iter = estim$Chain$Sigma
gam_iter = estim$Chain$Gamma
r_iter = estim$Chain$r
### generar Yt
Ygen = function(t,Y,U,i2...){
  r = r_iter[i2]
  rj = matrix(nrow = 2,ncol = l)
  if (l == 1) {
    rj[,1] = c(-Inf,Inf)
  }else{
    rj[,1] = c(-Inf,r[1])
    rj[,l] = c(rev(r)[1],Inf)
  }
  if (l > 2) {for (i2 in 2:{l - 1}) {rj[,i2] = c(r[i2 - 1],r[i2])}}
  # indicator variable for the regime
  for (j in 1:l) {
    if (U[1,t] > rj[1,j] & U[1,t] <= rj[2,j]) {Ind = j}
  }
  p = pj[Ind]
  q = qj[Ind]
  d = dj[Ind]
  maxj = max(p,q,d)
  # matrix wj =(1,lagY,lagX,lagZ)
  yti = c()
  for (w in 1:p) {yti = c(yti,Y[,t - w])}
  xti = c()
  for (w in 1:q) {xti = c(xti,U[-1,t - w])}
  zti = c()
  for (w in 1:d) {zti = c(zti,U[1,t - w])}
  if (q == 0 & d != 0) {
    wtj = c(1,yti,zti)
  }else if (d == 0 & q != 0) {
    wtj = c(1,yti,xti)
  }else if (d == 0 & q == 0) {
    wtj = c(1,yti)
  }else{
    wtj = c(1,yti,xti,zti)}
  mean = {t(wtj) %x% diag(k)} %*% (gam_iter[[Ind]][,i2]*theta_iter[[Ind]][,i2])
  cov = ks::invvec(sigma_iter[[Ind]][,i2])
  val = MASS::mvrnorm(1,mean,cov %*% cov)
  return(val)
}
### Cadenas para las estimaciones
ChainYt = matrix(ncol = k*h,nrow = niter)
ChainUt = matrix(ncol = (nu + 1)*h,nrow = niter)
Yp = Yt;Up = Ut
Ti = 5
Upi = matrix(nrow = nrow(Ut),ncol = h)
Ypi = matrix(nrow = nrow(Yt),ncol = h)
### Proceso iterativo
pb = txtProgressBar(min = 1, max = niter, style = 3)
for (i2 in 1:niter) {
  for (i3 in 1:h) {
    Upi[,i3] = Ugen(Ti + i3,Up)
    Up = cbind(Up,Upi[,i3])
    Ypi[,i3] = Ygen(Ti + i3,Yp,Up,i2)
    Yp = cbind(Yp,Ypi[,i3])
  }
  ChainYt[i2,] = ks::vec(Ypi)
  ChainUt[i2,] = ks::vec(Upi)
  setTxtProgressBar(pb,i2)
}
### Calculo de estimaciones y intervalos
estimYt = matrix(nrow = k*h,ncol = 4)
estimUt = matrix(nrow = (nu + 1)*h,ncol = 4)
colnames(estimUt) = colnames(estimYt) = 
  c(paste('lower limit ',(1 - level)/2*100,'%',sep = ""),'mean',paste('upper limit ',(1 + level)/2*100,'%',sep = ""),"RVPD")
#rchain = matrix(r_iter[,-c(1:burn)],ncol = niter,nrow = l - 1)
estimUt[,1] = apply(ChainUt,2,quantile,probs = (1 - level)/2)
estimUt[,2] = apply(ChainUt,2,mean)
estimUt[,3] = apply(ChainUt,2,quantile,probs = (1 + level)/2)
estimUt[,4] = apply(ChainUt,2,sd)
estimYt[,1] = apply(ChainYt,2,quantile,probs = (1 - level)/2)
estimYt[,2] = apply(ChainYt,2,mean)
estimYt[,3] = apply(ChainYt,2,quantile,probs = (1 + level)/2)
estimUt[,4] = apply(ChainYt,2,sd)
row.names(estimYt) = paste(paste("t",(1:h) %x% rep(1,k) ,sep = "+"),rep(1,h) %x% (1:k),sep = ".")
row.names(estimUt) = paste(paste("t",(1:h) %x% rep(1,nu + 1) ,sep = "+"),rep(1,h) %x% (1:{nu + 1}),sep = ".")
Yth = ks::invvec(estimYt[,2],ncol = h,nrow = k)
Uth = ks::invvec(estimUt[,2],ncol = h,nrow = nu + 1)

forecast::autoplot(ts(ChainUt),facets = TRUE)
forecast::autoplot(ts(ChainYt),facets = TRUE)

forecast::autoplot(ts(t(Yth)),facets = TRUE)
forecast::autoplot(ts(t(Uth)),facets = TRUE)

serielab = factor(rep(1,h) %x% 1:2)
time = ks::vec(matrix(1:h,k,h,T))
dataplot = data.frame(estim = estimYt[,2],Lw = estimYt[,3],Ul = estimYt[,1],time = time,lab = serielab)

ggplot(data = dataplot,aes(x = time,y = estim)) + geom_line() +
  facet_grid(lab~.) + geom_ribbon(aes(ymin = Lw,ymax = Ul),fill = "gray",alpha = 0.5)