Yzz = yt
# Simulacion de datos faltantes
forecast::autoplot(ts(yt$Yt),facets = TRUE)
posNA = sample(c(1:1000),4)
yt$Yt[c(posNA),] = c(NA,NA)
posNA = sample(c(1:1000),4)
Ut[c(posNA),] = c(NA,NA)
forecast::autoplot(ts(yt$Yt),facets = TRUE)
forecast::autoplot(ts(Ut),facets = TRUE)
# arguments of function mtarmiss
Yt = yt$Yt
Ut = Ut
pj = yt$pj
qj = yt$qj
dj = yt$dj
r = -0.22
l = 2
niter = 500
chain = T
level = 0.95
burn = 100
cU = 0.5
thetaini = sigmaini = b = NULL
# first entries
Yt = t(Yt)
Ut = t(Ut)
k = nrow(Yt)
N = ncol(Yt)
nu = nrow(Ut) - 1
b = ifelse(is.null(b),1,b)
PosNAMat = PosNAvec = PosNAvecT = vector(mode = 'list',2)
PosNAMat[[1]] = apply(Yt,2,is.na)
PosNAvec[[1]] = c(1:ncol(Yt))[apply(PosNAMat[[1]],2,any)]
PosNAvecT[[1]] = matrix(rep(c(1:Tlen),k),nrow = k,ncol = Tlen,byrow = T)[PosNAMat[[1]]]

PosNAMat[[2]] = apply(Ut,2,is.na)
PosNAvec[[2]] = c(1:ncol(Ut))[apply(PosNAMat[[2]],2,any)]
PosNAvecT[[2]] = matrix(rep(c(1:Tlen),nu + 1),nrow = nu + 1,ncol = Tlen,byrow = T)[PosNAMat[[2]]]

Zt = Ut[1,]
if (nu == 0) {
  Xt = matrix(0,ncol = N,nrow = 1)
  qjmax = rep(0,l)
}else{
  Xt = matrix(Ut[-1,],nrow = nu,ncol = N,byrow = TRUE)
}
etaj = 1 + k*pj + nu*qj + dj
#Completamos Ut faltantes con promedios Ut
if (length(PosNAvec[[2]]) != 0) {
  meanU = apply(Ut,1,mean,na.rm = TRUE)
  for (i in 1:nrow(Ut)) {
    Ut[i,PosNAMat[[2]][i,]] = meanU[i]
  }
}
#Completar datos faltantes con o en Yt y permutar en Yt y Kt(OJO)
modelU = mtarns(t(Ut),l = 1,orders = list(pj = b,qj = 0,dj = 0),niter = 2000,chain = FALSE,burn = 1000)
modelU = modelU$regime$R1
#functions
prodB = function(x){
  prod = 1
  for (a in 1:length(x)) {
    prod = prod*x[a]
  }
  return(prod)
}
dmnormB = function(x, mean, sigma){
  dist = Brobdingnag::as.brob(c(t(x - mean) %*% solve(sigma) %*% (x - mean)))
  cte = (2*pi)^{-nrow(sigma)/2}*determinant(sigma, logarithm = FALSE)$modulus^{-1/2}
  return(cte*exp(-1/2*dist))
}
lists = function(r, Yt, Ut,...){
  Zt = Ut[1,]
  if (nu == 0) {
    Xt = matrix(0,ncol = N,nrow = 1)
  }else{
    Xt = matrix(Ut[-1,],nrow = nu,ncol = N,byrow = TRUE)
  }
  rj = matrix(nrow = 2,ncol = l)
  if (l == 1) {
    rj[,1] = c(-Inf,Inf)
  }else{
    rj[,1] = c(-Inf,r[1])
    rj[,l] = c(rev(r)[1],Inf)
  }
  if (l > 2) {for (i2 in 2:{l - 1}) {rj[,i2] = c(r[i2 - 1],r[i2])}}
  # indicator variable for the regime
  Ind = c()
  for (j in 1:l) {
    for (w in 1:N) {
      if (Zt[w] > rj[1,j] & Zt[w] <= rj[2,j]) {
        Ind[w] = j
      }
    }
  }
  Nrg = c()
  listaWj = listaYj = list()
  length(listaWj) = length(listaYj) = l
  for (lj in 1:l) {
    p = pj[lj]
    q = qj[lj]
    d = dj[lj]
    maxj = max(p,q,d)
    Inj = which(Ind == lj)
    Inj = Inj[Inj > maxj]
    Nrg[lj] = length(Inj)
    Yj = matrix(Yt[,Inj],nrow = k,ncol = Nrg[lj])
    Wj = matrix(0,nrow = etaj[lj],ncol = Nrg[lj])
    count = 1
    for (ti in Inj) {
      yti = c()
      for (w in 1:p) {yti = c(yti,Yt[,ti - w])}
      xti = c()
      for (w in 1:q) {xti = c(xti,Xt[,ti - w])}
      zti = c()
      for (w in 1:d) {zti = c(zti,Zt[ti - w])}
      if (q == 0 & d != 0) {
        wtj = c(1,yti,zti)
      }else if (d == 0 & q != 0) {
        wtj = c(1,yti,xti)
      }else if (d == 0 & q == 0) {
        wtj = c(1,yti)
      }else{
        wtj = c(1,yti,xti,zti)}
      Wj[,count] = wtj
      count = count + 1
    }
    listaWj[[lj]] = Wj
    listaYj[[lj]] = Yj
  }
  return(list(Nrg = Nrg,listaW = listaWj,listaY = listaYj,Ind = Ind))
}
ker = function(t, Ut, ...){
  cs = modelU$cs
  At = as.matrix(as.data.frame(modelU$phi))
  Sig = as.matrix(modelU$sigma)
  EU = solve(diag(nu + 1) - At) %*% cs 
  vecVU = solve(diag(2*(nu + 1)) - At %x% At) %*% c(Sig %*% Sig)
  VU = ks::invvec(vecVU,ncol = nu + 1, nrow = nu + 1)
  val = dmnormB(Ut[,t], EU, VU)
  return(c(val))
}
transker = function(t, Ut, ...){
  p = length(modelU$phi)
  ## create matrix
  cs = modelU$cs
  At = as.matrix(as.data.frame(modelU$phi))
  Sig = as.matrix(modelU$sigma)
  ## make lags and calculate
  uti = c()
  for (w in 1:p) {uti = c(uti,Ut[,t - w])}
  val = dmnormB(Ut[,t], cs + At %*% uti, Sig %*% Sig)
  return(c(val))
}
kernU = Vectorize(ker,vectorize.args = 't')
transkernU = Vectorize(transker,vectorize.args = 't')
state_space = function(reg, iSS, theta, sigma, ...) {
  p = pj[reg]
  q = qj[reg]
  d = dj[reg]
  Aj = ks::invvec(theta[[reg]][,iSS],nrow = k, ncol = etaj[reg])
  if (p == 1) {
    Aj = cbind(Aj,Aj[,-1]*0,matrix(0,nrow = k,ncol = nu + 1))
  }
  paux = ifelse(p == 1,2,p)
  qaux = ifelse(q == 0,1,q)
  daux = ifelse(d == 0,1,d) 
  R_zt = t(cbind(expm::sqrtm(sigma[[reg]][[iSS]]),matrix(0,k,(paux - 1)*k + nu*qaux + daux)))
  L_zt = c(Aj[,1],rep(0,(paux - 1)*k + nu*qaux  + daux))
  hphi = cbind(diag(k*(paux - 1)),matrix(0,nrow = k*(paux - 1), ncol = k + qaux*nu + daux))
  hbeta = cbind(matrix(0,nrow = nu*(qaux - 1), ncol = k*paux),diag(nu*(qaux - 1)),matrix(0,nrow = nu*(qaux - 1),ncol = qaux*nu + daux))
  hdelta = cbind(matrix(0,nrow = daux - 1, ncol = k*paux + qaux*nu),diag(daux - 1),matrix(0,ncol = 1,nrow = daux - 1))
  H_zt = rbind(Aj[,-1], 
               hphi,
               matrix(0,nrow = nu, ncol = ncol(Aj) - 1),
               hbeta,
               matrix(0,nrow = 1, ncol = ncol(Aj) - 1),
               hdelta)
  K_zt = cbind(diag(k),matrix(0,k,(paux - 1)*k + nu*qaux + daux))
  M_zt = rbind(matrix(0,nrow = k*paux, ncol = nu + 1),
                     cbind(matrix(0,nrow = nu, ncol = 1),diag(nu)),
                     matrix(0,nrow = nu*(qaux - 1), ncol = nu + 1),
                     c(1,rep(0,nu)),
                     matrix(0,nrow = (daux - 1), ncol = nu + 1))
  return(list(K = K_zt, L = L_zt, H = H_zt, M = M_zt, R = R_zt))
}
alphacond = function(t, iA, Ut, Yt, theta, sigma, ...) {
  Zt = Ut[1,]
  if (nu == 0) {
    Xt = matrix(0,ncol = N,nrow = 1)
  }else{
    Xt = matrix(Ut[-1,],nrow = nu,ncol = N,byrow = TRUE)
  }
  rj = matrix(nrow = 2,ncol = l)
  if (l == 1) {
    rj[,1] = c(-Inf,Inf)
  }else{
    rj[,1] = c(-Inf,r[1])
    rj[,l] = c(rev(r)[1],Inf)
  }
  if (l > 2) {for (i2 in 2:{l - 1}) {rj[,i2] = c(r[i2 - 1],r[i2])}}
  # indicator variable for the regime
  Ind = c()
  for (j in 1:l) {
    for (w in 1:N) {
      if (Zt[w] > rj[1,j] & Zt[w] <= rj[2,j]) {
        Ind[w] = j
      }
    }
  }
  lj = Ind[t]
  p = pj[lj]
  q = qj[lj]
  d = dj[lj]
  Wj = matrix(0,nrow = etaj[lj],ncol = 1)
  yti = c()
  for (w in 1:p) {yti = c(yti,Yt[,t - w])}
  xti = c()
  for (w in 1:q) {xti = c(xti,Xt[,t - w])}
  zti = c()
  for (w in 1:d) {zti = c(zti,Zt[t - w])}
  if (q == 0 & d != 0) {
    wtj = c(1,yti,zti)
  }else if (d == 0 & q != 0) {
    wtj = c(1,yti,xti)
  }else if (d == 0 & q == 0) {
    wtj = c(1,yti)
  }else{
    wtj = c(1,yti,xti,zti)}
  Wj[,1] = wtj
  Hj = ks::invvec(theta[[lj]][,iA],nrow = k,ncol = etaj[lj])
  val = dmnormB(Yt[,t], {Hj %*% Wj}, sigma[[lj]][[iA]])
  return(val)
}
alphacond = Vectorize(alphacond,vectorize.args = 't')
#objects for each regimen and iterations
theta_iter = sigma_iter = list()
length(theta_iter) = length(sigma_iter) = l
itheta0j = isigma0j = list()
length(itheta0j) = length(isigma0j) = l
iS0j = inu0j = list()
length(iS0j) = length(inu0j) = l

Yt_iter = matrix(ncol = niter + burn,nrow = sum(ks::vec(PosNAMat[[1]])))
Ut_iter = matrix(ncol = niter + burn,nrow = sum(ks::vec(PosNAMat[[2]])))
Ytr = Yt #Yt que vamos a cambiar en el proceso
Utr = Ut #Yt que vamos a cambiar en el proceso
#Ytr[PosNAMat] = 0
Yt_iter[,1] = Ytr[PosNAMat[[1]]]
Ut_iter[,1] = Utr[PosNAMat[[2]]] 
#set initial values for each regime in each chain
#creacion de cadenas para sigma y theta
for (lj in 1:l) {
  theta_iter[[lj]] = matrix(ncol = niter + burn,nrow = k*etaj[lj])
  if (is.null(thetaini)) {
    itheta0j[[lj]] = rep(0,k*etaj[lj])
    isigma0j[[lj]] = diag(k*etaj[lj])
  }else{
    itheta0j[[lj]] = thetaini[[lj]][[1]]
    isigma0j[[lj]] = thetaini[[lj]][[2]]
  }
  theta_iter[[lj]][,1] = mvtnorm::rmvnorm(1,mean = itheta0j[[lj]],sigma = isigma0j[[lj]])
  sigma_iter[[lj]] = list()
  length(sigma_iter[[lj]]) = niter + burn
  if (is.null(sigmaini)) {
    iS0j[[lj]] = diag(k)
    inu0j[[lj]] = k
  }else{
    iS0j[[lj]] = sigmaini[[lj]][[1]]
    inu0j[[lj]] = sigmaini[[lj]][[2]]
  }
  sigma_iter[[lj]][[1]] = LaplacesDemon::rinvwishart(nu = inu0j[[lj]],S = iS0j[[lj]])
}
#state-space model
K_zti = K_zt = list()
R_zt = list()
L_zt = list()
H_zt = list()
M_zt = list()
#Primera permutaciones
for (lj in 1:l) {
  listmatrix = state_space(lj, 1, theta_iter,sigma_iter)
  R_zt[[lj]] = listmatrix$R
  L_zt[[lj]] = listmatrix$L  
  H_zt[[lj]] = listmatrix$H
  K_zt[[lj]] = listmatrix$K
  M_zt[[lj]] = listmatrix$M  
}
listj = lists(r, Ytr, Ut)
for (ij in 1:Tlen) {
  K_zti[[ij]] = K_zt[[listj$Ind[ij]]]
}
#permutaciones para cadena Yt:
for (ij in PosNAvec[[1]]) {
  posNAi = PosNAMat[[1]][,ij]
  if (!all(posNAi)) {
    K_zti[[ij]] = K_zti[[ij]][order(posNAi),]
    Ytr[,ij] = Ytr[,ij][order(posNAi)]
  }
  K_zti[[ij]][is.na(Ytr[,ij]),] = K_zti[[ij]][is.na(Ytr[,ij]),]*0
  Ytr[,ij][is.na(Ytr[,ij])] = 0
}
sersalY = Ytr
# Sampling
pb = txtProgressBar(min = 2, max = niter + burn, style = 3)
acepU = 0
for (i in 2:{niter + burn}) {
  #State space model
  PtC = AlphatC = vector('list',Tlen)
  QtC = ytC = vector('list',Tlen)
  Pt = Alphat = vector('list',Tlen + 1)
  Alphat[[1]] = matrix(0,nrow = k*max(pj) + nu*max(qj) + max(dj),ncol = 1)
  Pt[[1]] = 10*diag(k*max(pj) + nu*max(qj) + max(dj))
  # iteraciones
  Indi = listj$Ind
  for (i1 in 1:{Tlen}) {
    #Prediction Equations mt|t-1
    AlphatC[[i1]] = H_zt[[Indi[i1]]] %*% Alphat[[i1]] + L_zt[[Indi[i1]]] + M_zt[[Indi[i1]]] %*% Ut[,i1]
    R2 = R_zt[[Indi[i1]]] %*% diag(k) %*% t(R_zt[[Indi[i1]]])
    PtC[[i1]] = H_zt[[Indi[i1]]] %*% Pt[[i1]] %*% t(H_zt[[Indi[i1]]]) + R2
    ytC[[i1]] = K_zti[[i1]] %*% AlphatC[[i1]]
    QtC[[i1]] = K_zti[[i1]] %*% PtC[[i1]] %*% t(K_zti[[i1]])
    #Updating Equations mt
    St = PtC[[i1]] %*% t(K_zti[[i1]]) %*% MASS::ginv(QtC[[i1]])
    Alphat[[i1 + 1]] = AlphatC[[i1]] + St %*% {sersalY[,i1] - ytC[[i1]]}
    Pt[[i1 + 1]] = PtC[[i1]] - St %*% K_zti[[i1]] %*% PtC[[i1]]
  }
  #sampling for state vector (pg37)
  PT = AlphaT = vector('list',Tlen + 1)
  AlphaT[[Tlen + 1]] = Alphat[[Tlen + 1]]
  PT[[Tlen + 1]] = Pt[[Tlen + 1]]
  
  for (i1 in rev(1:{Tlen})) {
    Eig = eigen(Pt[[i1 + 1]])$values
    Eig = any(Mod(Eig) > exp(-6))
    if (Eig) {
      estUp = MASS::mvrnorm(1,AlphaT[[i1 + 1]],PT[[i1 + 1]])
    }else{
      estUp = AlphaT[[i1 + 1]]
    }
    R2 = R_zt[[Indi[i1]]][1:k,] %*% diag(k) %*% t(R_zt[[Indi[i1]]][1:k,])
    Qt = MASS::ginv(H_zt[[Indi[i1]]][1:k,] %*%  Pt[[i1]] %*% t(H_zt[[Indi[i1]]][1:k,]) + R2)
    Bt = Pt[[i1]] %*% t(H_zt[[Indi[i1]]][1:k,]) %*% Qt
    Gt = estUp[1:k] - M_zt[[Indi[i1]]][1:k,] %*% Ut[,i1] - L_zt[[Indi[i1]]][1:k] - H_zt[[Indi[i1]]][1:k,] %*% Alphat[[i1]]
    AlphaT[[i1]] = Alphat[[i1]] + Bt %*% Gt
    PT[[i1]] = Pt[[i1]] - Bt %*% H_zt[[Indi[i1]]][1:k,] %*% Pt[[i1]]
  }
  # Simulaion de datos faltantes en Yt
  for (i1 in PosNAvec[[1]]) {
    Ysim = as.matrix(MASS::mvrnorm(1,AlphaT[[i1 + 1]],PT[[i1 + 1]]))
    Yt_iter[,i][PosNAvecT[[1]] == i1] = Ysim[1:k][PosNAMat[[1]][,i1]]
    Ytr[,i1] = Ysim[1:k]
    AlphaT[[i1 + 1]] = Ysim
  }
  #random walk U
  for (i1 in PosNAvec[[2]]) {
    ek = mvtnorm::rmvnorm(1,mean = rep(0,nu + 1), sigma = cU*diag(nu + 1))
    # Simulacion de la propues
    Usim = Utr
    Usim[,i1] = Utr[,i1] + c(ek)
    Usim[!PosNAMat[[2]][,i1],i1] =  Utr[!PosNAMat[[2]][,i1],i1]
    # Calculo de las probabilidades segun sea el caso
    # Numerador
      if (i1 <= b) {
        prod1N = Reduce('*',kernU(1:b,Usim))
        prod2N = Reduce('*',alphacond(1:b,i - 1,Usim,Ytr,theta_iter,sigma_iter))
        prod3N = Reduce('*',transkernU({b + 1}:{2*b},Usim))
        prod1D = Reduce('*',kernU(1:b,Utr))
        prod2D = Reduce('*',alphacond(1:b,i - 1,Utr,Ytr,theta_iter,sigma_iter))
        prod3D = Reduce('*',transkernU({b + 1}:{2*b},Utr))
        val = (prod1N*prod2N*prod3N)/(prod1D*prod2D*prod3N)
        }else{
          prod1N = alphacond(i1,i - 1,Usim,Ytr,theta_iter,sigma_iter)[[1]]
          prod2N = Reduce('*',transkernU(i1:{i1 + b},Usim))
          prod1D = alphacond(i1,i - 1,Utr,Ytr,theta_iter,sigma_iter)[[1]]
          prod2D = Reduce('*',transkernU(i1:{i1 + b},Utr))
          val = (prod1N*prod2N)/(prod1D*prod2D)
          }
    if (val >= runif(1)) {
      Utr = Usim
      Ut_iter[,i][PosNAvecT[[2]] == i1] = Usim[PosNAMat[[2]][,i1],i1]
      }else{
        Utr = Utr
        Ut_iter[,i][PosNAvecT[[2]] == i1] = Utr[,i1][PosNAMat[[2]][,i1]]
      }
    }
  listj = lists(r, Ytr, Utr)
  for (lj in 1:l) {
    Wj = listj$listaW[[lj]]
    Yj = listj$listaY[[lj]]
    Nj = listj$Nrg[lj]
    yj = c(Yj)
    theta0j = itheta0j[[lj]]
    sigma0j = isigma0j[[lj]]
    S0j = iS0j[[lj]]
    nu0j = inu0j[[lj]]
    Vj = solve(Wj %*% t(Wj) %x% solve(sigma_iter[[lj]][[i - 1]]) + solve(sigma0j))
    thetaj = Vj %*% {(Wj %x% solve(sigma_iter[[lj]][[i - 1]])) %*% yj + solve(sigma0j) %*% theta0j}
    theta_iter[[lj]][,i] = mvtnorm::rmvnorm(1,mean = thetaj,sigma = Vj)
    Hj = ks::invvec(theta_iter[[lj]][,i],nrow = k,ncol = etaj[lj])
    Sj = (Yj - Hj %*% Wj) %*% t(Yj - Hj %*% Wj)
    sigma_iter[[lj]][[i]] = LaplacesDemon::rinvwishart(nu = Nj + nu0j,S = Sj + S0j)
  }
  #Actualizacion de las matrices espacio y estado
  for (lj in 1:l) {
    listmatrix = state_space(lj, i, theta_iter,sigma_iter)
    R_zt[[lj]] = listmatrix$R
    L_zt[[lj]] = listmatrix$L
    H_zt[[lj]] = listmatrix$H
    K_zt[[lj]] = listmatrix$K
    M_zt[[lj]] = listmatrix$M  
  }
  setTxtProgressBar(pb,i)
}
close(pb)