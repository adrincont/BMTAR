# arguments of function mtarnumreg
Yt = yt$Yt;Ut = Ut;l0 = 3
niter = 200;chain = T
level = 0.95
burn = 10
method = 'KUO'
kappa = 0.5
# first entries
Yt = t(Yt)
Ut = t(Ut)
k = nrow(Yt)
N = ncol(Yt)
nu = nrow(Ut) - 1
Zt = Ut[1,]
if (nu == 0) {
  Xt = matrix(0,ncol = N,nrow = 1)
}else{
  Xt = matrix(Ut[-1,],nrow = nu,ncol = N,byrow = TRUE)
}
a = quantile(Zt,probs = 0.2)
b = quantile(Zt,probs = 0.8)
#-------------------------------#
# independent functions #
prodB = function(x){
  prod = 1
  for (a in 1:length(x)) {
    prod = prod*x[a]
  }
  return(prod)
}
dmnormB = function(x, mean, sigma){
  dist = Brobdingnag::as.brob(t(x - mean) %*% solve(sigma) %*% (x - mean))
  cte = (2*pi)^{-nrow(sigma)/2}*determinant(sigma, logarithm = FALSE)$modulus^{-1/2}
  return(cte*exp(-1/2*dist))
}
dinvwishartB = function(x, nu, S){
  k = nrow(x)
  gamsum = 0
  for (i in 1:k) {
    gamsum = gamsum + lgamma((nu + 1 - i)/2)
  }
  dens = -(nu * k/2) * log(2) - ((k * (k - 1))/4) * log(pi) - 
    gamsum + (nu/2) * log(det(S)) - ((nu + k + 1)/2) * log(det(x)) - 
    0.5 * sum(diag(S %*% solve(x)))
  dens = exp(Brobdingnag::as.brob(dens))
  return(dens)
}
dwishartB = function(x, nu, S){
  k = nrow(x)
  gamsum = 0
  for (i in 1:k) {
    gamsum = gamsum + lgamma((nu + 1 - i)/2)
  }
  dens = -(nu * k/2) * log(2) - ((k * (k - 1))/4) * log(pi) - 
    gamsum - (nu/2) * log(det(S)) + ((nu - k - 1)/2) * log(det(x)) - 
    0.5 * sum(diag(solve(S) %*% x))
  dens = exp(Brobdingnag::as.brob(dens))
  return(dens)
}
rdunif = function(m,l0){
  sec = 2:l0
  sec = sec[sec != m]
  if (length(sec) == 1) {
    muestra = sec
  }else{
    muestra = sample(x = sec, size = 1)
  }
  return(muestra)
}
dmunif = function(l, r, a, b){
  l = length(r) + 1
  names(a) = names(b) = NULL
  volume = ((b - a)^{l - 1})/(factorial(l - 1))
  for (i in 1:{l - 1}) {
    if (r[i] >= a & r[i] <= b) {
      if (l <= 2) {prob = 1}
      if (l > 2) {
        prob = 1
        for (j in 1:{l - 2}) {
          if (r[j] < r[j + 1]) {prob = prob*1}else{prob = prob*0}
        }
      }
    }else{prob = 0}
  }
  rj = matrix(nrow = 2,ncol = l)
  rj[,1] = c(-Inf,r[1])
  rj[,l] = c(rev(r)[1],Inf)
  if (l > 2) {
    for (i2 in 2:{l - 1}) {rj[,i2] = c(r[i2 - 1],r[i2])}
  }
  Ind = c()
  for (j in 1:l) {
    for (w in 1:N) {
      if (Zt[w] > rj[1,j] & Zt[w] <= rj[2,j]) {
        Ind[w] = j
      }
    }
  }
  Nrg = c()
  for (lj in 1:l) {
    Nrg[lj] = length(Ind[Ind == lj])
  }
  if (sum(Nrg/sum(Nrg) > 0.2) == l) {prob = 1*prob}else{prob = 0*prob}
  return(prob/volume)
}
fycond = function(l, ir, listar, gamma, theta, sigma){
  acum = 0
  Nrg = listar$Nrg
  for (lj in 1:l) {
    yj = c(listar$listaY[[lj]])
    Xj = listar$listaX[[lj]]
    acum = acum + t(yj - Xj %*% diag(gamma[[lj]][,ir]) %*% theta[[lj]][,ir]) %*% 
      {diag(Nrg[lj]) %x% solve(sigma[[lj]][[ir]])} %*% 
      (yj - Xj %*% diag(gamma[[lj]][,ir]) %*% theta[[lj]][,ir])
  }
  sigmareg = lapply(sigma,function(x){x[[ir]]})
  cte = prodB(Brobdingnag::as.brob(sapply(sigmareg,function(x){
    return(c(determinant(x,logarithm = FALSE)$modulus))}))^{-Nrg/2})
  val = cte*exp(-1/2*Brobdingnag::as.brob(acum))
  return(val)
}
lists = function(l, r, pjmax, qjmax, djmax, ...){
  rj = matrix(nrow = 2,ncol = l)
  rj[,1] = c(-Inf,r[1])
  rj[,l] = c(rev(r)[1],Inf)
  if (l > 2) {for (j2 in 2:{l - 1}) {rj[,j2] = c(r[j2 - 1],r[j2])}}
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
  listaXj = listaYj = list()
  length(listaXj) = length(listaYj) = l
  etam = 1 + pjmax*k + qjmax*nu + djmax
  for (lj in 1:l) {
    p = pjmax[lj]
    q = qjmax[lj]
    d = djmax[lj]
    maxj = max(p,q,d)
    Inj = which(Ind == lj)
    Inj = Inj[Inj > maxj]
    Nrg[lj] = length(Inj)
    Yj = matrix(Yt[,Inj],nrow = k,ncol = Nrg[lj])
    Xj = matrix(1,nrow = k*Nrg[lj],ncol = k*etam[lj])
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
      Xj[count:{count + {k - 1}},] = diag(k) %x% t(wtj)
      count = count + k
    }
    listaXj[[lj]] = Xj
    listaYj[[lj]] = Yj
  }
  return(list(Nrg = Nrg,listaX = listaXj,listaY = listaYj))
}
# functions for each m #
fill = function(m, iter = 500, kappa = 0.5, ...){
  i = 1
  ordersm = list(pjmax = rep(3,m),qjmax = rep(2,m),djmax = rep(2,m))
  etam = 1 + ordersm$pjmax*k + ordersm$qjmax*nu + ordersm$djmax
  #Parametros priori Theta
  theta0Pm = lapply(1:m,function(x){rep(0,k*etam[x])})
  sigma0Pm = lapply(1:m,function(x){diag(k*etam[x])})
  #Parametros priori Sigma
  S0Pm = lapply(1:m,function(x){diag(k)})
  nu0Pm = lapply(1:m,function(x){5})
  #Parametros priori Gamma
  pij0Pm = lapply(1:m,function(x){rep(0.5,k*etam[x])})
  # 500 iterations
  par = mtarstr(Yt = t(Yt), Ut = t(Ut),l = m, orders = ordersm, niter = iter,rini = list(kappa = kappa), chain = TRUE)
  #Parametros seudo Theta
  theta0Sm = lapply(par$Theta,function(x){x[,2]})
  sigma0Sm = lapply(par$Chain$Theta,function(x){cov(t(x))})
  #Parametros seudo Sigma
  S0Sm = lapply(par$Regim,function(x){1/1000*x$sigma})
  nu0Sm = lapply(1:m,function(x){1000})
  #Parametros seudo Gamma
  pij0Sm = lapply(par$Chain$Gamma,function(x){apply(x,1,mean)})
  #Parametros seudo R
  rmean0Sm = par$r[,2]
  rcov0Sm = cov(t(par$Chain$r))
  # cadenas y primeros valores
  theta_iter = sigma_iter = gam_iter = Dj = Rj = list()
  length(theta_iter) = length(sigma_iter) = length(Rj) =
    length(gam_iter) = length(Dj) = m
  if (method == 'SSVS') {
    tauij = itauij = cij = list()
    length(tauij) = length(itauij) = length(cij) = m
  }
  r_iter = matrix(ncol = niter + burn,nrow = m - 1)
  for (lj in 1:m) {
    theta_iter[[lj]] = gam_iter[[lj]] = matrix(ncol = niter + burn,nrow = k*etam[lj])
    sigma_iter[[lj]] = list()
    length(sigma_iter[[lj]]) = niter + burn
    gam_iter[[lj]][,1] = rep(1,k*etam[lj])
    if (method == 'KUO') {
      theta_iter[[lj]][,1] = mvtnorm::rmvnorm(1,mean = theta0Pm[[lj]],sigma = sigma0Pm[[lj]])
    }
    if (method == 'SSVS') {
      cij[[lj]] = rep(25,k*etam[lj])
      tauij[[lj]] = ifelse(m == 2,rep(1.25,k*etam[lj]),rep(1.5,k*etam[lj]))
      itauij[[lj]] = cij[[lj]]*tauij[[lj]]
      Dj[[lj]] = diag(itauij[[lj]])
      Rj[[lj]] = diag(k*etam[lj])
      theta_iter[[lj]][,1] = mvtnorm::rmvnorm(1,mean = rep(0,k*etam[lj]),sigma = Dj[[lj]] %*% Rj[[lj]] %*% Dj[[lj]])
    }
    sigma_iter[[lj]][[1]] = LaplacesDemon::rinvwishart(nu = nu0Pm[[lj]],S = S0Pm[[lj]])
  }
  r_iter[,1] = c(quantile(Zt, probs = 1/m*(1:{m - 1})))
  # LISTS
  if (method == 'SSVS') {
    iniP = list(Theta = list(mean = theta0Pm, cov = sigma0Pm), Sigma = list(cov = S0Pm,gl = nu0Pm),
                Gamma = list(prob = pij0Pm, itauij = itauij, tauij = tauij, cij = cij, Rj = Rj))
  }else{
    iniP = list(Theta = list(mean = theta0Pm, cov = sigma0Pm), Sigma = list(cov = S0Pm,gl = nu0Pm),
                Gamma = list(prob = pij0Pm))
  }
  iniS = list(Theta = list(mean = theta0Sm,cov = sigma0Sm), Sigma = list(cov = S0Sm,gl = nu0Sm),
              Gamma = list(prob = pij0Sm), r = list(mean = rmean0Sm, cov = rcov0Sm))
  listchain = list(Theta = theta_iter, Sigma = sigma_iter,
                   Gamma = gam_iter, r = r_iter)
  listr = lists(m,r_iter[,1],ordersm$pjmax,ordersm$qjmax,ordersm$djmax)
  return(list(i = i,orders = ordersm,Priori = iniP,Pseudo = iniS,Chain = listchain,listr = listr))
}
updatelist = function(l, ...){
  rgamber = function(pos, reg, ig, ...){
    gam_j = gam_iter
    gam_j[[reg]][pos,ig] = 1
    pycond1 = fycond(l,ig,listaj,gam_j,theta_iter,sigma_iter)
    gam_j[[reg]][pos,ig] = 0
    pycond0 = fycond(l,ig,listaj,gam_j,theta_iter,sigma_iter)
    if (method == 'KUO') {
      aij = pycond1*pij[[reg]][pos]
      bij = pycond0*(1 - pij[[reg]][pos])
    }else if (method == 'SSVS') {
      Xj = listaj$listaX[[reg]]
      yj = c(listaj$listaY[[reg]])
      gam_j[[reg]][pos,ig] = 1
      itauij[[reg]][gam_j[[reg]][,ig] == 0] = tauij[[reg]][gam_j[[reg]][,ig] == 0]
      itauij[[reg]][gam_j[[reg]][,ig] == 1] = 
        cij[[reg]][gam_j[[reg]][,ig] == 1]*tauij[[reg]][gam_j[[reg]][,ig] == 1]
      Dj[[reg]] = diag(itauij[[reg]])
      pthetacond1 = dmnormB(x = theta_iter[[reg]][,ig],mean = rep(0,k*etam[reg]),sigma = Dj[[reg]] %*% Rj[[reg]] %*% Dj[[reg]])
      aij = pycond1*pthetacond1*pij[[reg]][pos]
      gam_j[[reg]][pos,ig] = 0
      itauij[[reg]][gam_j[[reg]][,ig] == 0] = tauij[[reg]][gam_j[[reg]][,ig] == 0]
      itauij[[reg]][gam_j[[reg]][,ig] == 1] = 
        cij[[reg]][gam_j[[reg]][,ig] == 1]*tauij[[reg]][gam_j[[reg]][,ig] == 1]
      Dj[[reg]] = diag(itauij[[reg]])
      pthetacond0 = dmnormB(x = theta_iter[[reg]][,ig],mean = rep(0,k*etam[reg]),sigma = Dj[[reg]] %*% Rj[[reg]] %*% Dj[[reg]])
      bij = pycond0*pthetacond0*(1 - pij[[reg]][pos])
    }
    return(rbinom(1,size = 1,prob = as.numeric((aij)/(aij + bij))))
  }
  listPr = listm[[paste0('m',l)]]
  Dj = list()
  Rj = list()
  i2 = listPr$i
  #creando temporales
  pjmax = listPr$orders$pjmax
  qjmax = listPr$orders$qjmax
  djmax = listPr$orders$djmax
  etam = 1 + pjmax*k + qjmax*nu + djmax
  theta_iter = listPr$Chain$Theta
  sigma_iter = listPr$Chain$Sigma
  gam_iter = listPr$Chain$Gamma
  r_iter = listPr$Chain$r
  theta0j = listPr$Priori$Theta$mean
  sigma0j = listPr$Priori$Theta$cov
  S0j = listPr$Priori$Sigma$cov
  nu0j = listPr$Priori$Sigma$gl
  pij = listPr$Priori$Gamma$prob
  if (method == 'SSVS') {
    itauij = listPr$Priori$Gamma$itauij
    tauij = listPr$Priori$Gamma$tauij
    cij = listPr$Priori$Gamma$cij
    Rj = listPr$Priori$Gamma$Rj
  }
  #iterations update
  listaj = lists(l,r_iter[,i2],pjmax,qjmax,djmax)
  for (lj in 1:l) {
    Xj = listaj$listaX[[lj]]
    Yj = listaj$listaY[[lj]]
    yj = c(Yj)
    Nj = listaj$Nrg[lj]
    if (method == 'SSVS') {
      itauij[[lj]][gam_iter[[lj]][,i2] == 0] = tauij[[lj]][gam_iter[[lj]][,i2] == 0]
      itauij[[lj]][gam_iter[[lj]][,i2] == 1] = cij[[lj]][gam_iter[[lj]][,i2] == 1]*tauij[[lj]][gam_iter[[lj]][,i2] == 1]
      Dj[[lj]] = diag(itauij[[lj]])
      theta0j = rep(0,k*etam[lj])
    }else if (method == 'KUO') {
      Dj[[lj]] = diag(k*etam[lj])
      Rj[[lj]] = sigma0j[[lj]]
    }
    Vj = solve(t(diag(gam_iter[[lj]][,i2])) %*% t(Xj) %*% {diag(Nj) %x% solve(sigma_iter[[lj]][[i2]])} %*% Xj %*% diag(gam_iter[[lj]][,i2]) + solve(Dj[[lj]] %*% Rj[[lj]] %*% Dj[[lj]]))
    thetaj = Vj %*% {t(diag(gam_iter[[lj]][,i2])) %*% t(Xj) %*% {diag(Nj) %x% solve(sigma_iter[[lj]][[i2]])} %*% yj + solve(sigma0j[[lj]]) %*% theta0j[[lj]]}
    theta_iter[[lj]][,i2 + 1] = mvtnorm::rmvnorm(1,mean = thetaj,sigma = Vj)
    Hj = ks::invvec({Xj %*% diag(gam_iter[[lj]][,i2]) %*% theta_iter[[lj]][,i2 + 1]},nrow = k,ncol = Nj)
    Sj = (Yj - Hj) %*% t(Yj - Hj)
    sigma_iter[[lj]][[i2 + 1]] = LaplacesDemon::rinvwishart(nu = Nj + nu0j[[lj]],S = Sj + S0j[[lj]])
    gam_iter[[lj]][,i2 + 1] = gam_iter[[lj]][,i2]
  }
  #gamma
  for (jj in 1:l) {
    for (ii in 1:{k*etam[jj]}) {
      gam_iter[[jj]][ii,i2 + 1] = rgamber(pos = ii,reg = jj,i = i2 + 1)
    }
  }
  #r
  ek = mvtnorm::rmvnorm(1,mean = rep(0,l - 1),sigma = kappa*diag(l - 1))
  rk = r_iter[,i2] + ek
  listark = lists(l,rk,pjmax,qjmax,djmax)
  pr = dmunif(l,rk,a,b)*fycond(l,i2 + 1,listark,gam_iter,theta_iter,sigma_iter)
  px = dmunif(l,r_iter[,i2],a,b)*fycond(l,i2 + 1,listaj,gam_iter,theta_iter,sigma_iter)
  alpha = min(1,as.numeric(pr/px))
  if (alpha >= runif(1)) {
    r_iter[,i2 + 1] = rk
    listr = listark
  }else{
    r_iter[,i2 + 1] = r_iter[,i2]
    listr = listaj
  }
  listPr$Chain = list(Theta = theta_iter, Sigma = sigma_iter,Gamma = gam_iter, r = r_iter)
  listPr$listr = listr
  listPr$i = i2 + 1
  return(listPr)
}
rpseudo = function(l,...){
  listPr = listm[[paste0('m',l)]]
  pjmax = listPr$orders$pjmax
  qjmax = listPr$orders$qjmax
  djmax = listPr$orders$djmax
  etam = 1 + pjmax*k + qjmax*nu + djmax
  theta_iter = sigma_iter = gam_iter = list()
  length(theta_iter) = length(sigma_iter) = length(gam_iter) = l
  for (lj in 1:l) {
    theta_iter[[lj]] = gam_iter[[lj]] = matrix(ncol = 1,nrow = k*etam[lj])
    sigma_iter[[lj]] = list()
    length(sigma_iter[[lj]]) = 1
  }
  r_iter = matrix(nrow = l - 1,ncol = 1)
  theta0jS = listPr$Pseudo$Theta$mean
  sigma0jS = listPr$Pseudo$Theta$cov
  S0jS = listPr$Pseudo$Sigma$cov
  nu0jS = listPr$Pseudo$Sigma$gl
  pijS = listPr$Pseudo$Gamma$prob
  rmeanS = listPr$Pseudo$r$mean
  rcovS = listPr$Pseudo$r$cov
  for (lj in 1:l) {
    theta_iter[[lj]][,1] = mvtnorm::rmvnorm(1,mean = theta0jS[[lj]],sigma = sigma0jS[[lj]])
    sigma_iter[[lj]][[1]] = LaplacesDemon::rwishart(nu = nu0jS[[lj]],S = S0jS[[lj]])
    #sigma_iter[[lj]][[1]] = LaplacesDemon::rinvwishart(nu = nu0jS[[lj]],S = S0jS[[lj]])
    for (iga in 1:{k*etam[lj]}) {
      gam_iter[[lj]][iga,1] = rbinom(n = 1,size = 1,prob = pijS[[lj]][iga])
    }
  }
  r_iter[,1] = mvtnorm::rmvnorm(1,mean = rmeanS, sigma = rcovS)
  listPr2 = list()
  listPr2$Chain = list(Theta = theta_iter, Sigma = sigma_iter,Gamma = gam_iter, r = r_iter)
  listPr2$listr = lists(l,r_iter[,1],pjmax,qjmax,djmax)
  return(listPr2)
}
prodA = function(thetaym, thetaymp, basem, basemp){
  pgammaPn = pthetaPn = psigmaPn = Brobdingnag::as.brob(1)
  pgammaPd = pthetaPd = psigmaPd = Brobdingnag::as.brob(1)
  pgammaSn = pthetaSn = psigmaSn = Brobdingnag::as.brob(1)
  pgammaSd = pthetaSd = psigmaSd = Brobdingnag::as.brob(1)
  lm = length(basem$listr$Nrg)
  lmp = length(basemp$listr$Nrg)
  theta_iterm = thetaym$Chain$Theta
  sigma_iterm = thetaym$Chain$Sigma
  gam_iterm = thetaym$Chain$Gamma
  r_iterm = thetaym$Chain$r
  iA = thetaym$i
  theta_itermp = thetaymp$Chain$Theta
  sigma_itermp = thetaymp$Chain$Sigma
  gam_itermp = thetaymp$Chain$Gamma
  r_itermp = thetaymp$Chain$r
  
  theta0jmp = basemp$Priori$Theta$mean
  sigma0jmp = basemp$Priori$Theta$cov
  S0jmp = basemp$Priori$Sigma$cov
  nu0jmp = basemp$Priori$Sigma$gl
  pijmp = basemp$Priori$Gamma$prob
  theta0jm = basem$Priori$Theta$mean
  sigma0jm = basem$Priori$Theta$cov
  S0jm = basem$Priori$Sigma$cov
  nu0jm = basem$Priori$Sigma$gl
  pijm = basem$Priori$Gamma$prob
  
  theta0jSm = basem$Pseudo$Theta$mean
  sigma0jSm = basem$Pseudo$Theta$cov
  S0jSm = basem$Pseudo$Sigma$cov
  nu0jSm = basem$Pseudo$Sigma$gl
  pijSm = basem$Pseudo$Gamma$prob
  rmeanSm = basem$Pseudo$r$mean
  rcovSm = basem$Pseudo$r$cov
  theta0jSmp = basemp$Pseudo$Theta$mean
  sigma0jSmp = basemp$Pseudo$Theta$cov
  S0jSmp = basemp$Pseudo$Sigma$cov
  nu0jSmp = basemp$Pseudo$Sigma$gl
  pijSmp = basemp$Pseudo$Gamma$prob
  rmeanSmp = basemp$Pseudo$r$mean
  rcovSmp = basemp$Pseudo$r$cov
  for (lj in 1:lmp) {
    pgammaPn = pgammaPn*prodB(Brobdingnag::as.brob(dbinom(gam_itermp[[lj]][,1],size = 1,prob = pijmp[[lj]])))
    pthetaPn = pthetaPn*dmnormB(theta_itermp[[lj]][,1],mean = theta0jmp[[lj]], sigma = sigma0jmp[[lj]])
    psigmaPn = psigmaPn*dinvwishartB(sigma_itermp[[lj]][[1]], nu = nu0jmp[[lj]],S = S0jmp[[lj]])
  }
  prPn = dmunif(lmp,r_itermp[,1],a,b)
  for (lj in 1:lmp) {
    pgammaPd = pgammaPd*prodB(Brobdingnag::as.brob(dbinom(gam_iterm[[lj]][,iA],size = 1,prob = pijm[[lj]])))
    pthetaPd = pthetaPd*dmnormB(theta_iterm[[lj]][,iA],mean = theta0jm[[lj]], sigma = sigma0jm[[lj]])
    psigmaPd = psigmaPd*dinvwishartB(sigma_iterm[[lj]][[iA]], nu = nu0jm[[lj]],S = S0jm[[lj]])
  }
  prPd = dmunif(lm,r_iterm[,iA],a,b)
  fn = fycond(lmp,1,thetaymp$listr,
              thetaymp$Chain$Gamma,thetaymp$Chain$Theta,thetaymp$Chain$Sigma)
  fd = fycond(lm,thetaym$i,thetaym$listr,
              thetaym$Chain$Gamma,thetaym$Chain$Theta,thetaym$Chain$Sigma)
  if ({lm - lmp} > 0) {
    for (lj in 1:lmp) {
      pijSmp[[lj]][pijSmp[[lj]] == 1] = 0.995
      pijSmp[[lj]][pijSmp[[lj]] == 0] = 0.005
      pgammaSn = pgammaSn*prodB(Brobdingnag::as.brob(dbinom(gam_iterm[[lj]][,iA],size = 1,prob = pijSmp[[lj]])))
      pthetaSn = pthetaSn*dmnormB(theta_iterm[[lj]][,iA],mean = theta0jSmp[[lj]], sigma = sigma0jSmp[[lj]])
      psigmaSn = psigmaSn*dwishartB(sigma_iterm[[lj]][[iA]], nu = nu0jSmp[[lj]],S = S0jSmp[[lj]])
    }
    prSn = dmnormB(r_iterm[1:(lmp - 1),iA],mean = rmeanSmp, sigma = rcovSmp)
    for (lj in 1:lmp) {
      pijSm[[lj]][pijSm[[lj]] == 1] = 0.995
      pijSm[[lj]][pijSm[[lj]] == 0] = 0.005
      pgammaSd = pgammaSd*prodB(Brobdingnag::as.brob(dbinom(gam_itermp[[lj]][,1],size = 1,prob = pijSm[[lj]])))
      pthetaSd = pthetaSd*dmnormB(theta_itermp[[lj]][,1],mean = theta0jSm[[lj]], sigma = sigma0jSm[[lj]])
      psigmaSd = psigmaSd*dwishartB(sigma_itermp[[lj]][[1]], nu = nu0jSm[[lj]],S = S0jSm[[lj]])
    }
    prSd = dmnormB(r_itermp[,1],mean = rmeanSm[1:(lmp - 1)], sigma = as.matrix(rcovSm[1:(lmp - 1),1:(lmp - 1)]))
  }
  if ({lm - lmp} < 0) {
    for (lj in 1:lm) {
      pijSmp[[lj]][pijSmp[[lj]] == 1] = 0.995
      pijSmp[[lj]][pijSmp[[lj]] == 0] = 0.005
      pgammaSn = pgammaSn*prodB(Brobdingnag::as.brob(dbinom(gam_iterm[[lj]][,iA],size = 1,prob = pijSmp[[lj]])))
      pthetaSn = pthetaSn*dmnormB(theta_iterm[[lj]][,iA],mean = theta0jSmp[[lj]], sigma = sigma0jSmp[[lj]])
      psigmaSn = psigmaSn*dwishartB(sigma_iterm[[lj]][[iA]], nu = nu0jSmp[[lj]],S = S0jSmp[[lj]])
    }
    prSn = dmnormB(r_iterm[,iA],mean = rmeanSmp[1:(lm - 1)], sigma = as.matrix(rcovSmp[1:(lm - 1),1:(lm - 1)]))
    for (lj in 1:lm) {
      pijSm[[lj]][pijSm[[lj]] == 1] = 0.995
      pijSm[[lj]][pijSm[[lj]] == 0] = 0.005
      pgammaSd = pgammaSd*prodB(Brobdingnag::as.brob(dbinom(gam_itermp[[lj]][,1],size = 1,prob = pijSm[[lj]])))
      pthetaSd = pthetaSd*dmnormB(theta_itermp[[lj]][,1],mean = theta0jSm[[lj]], sigma = sigma0jSm[[lj]])
      psigmaSd = psigmaSd*dwishartB(sigma_itermp[[lj]][[1]], nu = nu0jSm[[lj]],S = S0jSm[[lj]])
    }
    prSd = dmnormB(r_itermp[1:(lm - 1),1],mean = rmeanSm, sigma = as.matrix(rcovSm))
  }
  valn = fn*pgammaPn*pthetaPn*psigmaPn*prPn*pgammaSn*pthetaSn*psigmaSn*prSn
  vald = fd*pgammaPd*pthetaPd*psigmaPd*prPd*pgammaSd*pthetaSd*psigmaSd*prSd
  val = valn/vald 
  return(val)
}
# first values
listm = list()
listm[[paste0('m',2)]] = fill(m = 2,iter = 500)
listm[[paste0('m',3)]] = fill(m = 3,iter = 500)
listm[[paste0('m',4)]] = fill(m = 4,iter = 500, kappa = 0.1)
#
# iterations
m_iter = c()
m_iter[1] = 3
acepm = 0
pb = txtProgressBar(min = 2, max = niter + burn,style = 3)
pm_im = matrix(nrow = l0,ncol = niter + burn - 1)
for (im in 2:{niter + burn}) {
  #simulate m' different of m
  m_iter[im] = rdunif(m_iter[im - 1],l0)
  print(m_iter)
  #Creation m'
  if (!(paste0('m',m_iter[im]) %in% names(listm))) {
    listm[[paste0('m',m_iter[im])]] = fill(m = m_iter[im],kappa = 0.1)
  }
  #m posterior
  listTemp = updatelist(l = m_iter[im - 1])
  #random pseudo m'
  listPseudo = rpseudo(l = m_iter[im])
  #Compute alpha
  val = prodA(listTemp,listPseudo,listm[[paste0('m',m_iter[im - 1])]],listm[[paste0('m',m_iter[im])]])
  alpham = min(1,as.numeric(val))
  if (alpham >= runif(1)) {
    m_iter[im] = m_iter[im]
    acepm = acepm + 1
  }else{
    m_iter[im] = m_iter[im - 1]
    listm[[paste0('m',m_iter[im])]] = listTemp
    acepm = acepm
  }
  pm_im[,im - 1] = table(factor(m_iter,levels = 1:l0))/im
  setTxtProgressBar(pb,im)
}
close(pb)
forecast::autoplot(ts(t(pm_im)),facets = TRUE)

fn;pgammaPn;pthetaPn;psigmaPn;prPn;pgammaSn;pthetaSn;psigmaSn;prSn
fd;pgammaPd;pthetaPd;psigmaPd;prPd;pgammaSd;pthetaSd;psigmaSd;prSd
