# Date: 14/04/2020
# Description:
# Function:
mtarforecast = function(regimemodel,niter,newdata,level = 0.95, chain = FALSE, modelU){
  # validaciones
  if (class(regimemodel) != 'regime-model') {
    stop('regimemodel must be an object of type (regime-model)')
  }
  # objetos necesarios
  Yt = t(regimemodel$data$yt)
  Ut = t(regimemodel$data$Ut)
  nu = nrow(Ut) - 1
  l = length(regimemodel$regime)
  k = nrow(Yt)
  N = ncol(Yt)
  if (is.null(regimemodel$estimates$r)) {
    pj = regimemodel$orders$pj
    qj = regimemodel$orders$qj
    dj = regimemodel$orders$dj
  } else{
    pj = regimemodel$orders$pj$pmax
    qj = regimemodel$orders$qj$qmax
    dj = regimemodel$orders$dj$dmax
  }
  eta = 1 + pj*k + qj*nu + dj
  ## cadenas
  theta_iter = regimemodel$Chain$Theta
  if (!is.null(regimemodel$Chain$Sigma)) {sigma_iter = regimemodel$Chain$Sigma}
  if (!is.null(regimemodel$Chain$r)) {r_iter = regimemodel$Chain$r}
  if (!is.null(regimemodel$Chain$Gamma)) {gam_iter = regimemodel$Chain$Gamma}
  # funciones necesarias
  ### generar Yt
  Ygen = function(t,Y,U,i2...){
    if (!is.null(regimemodel$Chain$r)) {
      r = r_iter[i2]
    }else {r = regimemodel$r}
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
    if (!is.null(regimemodel$Chain$Sigma)) {
      cov = ks::invvec(sigma_iter[[Ind]][,i2])
    }else{cov = regimemodel$regime[[Ind]]$sigma}
    if (!is.null(regimemodel$Chain$Gamma)) {
      mean = {t(wtj) %x% diag(k)} %*% (gam_iter[[Ind]][,i2]*theta_iter[[Ind]][,i2]) 
    }else{
      mean = {t(wtj) %x% diag(k)} %*% (theta_iter[[Ind]][,i2])
    }
    val = MASS::mvrnorm(1,mean,cov %*% cov)
    return(val)
  }
  ### generador de Ut
  b = modelU$orders$pj
  modelU = modelU$regime$R1
  Ugen = function(t,U,...){
    cs = modelU$cs
    At = as.matrix(as.data.frame(modelU$phi))
    Sig = as.matrix(modelU$sigma)
    uti = c()
    for (w in 1:b) {uti = c(uti,U[,t - w])}
    val = MASS::mvrnorm(1, cs + At %*% uti,Sig %*% Sig)
    return(c(val))
  }
  ### Tratamiento para los tiempos
  maxj = max(pj,qj,dj,b)
  if (min(newdata$time) < maxj) {
    Ti = maxj
    Yp = cbind(matrix(0,nrow = k,ncol = maxj),Yt)
    Up = cbind(matrix(0,nrow = nu + 1,ncol = maxj),Ut)
    h = {max(newdata) + maxj} - Ti
    namesYT = paste((1:h) %x% rep(1,k),rep(1,h) %x% (1:k),sep = '.')
    namesUT = paste((1:h) %x% rep(1,nu + 1),rep(1,h) %x% (1:{nu + 1}),sep = '.')
  }else{
    Ti = ifelse(min(newdata$time) < ncol(Yt),min(newdata$time) - 1,ncol(Yt)) 
    h =  max(newdata) - Ti
    Yp = Yt
    Up = Ut
    namesYT = paste(Ti + (1:h) %x% rep(1,k),rep(1,h) %x% (1:k),sep = '.')
    namesUT = paste(Ti + (1:h) %x% rep(1,nu + 1),rep(1,h) %x% (1:{nu + 1}),sep = '.')
  }
  ### Cadenas para las estimaciones
  ChainYt = matrix(ncol = k*h,nrow = niter)
  ChainUt = matrix(ncol = (nu + 1)*h,nrow = niter)
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
    c(paste('lower limit ',(1 - level)/2*100,'%',sep = ''),'mean',paste('upper limit ',(1 + level)/2*100,'%',sep = ''),'RVPD')
  estimUt[,1] = apply(ChainUt,2,quantile,probs = (1 - level)/2)
  estimUt[,2] = apply(ChainUt,2,mean)
  estimUt[,3] = apply(ChainUt,2,quantile,probs = (1 + level)/2)
  estimUt[,4] = apply(ChainUt,2,sd)
  estimYt[,1] = apply(ChainYt,2,quantile,probs = (1 - level)/2)
  estimYt[,2] = apply(ChainYt,2,mean)
  estimYt[,3] = apply(ChainYt,2,quantile,probs = (1 + level)/2)
  estimYt[,4] = apply(ChainYt,2,sd)
  row.names(estimYt) = namesYT
  row.names(estimUt) = namesUT
  ## Salidas
  if (min(newdata$time) >= maxj & length(newdata$time) != ncol(Yt)) {
    estimUt = estimUt[paste(newdata$time %x% rep(1,k), (1:k),sep = '.'),]
    estimYt = estimYt[paste(newdata$time %x% rep(1,k), (1:k),sep = '.'),] 
  }
  Yth = ks::invvec(estimYt[,2],ncol = length(newdata$time),nrow = k)
  Uth = ks::invvec(estimUt[,2],ncol = length(newdata$time),nrow = nu + 1)
  colnames(Yth) = colnames(Uth) = newdata$time
  results = NULL
  results$data$Yt = t(Yt)
  results$data$Ut = t(Ut)
  results$forecast$estim$Ut = estimUt
  results$forecast$estim$Yt = estimYt
  results$forecast$Uth = Uth
  results$forecast$Yth = Yth
  if (chain) {
    results$forecast$chain$Ut = ChainUt
    results$forecast$chain$Yt = ChainYt
  }
  return(results)
}

# Example data(yt,ut):
yt = datasim
estimyt = mtarns(Yt = yt$Yt,Ut = Ut,l = 2,orders = list(pj = yt$pj,dj = yt$dj,qj = yt$qj)
                 ,r = qnorm(0.4),niter = 1000,chain = TRUE)

estimyt = mtarstr(Yt = yt$Yt,Ut = Ut,l = 2,orders = list(pjmax = c(2,2),
                 qjmax = c(1,1),djmax = c(1,1)),niter = 1000,method = 'KUO',chain = T)

Ut = Ut
estimut = mtarns(Yt = Ut,orders = list(pj = 1,qj = 0,dj = 0),niter = 1000)

newdata = data.frame(time = 1:1000)
pred1 = mtarforecast(regimemodel = estimyt,niter = 500,newdata = newdata,modelU = estimut)
