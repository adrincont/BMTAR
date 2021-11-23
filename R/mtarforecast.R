# Date: 18/07/2019
# Description:
# Function:
mtarforecast = function(object, ...) UseMethod("mtarforecast")
mtarforecast.regime_model = function(regimemodel,h,level = 0.95, chain = TRUE,  b = NULL){
  niter = dim(regimemodel$Chain$Theta$R1)[2] # Toma las iteraciones de regimemodel
  burn = round(0.3*niter)
  # validaciones
  if (class(regimemodel) != "regime_model") {
    stop("regimemodel must be an object of type (regime_model)")
  }
  if (!is.numeric(h)){
    stop("h must be a numeric object")
  }else{
    if (!(round(h) == h)){stop("h must be a integer numeric object")}
  }
  if (!is.numeric(b)){
    stop("b must be a numeric object")
  }else{
    if (!(round(b) == b)){stop("b must be a integer numeric object")}
  }
  if (!is.numeric(level)){stop("level must be a numeric object between (0,1)")}
  if (!(identical(chain,FALSE) | identical(chain,TRUE))){stop("chain must be TRUE or FALSE")}
  # Varibles auxiliares
  b = ifelse(is.null(b),1,b)
  newdata = data.frame(time = seq(regimemodel$data$N + 1,regimemodel$data$N + h))
  Yt = t(regimemodel$data$Yt)
  Ut = t(cbind(regimemodel$data$Zt,regimemodel$data$Xt))
  nu = nrow(Ut) - 1
  l = length(regimemodel$regime)
  k = nrow(Yt)
  N = ncol(Yt)
  if (is.null(regimemodel$estimates$Gamma)) {
    pj = regimemodel$orders$pj
    qj = regimemodel$orders$qj
    dj = regimemodel$orders$dj
  } else{
    pj = regimemodel$initial$orders$pj
    qj = regimemodel$initial$orders$qj
    dj = regimemodel$initial$orders$dj
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
      r = r_iter[,i2]
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
    yti = vector(mode = "numeric")
    for (w in 1:p) {
      if (t - w > 0) {yti = c(yti,Y[,t - w])
      }else{yti = c(yti,rep(0,k))}}
    if (identical(yti,numeric(0))) {yti = rep(0,k*p)}
    xti = vector(mode = "numeric")
    for (w in 1:q) {
      if (t - w > 0) {xti = c(xti,U[-1,t - w])
      }else{xti = c(xti,rep(0,nu))}}
    if (identical(xti,numeric(0))) {xti = rep(0,nu*q)}
    zti = vector(mode = "numeric")
    for (w in 1:d) {
      if (t - w > 0) {zti = c(zti,U[1,t - w])
      }else{zti = c(zti,0)}}
    if (identical(zti,numeric(0))) {zti = rep(0,d)}
    if (q == 0 & d != 0) {
      wtj = c(1,yti,zti)
    }else if (d == 0 & q != 0) {
      wtj = c(1,yti,xti)
    }else if (d == 0 & q == 0) {
      wtj = c(1,yti)
    }else{
      wtj = c(1,yti,xti,zti)}
    if (!is.null(regimemodel$Chain$Sigma)) {
      if (k == 1){
        cov = ks::invvec(sigma_iter[[Ind]][i2])
      }else{cov = ks::invvec(sigma_iter[[Ind]][,i2])}
    }else{cov = regimemodel$regime[[Ind]]$sigma}
    Xj = t(wtj) %x% diag(k)[1,]
    if (k != 1) {for (s in 2:k) {Xj = cbind(Xj,t(wtj) %x% diag(k)[s,])}}
    if (!is.null(regimemodel$Chain$Gamma)) {
      mean = Xj %*% diag(gam_iter[[Ind]][,i2]) %*% theta_iter[[Ind]][,i2]
    }else{
      mean = Xj %*% (theta_iter[[Ind]][,i2])
    }
    val = MASS::mvrnorm(1,mean,cov %*% cov)
    return(val)
  }
  ### generador de Ut
  initialU = mtarinipars(tsregime(t(Ut)),list_model = list(pars = list(l = 1,orders = list(pj = b,qj = 0,dj = 0))))
  message('Estimating model Ut = (Zt,Xt) \n')
  modelU = mtarns(ini_obj = initialU,niter = niter + 100,chain = FALSE,burn = 1000)
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
  if (min(newdata$time) <= maxj) {
    Ti = maxj
    Yp = cbind(matrix(0,nrow = k,ncol = maxj),Yt)
    Up = cbind(matrix(0,nrow = nu + 1,ncol = maxj),Ut)
    h = {max(newdata) + maxj} - Ti
    namesYT = paste((1:h) %x% rep(1,k),rep(1,h) %x% (1:k),sep = ".")
    if(is.null(regimemodel$data$Xt)) {
      namesUT = paste(1:h,1,sep = ".")
      namesZT = namesUT
    }else{
      namesUT = paste((1:h) %x% rep(1,nu + 1),rep(1,h) %x% (1:{nu + 1}),sep = ".")
      namesZT = paste(1:h,1,sep = ".")
      namesXt = paste((1:h) %x% rep(1,nu),rep(1,h) %x% (2:{nu + 1}),sep = ".")
      namesXt_def = paste((1:h) %x% rep(1,nu),rep(1,h) %x% 1:nu,sep = ".")
    }
  }else{
    Ti = ifelse(min(newdata$time) < ncol(Yt),min(newdata$time) - 1,ncol(Yt))
    h =  max(newdata) - Ti
    Yp = Yt
    Up = Ut
    namesYT = paste(Ti + (1:h) %x% rep(1,k),rep(1,h) %x% (1:k),sep = ".")
    if(is.null(regimemodel$data$Xt)) {
      namesUT = paste(Ti + (1:h),1,sep = ".")
      namesZT = namesUT
    }else{
      namesUT = paste(Ti + (1:h) %x% rep(1,nu + 1),rep(1,h) %x% (1:{nu + 1}),sep = ".")
      namesZT = paste(Ti + (1:h),1,sep = ".")
      namesXt = paste(Ti + (1:h) %x% rep(1,nu),rep(1,h) %x% (2:{nu + 1}),sep = ".")
      namesXt_def = paste(Ti + (1:h) %x% rep(1,nu),rep(1,h) %x% 1:nu,sep = ".")
    }
  }
  ### Cadenas para las estimaciones
  ChainYt = matrix(ncol = k*h,nrow = niter)
  ChainUt = matrix(ncol = (nu + 1)*h,nrow = niter)
  Upi = matrix(nrow = nrow(Ut),ncol = h)
  Ypi = matrix(nrow = nrow(Yt),ncol = h)
  ### Proceso iterativo
  pb = pbapply::timerProgressBar(min = 1, max = niter,style = 1)
  for (i2 in 1:niter) {
    for (i3 in 1:h) {
      Upi[,i3] = Ugen(Ti + i3,Up)
      Up = cbind(Up,Upi[,i3])
      Ypi[,i3] = Ygen(Ti + i3,Yp,Up,i2)
      Yp = cbind(Yp,Ypi[,i3])
    }
    ChainYt[i2,] = ks::vec(Ypi)
    ChainUt[i2,] = ks::vec(Upi)
    pbapply::setTimerProgressBar(pb,i2)
  }
  ### Calculo de estimaciones y intervalos
  # Estimaciones de Yt
  estimYt = matrix(nrow = k*h,ncol = 4)
  estimYt[,1] = apply(ChainYt[-c(1:burn),],2,quantile,probs = (1 - level)/2)
  estimYt[,2] = apply(ChainYt[-c(1:burn),],2,mean)
  estimYt[,3] = apply(ChainYt[-c(1:burn),],2,quantile,probs = (1 + level)/2)
  estimYt[,4] = apply(ChainYt[-c(1:burn),],2,sd)
  j = 1
  UF = YF = vector('numeric',h)
  for (i5 in seq(1,k*h,k)) {
    YF[j] = norm(cov(as.matrix(ChainYt[-c(1:burn),i5:(i5 + k - 1)])),type = 'F')
    j = j + 1
  }
  row.names(estimYt) = namesYT
  # Estimaciones de Ut
  estimUt = matrix(nrow = (nu + 1)*h,ncol = 4)
  estimUt[,1] = apply(ChainUt[-c(1:burn),],2,quantile,probs = (1 - level)/2)
  estimUt[,2] = apply(ChainUt[-c(1:burn),],2,mean)
  estimUt[,3] = apply(ChainUt[-c(1:burn),],2,quantile,probs = (1 + level)/2)
  estimUt[,4] = apply(ChainUt[-c(1:burn),],2,sd)
  j = 1
  for (i4 in seq(1,(nu + 1)*h,nu + 1)) {
    if (nu == 0){
      UF[j] = norm(cov(as.matrix(ChainUt[-c(1:burn),i4:(i4 + nu)])),type = 'F')
    }else{
      UF[j] = norm(cov(ChainUt[-c(1:burn),i4:(i4 + nu)]),type = 'F')
    }
    j = j + 1
  }
  row.names(estimUt) = namesUT
  colnames(estimUt) = colnames(estimYt) =
    c(paste('lower limit ',(1 - level)/2*100,'%',sep = ""),'mean',paste('upper limit ',(1 + level)/2*100,'%',sep = ""),"FNDP")
  estimUt[,'FNDP'] = UF %*% rep(1,nu)
  estimYt[,'FNDP'] = UF %*% rep(1,k)
  if(is.null(regimemodel$data$Xt)) {
    estimZt = estimUt
  }else{
    # Estimaciones de Zt
    estimZt = estimUt[namesZT,]
    # Estimaciones de Xt
    if(nu != 0) {
      estimXt = estimUt[namesXt,]
      row.names(estimXt) = namesXt_def
      }
  }
  ## Salidas
  if (min(newdata$time) >= maxj & length(newdata$time) != ncol(Yt)) {
    estimUt = estimUt[paste(newdata$time %x% rep(1,(nu + 1)), (1:(nu + 1)),sep = "."),]
    estimYt = estimYt[paste(newdata$time %x% rep(1,k), (1:k),sep = "."),]
    if(is.null(regimemodel$data$Xt)) {
      estimZt = estimUt
    }else{
      estimZt = estimUt[namesZT,]
      if (nu != 0){
        estimXt = estimUt[namesXt,]
        row.names(estimXt) = namesXt_def
        }
    }
  }
  Yth = ks::invvec(estimYt[,2],ncol = length(newdata$time),nrow = k)
  Uth = ks::invvec(estimUt[,2],ncol = length(newdata$time),nrow = nu + 1)
  colnames(Yth) = colnames(Uth) = newdata$time
  results = NULL
  results$forecast$Yth = Yth
  results$forecast$estim$Yt = estimYt
  if(!is.null(regimemodel$data$Zt)){
    if(nu != 0){
      Zth = Uth[1,]
      Xth = Uth[-1,]
      results$forecast$estim$Xt = estimXt
      results$forecast$estim$Zt = estimZt
      results$forecast$Zth = Zth
      results$forecast$Xth = Xth
      if(nu == 1){
        ts_reg_ht = tsregime(Yt = t(cbind(Yt,Yth)),Zt = c(Ut[1,],Uth[1,]),Xt = c(Ut[-1,],Uth[-1,]),r = regimemodel$r)
      }else{
          ts_reg_ht = tsregime(Yt = t(cbind(Yt,Yth)),Zt = c(Ut[1,],Uth[1,]),Xt = t(cbind(Ut[-1,],Uth[-1,])),r = regimemodel$r)
      }
    }else{
      Zth = Uth
      results$forecast$estim$Zt = estimUt
      results$forecast$Zth = Uth
      ts_reg_ht = tsregime(Yt = t(cbind(Yt,Yth)),Zt = t(cbind(Ut,Uth)),r = regimemodel$r)
    }
  }else{
    ts_reg_ht = tsregime(Yt = t(cbind(Yt,Yth)))
  }
  results$forecast$estim_Ut = estimUt
  results$tsregime = ts_reg_ht
  results$FNDP$Yt = YF
  results$FNDP$Ut = UF
  if (chain) {
    results$forecast$chain$Ut = ChainUt[-c(1:burn),]
    results$forecast$chain$Yt = ChainYt[-c(1:burn),]
  }
  class(results) = 'regime_forecast'
  return(results)
}
