#==================================================================================================#
# Date: 17/04/2019
# Description:
#-> When r_init is NULL default is the proportion that separate the observations in l equal parts
#-> Bayesia estimation with MCMC in order: theta gibbs sampling, sigma gibbs sampling, gamma gibs samplind and
# finally thershold metropolis-hasting with uniform proposal
# Function:
#==================================================================================================#
mtarstr = function(ini_obj, level = 0.95, niter = 1000, burn = NULL, chain = FALSE, r_init = NULL){
  # Just-In-Time (JIT)
  compiler::enableJIT(3)
  # checking
  if (!inherits(ini_obj, 'regim_inipars')) {
    stop('ini_obj must be a regim_inipars object')
  }
  # data
  Yt = ini_obj$tsregim_obj$Yt
  Ut = cbind(ini_obj$tsregim_obj$Zt,ini_obj$tsregim_obj$Xt)
  k = ini_obj$tsregim_obj$k
  N = ini_obj$tsregim_obj$N
  nu = ini_obj$tsregim_obj$nu
  if (is.null(nu)) {nu = 0}
  # parameters
  l = ini_obj$pars$l
  # method
  method = ini_obj$method
  # unknown
  orders = ini_obj$orders
  # code
  burn = ifelse(is.null(burn),round(0.3*niter),burn)
  other = 100
  pjmax = orders$pj
  qjmax = orders$qj
  djmax = orders$dj
  Yt = t(Yt)
  if (l == 1) {
    if (is.null(Ut)) {
      Ut = matrix(0, ncol = N,nrow = 1)
    }else{
      Ut = rbind(matrix(0, ncol = N,nrow = 1),t(Ut)) # only for covariable
    }
  }else{Ut = t(Ut)}
  Zt = Ut[1,]
  if (nu == 0) {
    Xt = matrix(0,ncol = N,nrow = 1)
    qjmax = rep(0,l)
  }else{
    Xt = t(ini_obj$tsregim_obj$Xt)
  }
  eta = 1 + pjmax*k + qjmax*nu + djmax
  # objects for each regimen and iterations
  theta_iter = sigma_iter = gam_iter = pij =  Dj = Rj = vector('list', l)
  if (method == 'SSVS') {
    tauij = itauij = cij = vector('list', l)
  }
  if (l != 1) {r_iter = matrix(ncol = niter + burn + other,nrow = l - 1)}
  # set initial values for each regime in each chain
  itheta0j = isigma0j = iS0j = inu0j = vector('list',l)
  thetaini = ini_obj$init$Theta
  sigmaini = ini_obj$init$Sigma
  gammaini = ini_obj$init$Gamma
  for (lj in 1:l) {
    theta_iter[[lj]] = gam_iter[[lj]] = matrix(ncol = niter + burn + other,nrow = k*eta[lj])
    sigma_iter[[lj]] = vector('list',length = niter + burn + other)
    itheta0j[[lj]] = thetaini[[paste0('R',lj)]]$theta0j
    isigma0j[[lj]] = thetaini[[paste0('R',lj)]]$cov0j
    iS0j[[lj]] = sigmaini[[paste0('R',lj)]]$S0j
    inu0j[[lj]] = sigmaini[[paste0('R',lj)]]$nu0j
    pij[[lj]] = gammaini[[paste0('R',lj)]]
    gam_iter[[lj]][,1] = rep(1,k*eta[lj])
    if (method == 'KUO') {
      theta_iter[[lj]][,1] = mvtnorm::rmvnorm(1,mean = itheta0j[[lj]],sigma = isigma0j[[lj]])
    }
    if (method == 'SSVS') {
      cij[[lj]] = ini_obj$init$Theta[[paste0('R',lj)]]$Cij
      tauij[[lj]] = ini_obj$init$Theta[[paste0('R',lj)]]$Tauij
      Rj[[lj]] = ini_obj$init$Theta[[paste0('R',lj)]]$R
      itauij[[lj]] = cij[[lj]]*tauij[[lj]]
      Dj[[lj]] = diag(itauij[[lj]])
      theta_iter[[lj]][,1] = mvtnorm::rmvnorm(1,mean = rep(0,k*eta[lj]),sigma = Dj[[lj]] %*% Rj[[lj]] %*% Dj[[lj]])
    }
    sigma_iter[[lj]][[1]] = MCMCpack::riwish(v = inu0j[[lj]],S = iS0j[[lj]])
  }
  #
  # objects for save regimes, chains and credibility intervals
  Rest = thetaest = thetachain =
    gamest = gamchain = sigmaest = sigmachain = vector('list', l)
  names(Rest) = names(thetaest) = names(thetachain) = names(gamest) = names(gamchain) =
    names(sigmaest) = names(sigmachain) = paste0('R',1:l)
  # necesary functions
  fycond = function(i2,listr,gamma,theta_iter,sigma_iter,l){
    acum = 0
    Nrg = listr$Nrg
    for (lj in 1:l) {
      yj = c(listr$listaY[[lj]])
      Xj = listr$listaX[[lj]]
      acum = acum + t(yj - Xj %*% diag(gamma[[lj]][,i2]) %*% theta_iter[[lj]][,i2]) %*%
        {diag(Nrg[lj]) %x% solve(sigma_iter[[lj]][[i2]])} %*%
        (yj - Xj %*% diag(gamma[[lj]][,i2]) %*% theta_iter[[lj]][,i2])
    }
    sigmareg = lapply(sigma_iter,function(x){x[[i2]]})
    cte = prodB(Brobdingnag::as.brob(sapply(sigmareg,function(x){
      return(c(determinant(x,logarithm = FALSE)$modulus))}))^{-Nrg/2})
    val = cte*exp(-1/2*Brobdingnag::as.brob(acum))
    return(val)
  }
  fycond = compiler::cmpfun(fycond)
  dmunif = function(r,a,b){
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
  dmunif = compiler::cmpfun(dmunif)
  lists = function(r,...){
    rj = matrix(nrow = 2,ncol = l)
    if (l == 1) {
      rj[,1] = c(-Inf,Inf)
    }else{
      rj[,1] = c(-Inf,r[1])
      rj[,l] = c(rev(r)[1],Inf)
    }
    if (l > 2) {for (i2 in 2:{l - 1}) {rj[,i2] = c(r[i2 - 1],r[i2])}}
    # indicator variable for the regime
    Ind = vector(mode = 'numeric',length = N)
    for (j in 1:l) {
      Ind[Zt > rj[1,j] & Zt <= rj[2,j]] = j
    }
    Nrg = vector(mode = 'numeric')
    listaXj = listaYj = vector('list', l)
    Inj_X = function(ti,Yt,Zt,Xt,p,q,d){
      yti = vector(mode = "numeric")
      for (w in 1:p) {yti = c(yti,Yt[,ti - w])}
      xti = vector(mode = "numeric")
      for (w in 1:q) {xti = c(xti,Xt[,ti - w])}
      zti = vector(mode = "numeric")
      for (w in 1:d) {zti = c(zti,Zt[ti - w])}
      if (q == 0 & d != 0) {
        wtj = c(1,yti,zti)
      }else if (d == 0 & q != 0) {
        wtj = c(1,yti,xti)
      }else if (d == 0 & q == 0) {
        wtj = c(1,yti)
      }else{
        wtj = c(1,yti,xti,zti)}
      return(wtj)
    }
    Inj_X = Vectorize(Inj_X,vectorize.args = "ti")
    for (lj in 1:l) {
      p = pjmax[lj]
      q = qjmax[lj]
      d = djmax[lj]
      maxj = max(p,q,d)
      Inj = which(Ind == lj)
      Inj = Inj[Inj > maxj]
      Nrg[lj] = length(Inj)
      Yj = matrix(Yt[,Inj],nrow = k,ncol = Nrg[lj])
      yj = c(Yj)
      Wj = sapply(Inj,Inj_X,Yt = Yt,Zt = Zt,Xt = Xt,p = p,q = q,d = d)
      Xj = t(Wj) %x% diag(k)[1,]
      for (s in 2:k) {Xj = cbind(Xj,t(Wj) %x% diag(k)[s,])}
      listaXj[[lj]] = Xj
      listaYj[[lj]] = Yj
    }
    return(list(Nrg = Nrg,listaX = listaXj,listaY = listaYj))
  }
  lists = compiler::cmpfun(lists)
  rgamber = function(pos,reg,i,listj,theta_iter,sigma_iter,gam_iter,...){
    gam_j = gam_iter
    gam_j[[reg]][pos,i] = 1
    pycond1 = fycond(i,listj,gam_j,theta_iter,sigma_iter,l)
    gam_j[[reg]][pos,i] = 0
    pycond0 = fycond(i,listj,gam_j,theta_iter,sigma_iter,l)
    if (method == 'KUO') {
      aij = pycond1*pij[[reg]][pos]
      bij = pycond0*(1 - pij[[reg]][pos])
    }else if (method == 'SSVS') {
      gam_j[[reg]][pos,i] = 1
      itauij[[reg]][gam_j[[reg]][,i] == 0] = tauij[[reg]][gam_j[[reg]][,i] == 0]
      itauij[[reg]][gam_j[[reg]][,i] == 1] =
        cij[[reg]][gam_j[[reg]][,i] == 1]*tauij[[reg]][gam_j[[reg]][,i] == 1]
      Dj[[reg]] = diag(itauij[[reg]])
      pthetacond1 = dmnormB(x = theta_iter[[reg]][,i],mean = rep(0,k*eta[reg]),sigma = Dj[[reg]] %*% Rj[[reg]] %*% Dj[[reg]])
      aij = pycond1*pthetacond1*pij[[reg]][pos]
      gam_j[[reg]][pos,i] = 0
      itauij[[reg]][gam_j[[reg]][,i] == 0] = tauij[[reg]][gam_j[[reg]][,i] == 0]
      itauij[[reg]][gam_j[[reg]][,i] == 1] =
        cij[[reg]][gam_j[[reg]][,i] == 1]*tauij[[reg]][gam_j[[reg]][,i] == 1]
      Dj[[reg]] = diag(itauij[[reg]])
      pthetacond0 = dmnormB(x = theta_iter[[reg]][,i],mean = rep(0,k*eta[reg]),sigma = Dj[[reg]] %*% Rj[[reg]] %*% Dj[[reg]])
      bij = pycond0*pthetacond0*(1 - pij[[reg]][pos])
    }
    return(rbinom(1,size = 1,prob = as.numeric((aij)/(aij + bij))))
  }
  rgamber = compiler::cmpfun(Vectorize(rgamber,"pos"))
  # edges for rmunif
  rini = ini_obj$init$r
  a = ifelse(is.null(rini$za),min(Zt),quantile(Zt,probs = rini$za))
  b = ifelse(is.null(rini$zb),max(Zt),quantile(Zt,probs = rini$zb))
  #
  # last check
  if (!is.null(r_init)) {
    if (is.numeric(r_init) & length(r_init) == {l - 1}) {
      if (dmunif(r_init,a,b) == 0) {
        stop('r_init must be in Zt range and for l >= 2, r[i] < r[i+1]')
      }
    }else{stop('r_init must be vector of length l - 1')}
  }
  if (l != 1) {
    if (is.null(r_init)) {
      r_iter[,1] = c(quantile(Zt, probs = 1/l*(1:{l - 1})))
    }else{
      r_iter[,1] = r_init
    }
  }
  list_m = list(r_iter = r_iter,theta_iter = theta_iter,sigma_iter = sigma_iter,gam_iter = gam_iter)
  # iterations function
  iter_str = function(i, list_m, ...){
    r_iter = list_m$r_iter
    theta_iter = list_m$theta_iter
    sigma_iter = list_m$sigma_iter
    gam_iter = list_m$gam_iter
    if (l != 1) {
      listj = lists(r_iter[,i - 1])
    }else{
      listj = lists(0)
    }
    for (lj in 1:l) {
      Xj = listj$listaX[[lj]]
      Yj = listj$listaY[[lj]]
      yj = c(Yj)
      Nj = listj$Nrg[lj]
      theta0j = itheta0j[[lj]]
      sigma0j = isigma0j[[lj]]
      S0j = iS0j[[lj]]
      nu0j = inu0j[[lj]]
      if (method == 'SSVS') {
        itauij[[lj]][gam_iter[[lj]][,i - 1] == 0] = tauij[[lj]][gam_iter[[lj]][,i - 1] == 0]
        itauij[[lj]][gam_iter[[lj]][,i - 1] == 1] = cij[[lj]][gam_iter[[lj]][,i - 1] == 1]*tauij[[lj]][gam_iter[[lj]][,i - 1] == 1]
        Dj[[lj]] = diag(itauij[[lj]])
        theta0j = rep(0,k*eta[lj])
      }else if (method == 'KUO') {
        Dj[[lj]] = diag(k*eta[lj])
        Rj[[lj]] = sigma0j
      }
      Vj = solve(t(diag(gam_iter[[lj]][,i - 1])) %*% t(Xj) %*% {diag(Nj) %x% solve(sigma_iter[[lj]][[i - 1]])} %*% Xj %*% diag(gam_iter[[lj]][,i - 1]) + solve(Dj[[lj]] %*% Rj[[lj]] %*% Dj[[lj]]))
      thetaj = Vj %*% {t(diag(gam_iter[[lj]][,i - 1])) %*% t(Xj) %*% {diag(Nj) %x% solve(sigma_iter[[lj]][[i - 1]])} %*% yj + solve(sigma0j) %*% theta0j}
      theta_iter[[lj]][,i] = mvtnorm::rmvnorm(1,mean = thetaj,sigma = Vj)
      Hj = ks::invvec({Xj %*% diag(gam_iter[[lj]][,i - 1]) %*% theta_iter[[lj]][,i]},nrow = k,ncol = Nj)
      Sj = (Yj - Hj) %*% t(Yj - Hj)
      sigma_iter[[lj]][[i]] = MCMCpack::riwish(v = Nj + nu0j,S = Sj + S0j)
      gam_iter[[lj]][,i] = gam_iter[[lj]][,i - 1]
    }
    for (jj in 1:l) {
      gam_iter[[jj]][,i] = rgamber(pos = 1:{k*eta[jj]},reg = jj,i = i,listj = listj,theta_iter,sigma_iter,gam_iter)
    }
    if (l != 1) {
      if (i < other) {
        ek = mvtnorm::rmvnorm(1,mean = rep(0,l - 1),sigma = 0.5*diag(l - 1))
      }else{
        ek = runif(l - 1,-abs(rini$val_rmh),abs(rini$val_rmh))
      }
      rk = r_iter[,i - 1] + ek
      listrk = lists(rk)
      pr = dmunif(rk,a,b)*fycond(i,listrk,gam_iter,theta_iter,sigma_iter,l)
      px = dmunif(r_iter[,i - 1],a,b)*fycond(i,listj,gam_iter,theta_iter,sigma_iter,l)
      alpha = min(1,as.numeric(pr/px))
      if (alpha >= runif(1)) {
        r_iter[,i] = rk
        acep = acep + 1
      }else{
        r_iter[,i] = r_iter[,i - 1]
        acep = acep
      }
    }
    list_m$r_iter = r_iter
    list_m$theta_iter = theta_iter
    list_m$sigma_iter = sigma_iter
    list_m$gam_iter = gam_iter
    return(list_m)
  }
  iter_str = compiler::cmpfun(iter_str)
  # iterations: gibbs and metropolis sampling
  acep = 0
  cat('Estimating threshold(s), structural and non-structural parameters ...','\n')
  pb = txtProgressBar(min = 2, max = niter + burn + other,style = 3)
  for (i in 2:{niter + burn + other}) {
    list_m = iter_str(i,list_m)
    setTxtProgressBar(pb,i)
  }
  close(pb)
  cat('Saving results ... \n')
  r_iter = list_m$r_iter[-c(1:other)]
  theta_iter = list_m$theta_iter
  sigma_iter = list_m$sigma_iter
  gam_iter = list_m$gam_iter
  # exits
  if (l != 1) {
    rest = matrix(nrow = l - 1,ncol = 3)
    colnames(rest) = colnames(rest) =
      c(paste('lower limit ',(1 - level)/2*100,'%',sep = ''),'mean',paste('upper limit ',(1 + level)/2*100,'%',sep = ''))
    rchain = matrix(r_iter[-c(1:burn)],ncol = niter,nrow = l - 1)
    rest[,1] = apply(rchain,1,quantile,probs = (1 - level)/2)
    rest[,3] = apply(rchain,1,quantile,probs = (1 + level)/2)
    rest[,2] = apply(rchain,1,mean)
    rvec = c('mean' = rest[,2],'prop %' = acep/niter*100)
  }
  # save chains
  # logLik
  listj = lists(rest[,2])
  SigmaPrep = function(x){return(c(expm::sqrtm(matrix(x,k,k))))}
  logLikj = vector(mode = "numeric")
  for (lj in 1:l) {
    theta_iter[[lj]] = theta_iter[[lj]][,-c(1:other)]
    gam_iter[[lj]] = gam_iter[[lj]][,-c(1:other)]
    sigma_iter[[lj]] = sigma_iter[[lj]][-c(1:other)]
    #
    thetachain[[lj]] = theta_iter[[lj]][,-c(1:burn)]
    gamchain[[lj]] = gam_iter[[lj]][,-c(1:burn)]
    sigmachain[[lj]] = sapply(sigma_iter[[lj]][-c(1:burn)], ks::vec)
    # credibility intervals
    vecgamma = vectheta = matrix(nrow = k*eta[lj],ncol = 3)
    vecsigma = matrix(nrow = k*k,ncol = 3)
    colnames(vectheta) = colnames(vecsigma) =
      c(paste0('lower limit ',(1 - level)/2*100,'%'),'mean',paste0('upper limit ',(1 + level)/2*100,'%'))
    rownames(vecsigma) = c(sapply(1:k, function(x){paste0(1:k,x)}))
    if (nu != 0 & qjmax[lj] != 0 & djmax[lj] != 0) {
      rownames(vectheta) = rownames(vecgamma) =
        rep(c('phi0',rep(paste0('phi',1:pjmax[lj]),each = k),rep(paste0('beta',1:qjmax[lj]),each = nu),paste0('delta',1:djmax[lj])),k)
    }else if (nu != 0 & qjmax[lj] != 0 & djmax[lj] == 0) {
      rownames(vectheta) = rownames(vecgamma) =
        rep(c('phi0',rep(paste0('phi',1:pjmax[lj]),each = k),rep(paste0('beta',1:qjmax[lj]),each = nu)),k)
    }else if (qjmax[lj] == 0 & djmax[lj] != 0) {
      rownames(vectheta) = rownames(vecgamma) =
        rep(c('phi0',rep(paste0('phi',1:pjmax[lj]),each = k),paste0('delta',1:djmax[lj])),k)
    }else if (qjmax[lj] == 0 & djmax[lj] == 0) {
      rownames(vectheta) = rownames(vecgamma) =
        rep(c('phi0',rep(paste0('phi',1:pjmax[lj]),each = k)),k)
    }
    rownames(vecsigma) = c(sapply(1:k, function(x){paste0(1:k,x)}))
    vectheta[,1] = apply(thetachain[[lj]],1,quantile,probs = (1 - level)/2)
    vectheta[,3] = apply(thetachain[[lj]],1,quantile,probs = (1 + level)/2)
    vectheta[,2] = apply(thetachain[[lj]],1,mean)
    thetaest[[lj]] = vectheta
    if (k == 1) {
      vecsigma[,1] = sqrt(quantile(sigmachain[[lj]],probs = (1 - level)/2))
      vecsigma[,3] = sqrt(quantile(sigmachain[[lj]],probs = (1 + level)/2))
      vecsigma[,2] = sqrt(mean(sigmachain[[lj]]))
    }else{
      vecsigma[,1] = SigmaPrep(apply(sigmachain[[lj]],1,quantile,probs = (1 - level)/2))
      vecsigma[,3] = SigmaPrep(apply(sigmachain[[lj]],1,quantile,probs = (1 + level)/2))
      vecsigma[,2] = SigmaPrep(apply(sigmachain[[lj]],1,mean))
    }
    sigmaest[[lj]] = vecsigma
    vec = apply(gamchain[[lj]],2,paste,collapse = '')
    vecs = sort(table(vec), decreasing = TRUE)[1:3]
    colnames(vecgamma) = c(paste('first freq',paste0(vecs[1]/niter*100,'%')),
                           paste('second freq',paste0(vecs[2]/niter*100,'%')),
                           paste('third freq',paste0(vecs[3]/niter*100,'%')))
    vecgamma[,1] = gamchain[[lj]][,which(vec == names(vecs[1]))[1]]
    vecgamma[,2] = gamchain[[lj]][,which(vec == names(vecs[2]))[1]]
    vecgamma[,3] = gamchain[[lj]][,which(vec == names(vecs[3]))[1]]
    gamest[[lj]] = vecgamma
    # creation of the 'regime' type object
    p = pjmax[lj]
    q = qjmax[lj]
    d = djmax[lj]
    if (q == 0 & d == 0) {
      thetaind = c(90,(10 + (1:p)) %x% rep(1,k))
    }else if (q != 0 & d == 0) {
      thetaind = c(90,(10 + (1:p)) %x% rep(1,k),(20 + (1:q)) %x% rep(1,nu))
    }else if (q == 0 & d != 0) {
      thetaind = c(90,(10 + (1:p)) %x% rep(1,k),30 + (1:d))
    }else{
      thetaind = c(90,(10 + (1:p)) %x% rep(1,k),(20 + (1:q)) %x% rep(1,nu),30 + (1:d))
    }
    Thetaj = t(ks::invvec(thetaest[[lj]][,2],ncol = k,nrow = eta[lj]))*t(ks::invvec(gamest[[lj]][,1],ncol = k,nrow = eta[lj]))
    Ri = list()
    cs = matrix(Thetaj[,thetaind == 90],nrow = k,ncol = 1)
    if (sum(cs == 0) != k) {Ri$cs = cs}
    phiest = vector('list', p)
    names(phiest) = paste0('phi',1:p)
    for (j in 1:p) {phiest[[j]] = matrix(Thetaj[,thetaind == (10 + j)],nrow = k,ncol = k)}
    pest = sapply(phiest,function(x){sum(x == 0) != k*k})
    if (sum(pest) != 0) {
      Ri$phi = phiest[pest]
    }
    if (q != 0) {
      betaest = vector('list', q)
      names(betaest) = paste0('beta',1:q)
      for (j in 1:q) {betaest[[j]] = matrix(Thetaj[,thetaind == (20 + j)],nrow = k,ncol = nu)}
      qest = sapply(betaest,function(x){sum(x == 0) != k*nu})
      if (sum(qest) != 0) {
        Ri$beta = betaest[qest]
      }
    }
    if (d != 0) {
      deltaest = vector('list', d)
      names(deltaest) = paste0('delta',1:d)
      for (j in 1:d) {deltaest[[j]] = matrix(Thetaj[,thetaind == (30 + j)],nrow = k,ncol = 1)}
      dest = sapply(deltaest,function(x){sum(x == 0) != k})
      if (sum(dest) != 0) {
        Ri$delta = deltaest[dest]
      }
    }
    Ri$sigma = ks::invvec(sigmaest[[lj]][,2],ncol = k,nrow = k)
    Rest[[lj]] = mtaregim(orders = list(p = length(Ri$phi),q = length(Ri$beta),d = length(Ri$delta)),cs = Ri$cs,
                          Phi = Ri$phi,Beta = Ri$beta,Delta = Ri$delta,
                          Sigma = Ri$sigma)
    Xj = listj$listaX[[lj]]
    Yj = listj$listaY[[lj]]
    Nj = listj$Nrg[lj]
    Hj = ks::invvec({Xj %*% diag(gamest[[lj]][,1]) %*% thetaest[[lj]][,2]},nrow = k,ncol = Nj)
    Sj = (Yj - Hj) %*% t(Yj - Hj)
    logLikj[lj] = log(det(Sj/Nj))
  }
  # exits
  ## chain
  if (chain) {
    Chain = vector('list')
    Chain$Theta = thetachain
    Chain$Gamma = gamchain
    Chain$Sigma = sigmachain
    if (l != 1) {Chain$r = rchain}
  }
  ## estimates
  estimates = vector('list')
  estimates$Theta = thetaest
  estimates$Gamma = gamest
  estimates$Sigma = sigmaest
  data = vector('list')
  data$yt = t(Yt)
  data$Ut = t(Ut)
  orders = vector('list')
  orders$pj = sapply(Rest,function(x){length(x$phi)})
  orders$qj = sapply(Rest,function(x){length(x$beta)})
  orders$dj = sapply(Rest,function(x){length(x$delta)})
  if (l != 1) {estimates$r = rest}
  if (chain) {
    results = list(Nj = listj$Nrg,estimates = estimates,regime = Rest,Chain = Chain,logLikj = logLikj,data = data,r = rvec,orders = orders)
    class(results) = 'regim_model'
  }else{
    results = list(Nj = listj$Nrg,estimates = estimates,regime = Rest,logLikj = logLikj,data = data,r = rvec,orders = orders)
    class(results) = 'regim_model'
  }
  compiler::enableJIT(0)
  return(results)
}
