# Date: 07/04/2019
# Description:
#-> Comentar como iniciamos r si no es dado
#-> Comentar que hacemos si r desconocido en este caso
# Function:
mtarns = function(Yt, Ut = NULL, l = 1, orders, r = NULL, Sigma = NULL, level = 0.95, burn = NULL, niter,
                  chain = FALSE, thetaini = NULL, sigmaini = NULL,
                  rini = list(za = NULL, zb = NULL, kappa = NULL, ini = NULL)){
  if (l != 1 & is.null(Ut)) {stop('Threshold process Ut (optional covariable) must be enter')}
  if (l != 1 & !is.matrix(Ut)) {stop('Ut must be a matrix type object')}
  if (!is.matrix(Yt)) {stop('Yt must be a matrix type object')}
  if (!is.list(orders) | length(orders) != 3) {
    stop('orders must be a list of length 3 list(pj, qj, dj)')
  }else if (sum(names(orders) %in% c("pj","qj","dj")) != 3) {
    stop('orders must be a list of length 3 list(pj, qj, dj)')
  }
  burn = ifelse(is.null(burn),0.1*niter,burn)
  pj = orders$pj
  qj = orders$qj
  dj = orders$dj
  Yt = t(Yt)
  k = nrow(Yt)
  N = ncol(Yt)
  if (l == 1) {
    if (is.null(Ut)) {
      Ut = matrix(0, ncol = N,nrow = 1)
    }else{
      Ut = rbind(matrix(0, ncol = N,nrow = 1),t(Ut)) # only for covariable
    }
  }else{
    if (is.null(Ut)) {
      stop('Ut must be enter with threshold process (and optional covariate process)')
    }else{
      Ut = t(Ut)
    }
  }
  nu = nrow(Ut) - 1
  Zt = Ut[1,]
  if (nu == 0) {
    Xt = matrix(0,ncol = N,nrow = 1)
    qj = rep(0,l)
  }else{
    Xt = matrix(Ut[-1,],nrow = nu,ncol = N,byrow = TRUE)
  }
  eta = 1 + pj*k + qj*nu + dj
  # Checking
  if (!is.numeric(Yt)) {stop('Yt must be a real matrix of dimension Nxk')}
  if (l != 1) {
    if (!is.numeric(Ut)) {stop('Ut must be a real matrix of dimension Nx(nu+1)')}
    if (ncol(Ut) != ncol(Yt)) {stop("Ut and Yt number of rows must match")}
  }
  if (!is.null(r)) {
    if (!length(r) == {l - 1}) {stop('r must have length l-1')}
    rini = NULL
    if (l == 1) {stop('One regime must not have threshold')}
  }
  if (is.vector(pj) & is.vector(qj) & is.vector(dj)) {
    if (!{length(pj) == l & length(qj) == l & length(dj) == l}) {
      stop('pj qj and dj must have length l')
    }else{
      for (lj in 1:l) {
        if (!{round(pj[lj]) == pj[lj] & pj[lj] >= 0}) {stop('pj[j] must be a positive integer or 0')}
        if (!{round(qj[lj]) == qj[lj] & qj[lj] >= 0}) {stop('qj[j] must be a positive integer or 0')}
        if (!{round(dj[lj]) == dj[lj] & dj[lj] >= 0}) {stop('dj[j] must be a positive integer or 0')}
      }
    }
  }else{stop('pj qj and dj must be of type vector')}
  if (!is.null(Sigma)) {
    if (!is.null(sigmaini)) {stop('Sigma known not need initial parameters')}
    sigmaini = NULL
    if (!is.list(Sigma)) {
      stop('Sigma must be a list of length l of real positive matrix of dimension kxk')
    }else{
      if (length(Sigma) != l) {
        stop('Sigma must be a list of length l of real positive matrix of dimension kxk')
      }else{
        for (i in 1:l) {
          if (is.numeric(Sigma[[i]])) {
            if (!is.matrix(Sigma[[i]])) {stop('Sigma[[i]] must be a matrix type object')}
            vl = sum(dim(Sigma[[i]]) == c(k,k))
            if (vl != 2) {stop('Sigma[[i]] must be a real positive matrix of dimension kxk')}
            vl = sum(eigen(Sigma[[i]])$values >= 0)
            if (vl != k) {stop('Sigma[[i]] must be a real positive matrix of dimension kxk')}
          }else{stop('Sigma[[i]] must be a real positive matrix of dimension kxk')}
        }
      }
    }
  }
  if (!is.null(thetaini)) {
    if (!is.list(thetaini)) {
      stop('thetaini must be a list type object')
    }else{
      if (length(thetaini) == l) {
        for (i in 1:l) {
          if (!is.list(thetaini[[i]])) {
            stop('thetaini[[i]] must be a list type object')
          }else{
            if (is.numeric(thetaini[[i]][[1]])) {
              if (!is.matrix(thetaini[[i]][[1]])) {stop('thetaini[[i]] first value must be a matrix type object')}
              vl = sum(dim(thetaini[[i]][[1]]) == c(k*eta[i],1))
              if (vl != 2) {stop('thetaini[[i]] first value must be a matrix of dimension (k*eta[i])x1')}
            }else{stop('thetaini[[i]] first value must be a real positive matrix of dimension (k*eta[i])x1')}
            if (is.numeric(thetaini[[i]][[2]])) {
              if (!is.matrix(thetaini[[i]][[2]])) {stop('thetaini[[i]] second value must be a matrix type object')}
              vl = sum(dim(thetaini[[i]][[2]]) == c(k*eta[i],k*eta[i]))
              if (vl != 2) {stop('thetaini[[i]] second value must be a matrix of dimension (k*eta[i])xk*eta[i]')}
              vl = sum(eigen(thetaini[[i]][[2]])$values >= 0)
              if (vl != (k*eta[i])) {stop('thetaini[[i]] second value must be a real positive matrix of dimension k*eta[i]xk*eta[i]')}
            }
          }
        }
      }else{stop('thetaini must be a list of length l')}
    }
  }
  if (!is.null(sigmaini)) {
    if (!is.list(sigmaini)) {
      stop('sigmaini must be a list type object')
    }else{
      if (length(sigmaini) == l) {
        for (i in 1:l) {
          if (!is.list(sigmaini[[i]])) {
            stop('sigmaini[[i]] must be a list type object')
          }else{
            if (is.numeric(sigmaini[[i]][[1]])) {
              if (!is.matrix(sigmaini[[i]][[1]])) {stop('sigmaini[[i]] must be a list with first value a matrix type object')}
              vl = sum(dim(sigmaini[[i]][[1]]) == c(k,k))
              if (vl != 2) {stop('sigmaini[[i]] must be a list with first value a real positive matrix of dimension kxk')}
              vl = sum(eigen(sigmaini[[i]][[1]])$values >= 0)
              if (vl != k) {stop('sigmaini[[i]] must be a list with first value a real positive matrix of dimension kxk')}
            }else{stop('sigmaini[[i]] must be a list with first value a real positive matrix of dimension kxk')}
            if (is.numeric(sigmaini[[i]][[2]]) & length(sigmaini[[i]][[2]]) == 1) {
              if (!{round(sigmaini[[i]][[2]]) == sigmaini[[i]][[2]] & sigmaini[[i]][[2]] >= k}) {
                stop('sigmaini[[i]] must be a list with second value an integer greater or equal k')
              }
            }else{stop('sigmaini[[i]] must be a list with second value an integer greater or equal k')}
          }
        }
      }else{stop('sigmaini must be a list of length l')}
    }
  }
  if (is.null(r)) {
    if (!is.list(rini)) {
      stop('rini must be a list type object with names za,zb,kappa,ini')
    }else{
      if (!all(names(rini) %in% c('za','zb','kappa','ini'))) {
        stop('rini must be a list type object with names za,zb,kappa,ini')
      }
      if (!is.null(rini$za)) {
        if (is.numeric(rini$za) & length(rini$za) == 1) {
          if (!{rini$za > 0 & rini$za < 1}) {
            stop('rini$za must be between 0 and 1, note: suggestion not less than 0.2')
          }
        }else{stop('rini$za must be a real number')}
      }
      if (!is.null(rini$zb)) {
        if (is.numeric(rini$zb) & length(rini$zb) == 1) {
          if (!{rini$zb > 0 & rini$zb < 1}) {
            stop('rini$zb must be between 0 and 1, note: suggestion not less than 0.8')
          }
        }else{stop('rini$zb must be a real number')}
      }
      if (!is.null(rini$kappa)) {
        if (is.numeric(rini$kappa) & length(rini$kappa) == 1) {
          #if (!{rini$kappa >= 0.5 & rini$kappa <= 0.8}) {
          #  stop('rini$kappa suggest between 0.5 and 0.8')
          #}
        }else{stop('rini$kappa must be a real number')}
      }
    }
  }
  # functions and values for r
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
  # initials values for r
  a = ifelse(is.null(rini$za),min(Zt),quantile(Zt,probs = rini$za))
  b = ifelse(is.null(rini$zb),max(Zt),quantile(Zt,probs = rini$zb))
  # function for product of Brobdingnag numbers
  prodB = function(x){
    prod = 1
    for (a in 1:length(x)) {
      prod = prod*x[a]
    }
    return(prod)
  }
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
    Ind = c()
    for (j in 1:l) {
      for (w in 1:N) {
        if (Zt[w] > rj[1,j] & Zt[w] <= rj[2,j]) {
          Ind[w] = j
        }
      }
    }
    Nrg = c()
    listaWj = listaYj = vector("list", l)
    for (lj in 1:l) {
      p = pj[lj]
      q = qj[lj]
      d = dj[lj]
      maxj = max(p,q,d)
      Inj = which(Ind == lj)
      Inj = Inj[Inj > maxj]
      Nrg[lj] = length(Inj)
      Yj = matrix(Yt[,Inj],nrow = k,ncol = Nrg[lj])
      # matrix Wj =(1,lagY,lagX,lagZ)
      Wj = matrix(0,nrow = eta[lj],ncol = Nrg[lj])
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
    return(list(Nrg = Nrg,listaW = listaWj,listaY = listaYj))
  }
  fycond = function(i2,listr,...){
    acum = 0
    Nrg = listr$Nrg
    for (lj in 1:l) {
      yj = c(listr$listaY[[lj]])
      Wj = listr$listaW[[lj]]
    if (is.null(Sigma)) {
      acum = acum + t(yj - {t(Wj) %x% diag(k)} %*% theta_iter[[lj]][,i2]) %*% {
        diag(Nrg[lj]) %x% solve(sigma_iter[[lj]][[i2]])} %*% (yj - {t(Wj) %x% diag(k)} %*% theta_iter[[lj]][,i2])
    }else{
      acum = acum + t(yj - {t(Wj) %x% diag(k)} %*% theta_iter[[lj]][,i2]) %*% {
        diag(Nrg[lj]) %x% solve(sigma[[lj]])} %*% (yj - {t(Wj) %x% diag(k)} %*% theta_iter[[lj]][,i2])
    }
  }
  if (is.null(Sigma)) {
    sigmareg = lapply(sigma_iter,function(x){x[[i2]]})
    val = prodB(Brobdingnag::as.brob(sapply(sigmareg,function(x){
      return(c(determinant(x,logarithm = FALSE)$modulus))}))^{-Nrg/2})*exp(-1/2*Brobdingnag::as.brob(acum))
  }else{
    val = prodB(Brobdingnag::as.brob(sapply(sigma,function(x){
      return(c(determinant(x,logarithm = FALSE)$modulus))}))^{-Nrg/2})*exp(-1/2*Brobdingnag::as.brob(acum))
  }
  return(val)
}
  #objects for each regimen and iterations
  theta_iter = itheta0j = isigma0j = vector("list", l)
  if (is.null(Sigma)) {
    sigma_iter = iS0j = inu0j = vector("list", l)
  }else{
    sigma = vector("list", l)
  }
  #set initial values for each regime in each chain
  for (lj in 1:l) {
    theta_iter[[lj]] = matrix(ncol = niter + burn,nrow = k*eta[lj])
    if (is.null(thetaini)) {
      itheta0j[[lj]] = rep(0,k*eta[lj])
      isigma0j[[lj]] = diag(k*eta[lj])
    }else{
      itheta0j[[lj]] = thetaini[[lj]][[1]]
      isigma0j[[lj]] = thetaini[[lj]][[2]]
    }
    theta_iter[[lj]][,1] = mvtnorm::rmvnorm(1,mean = itheta0j[[lj]],sigma = isigma0j[[lj]])
    if (is.null(Sigma)) {
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
    }else{
      sigma[[lj]] = Sigma[[lj]]
    }
  }
  # for r unknown
  if (is.null(r) & l != 1) {
    r_iter = matrix(ncol = niter + burn,nrow = l - 1)
    # use factor kappa to scale covariance matrix for proposal
    if (is.null(rini$kappa)) {
      rsig = 0.5*diag(l - 1)
    }else{
      rsig = rini$kappa*diag(l - 1)
    }
    #last check
    if (!is.null(rini$ini)) {
      if (is.numeric(rini$ini) & length(rini$ini) == {l - 1}) {
        if (dmunif(rini$ini,a,b) == 0) {
          stop('rini$ini must be in Zt range and for l > 3, r[i} < r[i+1]')
        }
      }else{stop('rini$ini must be vector of length l - 1')}
    }
    if (is.null(rini$ini)) {
      r_iter[,1] = c(quantile(Zt, probs = 1/l*(1:{l - 1})))
    }else{
      r_iter[,1] = rini$ini
    }
  }else if (is.null(r) & l == 1) {r = 0}
  # iterations gibbs and metropolis for r unknown
  if (is.null(r)) {
    pb = txtProgressBar(min = 2, max = niter + burn, style = 3)
    acep = 0
    for (i in 2:{niter + burn}) {
      listj = lists(r_iter[,i - 1])
      for (lj in 1:l) {
        Wj = listj$listaW[[lj]] 
        Yj = listj$listaY[[lj]]
        Nj = listj$Nrg[lj]
        yj = c(Yj)
        theta0j = itheta0j[[lj]]
        sigma0j = isigma0j[[lj]]
        if (!is.null(Sigma)) {
          Vj = solve(Wj %*% t(Wj) %x% solve(sigma[[lj]] %*% sigma[[lj]]) + solve(sigma0j))
          thetaj = Vj %*% {(Wj %x% solve(sigma[[lj]] %*% sigma[[lj]])) %*% yj + solve(sigma0j) %*% theta0j}
          theta_iter[[lj]][,i] = mvtnorm::rmvnorm(1,mean = thetaj,sigma = Vj)
        }else{
          S0j = iS0j[[lj]]
          nu0j = inu0j[[lj]]
          Vj = solve(Wj %*% t(Wj) %x% solve(sigma_iter[[lj]][[i - 1]]) + solve(sigma0j))
          thetaj = Vj %*% {(Wj %x% solve(sigma_iter[[lj]][[i - 1]])) %*% yj + solve(sigma0j) %*% theta0j}
          theta_iter[[lj]][,i] = mvtnorm::rmvnorm(1,mean = thetaj,sigma = Vj)
          Hj = ks::invvec(theta_iter[[lj]][,i],nrow = k,ncol = eta[lj])
          Sj = (Yj - Hj %*% Wj) %*% t(Yj - Hj %*% Wj)
          sigma_iter[[lj]][[i]] = LaplacesDemon::rinvwishart(nu = Nj + nu0j,S = Sj + S0j)
        }
      }
      # Use of metropolis with random walk
      #ek = mvtnorm::rmvnorm(1,mean = rep(0,l - 1),sigma = rsig)
      ek = runif(l - 1,-0.0375,0.0375)
      rk = r_iter[,i - 1] + ek
      listrk = lists(rk)
      pr = dmunif(rk,a,b)*fycond(i,listrk)
      px = dmunif(r_iter[,i - 1],a,b)*fycond(i,listj)
      alpha = min(1,as.numeric(pr/px))
      if (alpha >= runif(1)) {
        r_iter[,i] = rk
        acep = acep + 1
      }else{
        r_iter[,i] = r_iter[,i - 1]
        acep = acep
      }
      setTxtProgressBar(pb,i)
    }
    close(pb)
  }else{#r known
    listj = lists(r)
    for (lj in 1:l) {
      Wj = listj$listaW[[lj]] 
      Yj = listj$listaY[[lj]]
      Nj = listj$Nrg[lj]
      yj = c(Yj)
      theta0j = itheta0j[[lj]]
      sigma0j = isigma0j[[lj]]
      pb = txtProgressBar(min = 2, max = niter + burn, style = 3)
      for (i in 2:{niter + burn}) {
        if (!is.null(Sigma)) {
          Vj = solve(Wj %*% t(Wj) %x% solve(sigma[[lj]] %*% sigma[[lj]]) + solve(sigma0j))
          thetaj = Vj %*% {(Wj %x% solve(sigma[[lj]] %*% sigma[[lj]])) %*% yj + solve(sigma0j) %*% theta0j}
          theta_iter[[lj]][,i] = mvtnorm::rmvnorm(1,mean = thetaj,sigma = Vj)
        }else{
          S0j = iS0j[[lj]]
          nu0j = inu0j[[lj]]
          Vj = solve(Wj %*% t(Wj) %x% solve(sigma_iter[[lj]][[i - 1]]) + solve(sigma0j))
          thetaj = Vj %*% {(Wj %x% solve(sigma_iter[[lj]][[i - 1]])) %*% yj + solve(sigma0j) %*% theta0j}
          theta_iter[[lj]][,i] = mvtnorm::rmvnorm(1,mean = thetaj,sigma = Vj)
          Hj = ks::invvec(theta_iter[[lj]][,i],nrow = k,ncol = eta[lj])
          Sj = (Yj - Hj %*% Wj) %*% t(Yj - Hj %*% Wj)
          sigma_iter[[lj]][[i]] = LaplacesDemon::rinvwishart(nu = Nj + nu0j,S = Sj + S0j)
        }
        setTxtProgressBar(pb,i)
      }
    }
    close(pb)
  }
  # objects for chains and info in each regime
  Rest = thetaest = thetachain = vector("list", l)
  names(Rest) = names(thetaest) = names(thetachain) = paste0('R',1:l)
  if (is.null(Sigma)) {
    sigmaest = sigmachain = vector("list", l)
    names(sigmaest) = names(sigmachain) = paste0('R',1:l)
  }
  # save chains and creation of the 'regime' type object
  for (lj in 1:l) {
    # save chains of theta
    thetachain[[lj]] = theta_iter[[lj]][,-c(1:burn)]
    # credibility intervals for theta
    vectheta = matrix(nrow = k*eta[lj],ncol = 3)
    colnames(vectheta) = c(paste0('lower limit ',(1 - level)/2*100,'%'),'mean',paste0('upper limit ',(1 + level)/2*100,'%'))
    vectheta[,1] = apply(thetachain[[lj]],1,quantile,probs = (1 - level)/2)
    vectheta[,3] = apply(thetachain[[lj]],1,quantile,probs = (1 + level)/2)
    vectheta[,2] = apply(thetachain[[lj]],1,mean)
    thetaest[[lj]] = vectheta
    if (nu != 0 & qj[lj] != 0 & dj[lj] != 0) {
      rownames(vectheta) = 
        rep(c('phi0',rep(paste0('phi',1:pj[lj]),each = k),rep(paste0('beta',1:qj[lj]),each = nu),paste0('delta',1:dj[lj])),k)
    }else if (nu != 0 & qj[lj] != 0 & dj[lj] == 0) {
      rownames(vectheta) = 
        rep(c('phi0',rep(paste0('phi',1:pj[lj]),each = k),rep(paste0('beta',1:qj[lj]),each = nu)),k)
    }else if (qj[lj] == 0 & dj[lj] != 0) {
      rownames(vectheta) = 
        rep(c('phi0',rep(paste0('phi',1:pj[lj]),each = k),paste0('delta',1:dj[lj])),k)
    }else if (qj[lj] == 0 & dj[lj] == 0) {
      rownames(vectheta) = 
        rep(c('phi0',rep(paste0('phi',1:pj[lj]),each = k)),k)
    }
    if (is.null(Sigma)) {
      # save chains of sigma^1/2
      sigchain = matrix(nrow = k*k,ncol = (niter + burn))
      for (o in 1:(niter + burn)) {sigchain[,o] = c(expm::sqrtm(sigma_iter[[lj]][[o]]))}
      sigmachain[[lj]] = matrix(sigchain[,-c(1:burn)],ncol = niter,nrow = k*k)
      # credibility intervals for sigma^1/2
      vecsigma = matrix(nrow = k*k,ncol = 3)
      colnames(vecsigma) = c(paste0('lower limit ',(1 - level)/2*100,'%'),'mean',paste0('upper limit ',(1 + level)/2*100,'%'))
      a = paste0(1:k,1)
      for (i3 in 2:k) {a = c(a,paste0(1:k,i3))}
      rownames(vecsigma) = a
      vecsigma[,1] = apply(sigmachain[[lj]],1,quantile,probs = (1 - level)/2)
      vecsigma[,3] = apply(sigmachain[[lj]],1,quantile,probs = (1 + level)/2)
      vecsigma[,2] = apply(sigmachain[[lj]],1,mean)
      sigmaest[[lj]] = vecsigma
    }
    # creation of the 'regime' type object
    p = pj[lj]
    q = qj[lj]
    d = dj[lj]
    if (q == 0 & d == 0) {
      thetaind = c(0,(10 + (1:p)) %x% rep(1,k))
    }else if (q != 0 & d == 0) {
      thetaind = c(0,(10 + (1:p)) %x% rep(1,k),(20 + (1:q)) %x% rep(1,nu))
    }else if (q == 0 & d != 0) {
      thetaind = c(0,(10 + (1:p)) %x% rep(1,k),30 + (1:d))
    }else{
      thetaind = c(0,(10 + (1:p)) %x% rep(1,k),(20 + (1:q)) %x% rep(1,nu),30 + (1:d))
    }
    Thetaj = ks::invvec(thetaest[[lj]][,2],ncol = eta[lj],nrow = k)
    Ri = NULL
    Ri$cs = matrix(Thetaj[,thetaind == 0],nrow = k,ncol = 1)
    Ri$phi = vector("list", p)
    names(Ri$phi) = paste0('phi',1:p)
    for (j in 1:p) {
      Ri$phi[[j]] = matrix(Thetaj[,thetaind == (10 + j)],nrow = k,ncol = k)
    }
    if (q != 0) {
      Ri$beta = vector("list", q)
      names(Ri$beta) = paste0('beta',1:q)
      for (j in 1:q) {Ri$beta[[j]] = matrix(Thetaj[,thetaind == (20 + j)],nrow = k,ncol = nu)}
    }
    if (d != 0) {
      Ri$delta = vector("list", d)
      names(Ri$delta) = paste0('delta',1:d)
      for (j in 1:d) {Ri$delta[[j]] = matrix(Thetaj[,thetaind == (30 + j)],nrow = k,ncol = 1)}
    }
    if (is.null(Sigma)) {
      Ri$sigma = ks::invvec(sigmaest[[lj]][,2],ncol = k,nrow = k)
    }else{
      Ri$sigma = sigma[[lj]]
    }
    class(Ri) = 'regime'
    Rest[[lj]] = Ri
  }
  if (is.null(r)) {
    rest = matrix(nrow = l - 1,ncol = 3)
    colnames(rest) = colnames(rest) = 
      c(paste('lower limit ',(1 - level)/2*100,'%',sep = ""),'mean',paste('upper limit ',(1 + level)/2*100,'%',sep = ""))
    rchain = matrix(r_iter[,-c(1:burn)],ncol = niter,nrow = l - 1)
    rest[,1] = apply(rchain,1,quantile,probs = (1 - level)/2)
    rest[,3] = apply(rchain,1,quantile,probs = (1 + level)/2)
    rest[,2] = apply(rchain,1,mean)
  }
  # logLik
  if (is.null(r)) {rt = rest[,2]}else {rt = r}
  listj = lists(rt)
  logLikj = c()
  for (lj in 1:l) {
    Wj = listj$listaW[[lj]]
    Yj = listj$listaY[[lj]]
    Nj = listj$Nrg[lj]
    yj = c(Yj)
    Hj = ks::invvec(thetaest[[lj]][,2],nrow = k,ncol = eta[lj])
    Sj = (Yj - Hj %*% Wj) %*% t(Yj - Hj %*% Wj)
    logLikj[lj] = log(det(Sj/Nj))
  }
  # exits
  ## Chain
  if (chain) {
    Chain = NULL
    Chain$Theta = thetachain
    if (is.null(r) & l != 1) {Chain$r = rchain}
    if (is.null(Sigma)) {Chain$Sigma = sigmachain}
  }
  ## estimates and credibility interval
  estimates = NULL
  estimates$Theta = thetaest
  if (is.null(r) & l != 1) {estimates$r = rest}
  if (is.null(Sigma)) {estimates$Sigma = sigmaest}
  data = NULL
  data$yt = t(Yt)
  data$Ut = t(Ut)
  if (chain) {
    results = list(Nj = listj$Nrg,estimates = estimates,regime = Rest,Chain = Chain,logLikj = logLikj,data = data,r = r,orders = list(pj = pj,qj = qj,dj = dj))
    class(results) = "regime-model" 
  }else{
    results = list(Nj = listj$Nrg,estimates = estimates,regime = Rest,logLikj = logLikj,data = data,r = r,orders = list(pj = pj,qj = qj,dj = dj))
    class(results) = "regime-model"
  }
  return(results)
}
# Example data(yt,ut):
yt = datasim #(Buscar datos)
#r known
estim1 = mtarns(Yt = yt$Yt,Ut = Ut,l = 2,orders = list(pj = yt$pj,dj = yt$dj,qj = yt$qj),r = qnorm(0.4),niter = 500,chain = TRUE)
forecast::autoplot(ts(t(estim1$Chain$Theta$R1)),facets = TRUE)
forecast::autoplot(ts(t(estim1$Chain$Theta$R2)),facets = TRUE)
#r unknown
estim2 = mtarns(Yt = yt$Yt,Ut = Ut,l = 2,orders = list(pj = yt$pj,dj = yt$dj,qj = yt$qj),niter = 500,chain = TRUE)
forecast::autoplot(ts(t(estim2$Chain$Theta$R1)),facets = TRUE)
forecast::autoplot(ts(t(estim2$Chain$Sigma$R1)),facets = TRUE)
forecast::autoplot(ts(t(estim2$Chain$Theta$R2)),facets = TRUE)
forecast::autoplot(ts(t(estim2$Chain$Sigma$R2)),facets = TRUE)
forecast::autoplot(ts(t(estim2$Chain$r)))