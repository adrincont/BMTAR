# Date:
# Description:
# Coments:
# Function:
mtarmissing = function(Yt, Ut = NULL, l, orders, level = 0.95, method= c('KUO','SSVS')[1], niter, burn = NULL,
                       chain = FALSE, thetaini = NULL, sigmaini = NULL, gammaini = NULL, 
                       rini = list(za = NULL, zb = NULL, kappa = NULL, ini = NULL),
                       Cij = NULL, Tauij = NULL, R = NULL){
  if (!is.matrix(Yt)) {stop('Yt must be a matrix type object')}
  if (l != 1 & is.null(Ut)) {stop('Threshold process Ut (optional covariable) must be enter')}
  if (l != 1 & !is.matrix(Ut)) {stop('Ut must be a matrix type object')}
  if (!is.list(orders) | length(orders) != 3) {
    stop('orders must be a list of length 3 list(pjmax, qjmax, djmax)')
  }else if (sum(names(orders) %in% c('pjmax','qjmax','djmax')) != 3) {
    stop('orders must be a list of length 3 list(pjmax, qjmax, djmax)')
  }
  burn = ifelse(is.null(burn),round(0.1*niter),burn)
  pjmax = orders$pjmax
  qjmax = orders$qjmax
  djmax = orders$djmax
  Yt = t(Yt)
  k = nrow(Yt)
  N = ncol(Yt)
  # Indicadores de datos faltantes para Yt
  PosNAMat = PosNAvec = PosNAvecT = vector(mode = 'list',2)
  PosNAMat[[1]] = apply(Yt,2,is.na)
  PosNAvec[[1]] = c(1:ncol(Yt))[apply(PosNAMat[[1]],2,any)]
  PosNAvecT[[1]] = matrix(rep(c(1:Tlen),k),nrow = k,ncol = Tlen,byrow = T)[PosNAMat[[1]]]
  # -------------------------------------
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
  # Indicadores de datos faltantes para Ut
  
  Zt = Ut[1,]
  if (nu == 0) {
    Xt = matrix(0,ncol = N,nrow = 1)
    qjmax = rep(0,l)
  }else{
    Xt = matrix(Ut[-1,],nrow = nu,ncol = N,byrow = TRUE)
  }
  eta = 1 + pjmax*k + qjmax*nu + djmax
  #checking
  if (!(length(pjmax) == l & length(qjmax) == l & length(djmax) == l)) {stop('pjmax qjmax and djmax must be of length l')}
  for (lj in 1:l) {
    if (pjmax[lj] > 5) {stop('pjmax must be smaller or 5 for each regime')}
    if (qjmax[lj] > 5) {stop('qjmax must be smaller or 5 for each regime')}
    if (djmax[lj] > 5) {stop('djmax must be smaller or 5 for each regime')}
  }
  if (!is.numeric(Yt)) {stop('Yt must be a real matrix of dimension Nxk')}
  if (l != 1) {
    if (!is.numeric(Ut)) {stop('Ut must be a real matrix of dimension Nx(nu+1)')}
    if (ncol(Ut) != ncol(Yt)) {stop('Ut and Yt number of rows must match')}
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
  if (!is.null(gammaini)) {
    if (!is.list(gammaini)) {stop('gammaini must be a list type object')}
    if (length(gammaini) != l) {stop('gammaini must be a list of lenght l')
    }else{
      for (lj in 1:l) {
        if (length(gammaini[[lj]]) != k*eta[lj]) {stop('gammaini[[j]] must be a vector of lenght k*etaj')}
        if (sum(gammaini[[lj]] >= 0 & gammaini[[lj]] <= 1) != k*eta[lj]) {stop('gammaini[[j]] values must be between 0 and 1')}
      } 
    }
  }
  if (!is.list(rini)) {
    stop('rini must be a list type object with names za,zb,kappa,ini')
  }else{
    if (sum(names(rini) %in% c('za','zb','kappa','ini')) < 1) {
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
        if (!{rini$kappa <= 0.8}) {
          stop('rini$kappa suggest between 0.5 and 0.8')
        }
      }else{stop('rini$kappa must be a real number')}
    }
  }
  if (method == 'SSVS') {
    if (!is.null(Cij)) {
      if (!is.list(Cij)) {stop('Cij must be a list object type')}
      if (length(Cij) != l) {stop('Cij must be a list of length l')
      }else{
        for (lj in 1:l) {
          if (length(Cij[[lj]]) != k*eta[lj]) {stop('Cij[[j]] must be a vector of length k*etaj')}
        } 
      }
    }
    if (!is.null(Tauij)) {
      if (!is.list(Tauij)) {stop('Tauij must be a list object type')}
      if (length(Tauij) != l) {stop('Tauij must be a list of length l')
      }else{
        for (lj in 1:l) {
          if (length(Tauij[[lj]]) != k*eta[lj]) {stop('Tauij[[j]] must be a vector of length k*etaj')}
        } 
      }
    }
    if (!is.null(R)) {
      if (!is.list(R)) {stop('R must be a list object type')}
      if (length(R) != l) {stop('R must be a list of length l')
      }else{
        for (lj in 1:l) {
          if (is.numeric(R[[lj]])) {
            if (!is.matrix(R[[lj]])) {stop('R[[j]] must be a list with a matrix type object')}
            if (sum(dim(R[[lj]]) == c(k*eta[lj],k*eta[lj])) != 2) {stop('R[[j]] must be a matrix of dimension k*etaj')}
            vl = sum(eigen(R[[lj]])$values >= 0)
            if (vl != k*eta[lj]) {stop('R[[j]] must be a real positive matrix of dimension k*etaj')}
          }else{stop('R must be a list with l real positive matrix of dimension k*etaj')}
        }
      }
    }
  }
  #
  #objects for each regimen and iterations
  theta_iter = sigma_iter = gam_iter = pij =  Dj = Rj = vector('list', l)
  if (method == 'SSVS') {
    tauij = itauij = cij = vector('list', l)
  }
  if (l != 1) {r_iter = matrix(ncol = niter + burn,nrow = l - 1)}
  #set initial values for each regime in each chain
  itheta0j = isigma0j = iS0j = inu0j = vector('list',l)
  for (lj in 1:l) {
    theta_iter[[lj]] = gam_iter[[lj]] = matrix(ncol = niter + burn,nrow = k*eta[lj])
    sigma_iter[[lj]] = list()
    length(sigma_iter[[lj]]) = niter + burn
    if (is.null(thetaini)) {
      itheta0j[[lj]] = rep(0,k*eta[lj])
      isigma0j[[lj]] = diag(k*eta[lj])
    }else{
      itheta0j[[lj]] = thetaini[[lj]][[1]]
      isigma0j[[lj]] = thetaini[[lj]][[2]]
    }
    if (is.null(sigmaini)) {
      iS0j[[lj]] = diag(k)
      inu0j[[lj]] = k
    }else{
      iS0j[[lj]] = sigmaini[[lj]][[1]]
      inu0j[[lj]] = sigmaini[[lj]][[2]]
    }
    if (is.null(gammaini)) {
      pij[[lj]] = rep(0.5,k*eta[lj])
    }else{
      pij[[lj]] = gammaini[[lj]]
    }
    gam_iter[[lj]][,1] = rep(1,k*eta[lj])
    if (method == 'KUO') {
      theta_iter[[lj]][,1] = mvtnorm::rmvnorm(1,mean = itheta0j[[lj]],sigma = isigma0j[[lj]])
    }
    if (method == 'SSVS') {
      if (is.null(Cij)) {
        cij[[lj]] = rep(25,k*eta[lj])
      }else{cij[[lj]] = Cij[[lj]]}
      if (is.null(Tauij)) {
        if (l == 2) {
          tauij[[lj]] = rep(1.25,k*eta[lj])
        }else{tauij[[lj]] = rep(1.5,k*eta[lj])}
      }else{tauij[[lj]] = Tauij[[lj]]}
      if (is.null(R)) {
        Rj[[lj]] = diag(k*eta[lj])
      }else{Rj[[lj]] = R[[lj]]}
      itauij[[lj]] = cij[[lj]]*tauij[[lj]]
      Dj[[lj]] = diag(itauij[[lj]])
      theta_iter[[lj]][,1] = mvtnorm::rmvnorm(1,mean = rep(0,k*eta[lj]),sigma = Dj[[lj]] %*% Rj[[lj]] %*% Dj[[lj]])
    }
    sigma_iter[[lj]][[1]] = LaplacesDemon::rinvwishart(nu = inu0j[[lj]],S = iS0j[[lj]])
  }
  if (is.null(rini$kappa)) {
    rsig = 0.5*diag(l - 1)
  }else{
    rsig = rini$kappa*diag(l - 1)
  }
  #
  #objects for save regimes, chains and credibility intervals
  Rest = thetaest = thetachain = 
    gamest = gamchain = sigmaest = sigmachain = vector('list', l)
  names(Rest) = names(thetaest) = names(thetachain) = names(gamest) = names(gamchain) = 
    names(sigmaest) = names(sigmachain) = paste0('R',1:l)
  #
  #necesary functions
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
    listaXj = listaYj = listaWj = vector('list', l)
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
      Xj = matrix(1,nrow = k*Nrg[lj],ncol = k*eta[lj])
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
  fycond = function(i2,listr,gamma,...){
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
  rgamber = function(pos,reg,i,...){
    gam_j = gam_iter
    gam_j[[reg]][pos,i] = 1
    pycond1 = fycond(i,listj,gam_j)
    gam_j[[reg]][pos,i] = 0
    pycond0 = fycond(i,listj,gam_j)
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
  #
  #edges for rmunif
  a = ifelse(is.null(rini$za),quantile(Zt,probs = 0.2),quantile(Zt,probs = rini$za))
  b = ifelse(is.null(rini$zb),quantile(Zt,probs = 0.8),quantile(Zt,probs = rini$zb))
  #
  #last check
  if (!is.null(rini$ini)) {
    if (is.numeric(rini$ini) & length(rini$ini) == {l - 1}) {
      if (dmunif(rini$ini,a,b) == 0) {
        stop('rini$ini must be in Zt range and for l > 3, r[i] < r[i+1]')
      }
    }else{stop('rini$ini must be vector of length l - 1')}
  }
  if (l != 1) {
    if (is.null(rini$ini)) {
      r_iter[,1] = c(quantile(Zt, probs = 1/l*(1:{l - 1})))
    }else{
      r_iter[,1] = rini$ini
    }
  }
  #iterations: gibbs and metropolis sampling
  acep = 0
  pb = txtProgressBar(min = 2, max = niter + burn,style = 3)
  for (i in 2:{niter + burn}) {
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
      sigma_iter[[lj]][[i]] = LaplacesDemon::rinvwishart(nu = Nj + nu0j,S = Sj + S0j)
      gam_iter[[lj]][,i] = gam_iter[[lj]][,i - 1]
    }
    for (jj in 1:l) {
      for (ii in 1:{k*eta[jj]}) {
        gam_iter[[jj]][ii,i] = rgamber(pos = ii,reg = jj,i = i)
      }
    }
    if (l != 1) {
      #ek = mvtnorm::rmvnorm(1,mean = rep(0,l - 1),sigma = rsig)
      ek = runif(l - 1,-0.00375,0.00375)
      rk = r_iter[,i - 1] + ek
      listrk = lists(rk)
      pr = dmunif(rk,a,b)*fycond(i,listrk,gam_iter)
      px = dmunif(r_iter[,i - 1],a,b)*fycond(i,listj,gam_iter)
      alpha = min(1,as.numeric(pr/px))
      if (alpha >= runif(1)) {
        r_iter[,i] = rk
        acep = acep + 1
      }else{
        r_iter[,i] = r_iter[,i - 1]
        acep = acep
      }
    }
    setTxtProgressBar(pb,i)
  }
  close(pb)
  #exits 
  if (l != 1) {
    rest = matrix(nrow = l - 1,ncol = 3)
    colnames(rest) = colnames(rest) = 
      c(paste('lower limit ',(1 - level)/2*100,'%',sep = ''),'mean',paste('upper limit ',(1 + level)/2*100,'%',sep = ''))
    rchain = matrix(r_iter[,-c(1:burn)],ncol = niter,nrow = l - 1)
    rest[,1] = apply(rchain,1,quantile,probs = (1 - level)/2)
    rest[,3] = apply(rchain,1,quantile,probs = (1 + level)/2)
    rest[,2] = apply(rchain,1,mean)
  }
  #save chains
  for (lj in 1:l) {
    thetachain[[lj]] = theta_iter[[lj]][,-c(1:burn)]
    gamchain[[lj]] = gam_iter[[lj]][,-c(1:burn)]
    sigchain = matrix(nrow = k*k,ncol = (niter + burn))
    for (o in 1:(niter + burn)) {sigchain[,o] = c(expm::sqrtm(sigma_iter[[lj]][[o]]))}
    sigmachain[[lj]] = matrix(sigchain[,-c(1:burn)],ncol = niter,nrow = k*k)
  }
  for (lj in 1:l) {
    #credibility intervals
    vecgamma = vectheta = matrix(nrow = k*eta[lj],ncol = 3)
    vecsigma = matrix(nrow = k*k,ncol = 3)
    colnames(vectheta) = colnames(vecsigma) = 
      c(paste0('lower limit ',(1 - level)/2*100,'%'),'mean',paste0('upper limit ',(1 + level)/2*100,'%'))
    a = paste0(1:k,1)
    for (i3 in 2:k) {a = c(a,paste0(1:k,i3))}
    rownames(vecsigma) = a
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
    vectheta[,1] = apply(thetachain[[lj]],1,quantile,probs = (1 - level)/2)
    vectheta[,3] = apply(thetachain[[lj]],1,quantile,probs = (1 + level)/2)
    vectheta[,2] = apply(thetachain[[lj]],1,mean)
    thetaest[[lj]] = vectheta
    vecsigma[,1] = apply(sigmachain[[lj]],1,quantile,probs = (1 - level)/2)
    vecsigma[,3] = apply(sigmachain[[lj]],1,quantile,probs = (1 + level)/2)
    vecsigma[,2] = apply(sigmachain[[lj]],1,mean)
    sigmaest[[lj]] = vecsigma
    vec = apply(gamchain[[lj]],2,paste,collapse = '')
    vecs = sort(table(vec), decreasing = TRUE)[1:3]
    colnames(vecgamma) = c(paste('primer frec',paste0(vecs[1]/niter*100,'%')),
                           paste('segundo frec',paste0(vecs[2]/niter*100,'%')),
                           paste('tercer frec',paste0(vecs[3]/niter*100,'%')))
    vecgamma[,1] = gamchain[[lj]][,which(vec == names(vecs[1]))[1]]
    vecgamma[,2] = gamchain[[lj]][,which(vec == names(vecs[2]))[1]]
    vecgamma[,3] = gamchain[[lj]][,which(vec == names(vecs[3]))[1]]
    gamest[[lj]] = vecgamma
    #creation of the 'regime' type object
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
    class(Ri) = 'regime'
    Rest[[lj]] = Ri
  }
  # logLik
  listj = lists(rest[,2])
  logLikj = c()
  for (lj in 1:l) {
    Xj = listj$listaX[[lj]]
    Yj = listj$listaY[[lj]]
    Nj = listj$Nrg[lj]
    Hj = ks::invvec({Xj %*% diag(gamest[[lj]][,1]) %*% thetaest[[lj]][,2]},nrow = k,ncol = Nj)
    Sj = (Yj - Hj) %*% t(Yj - Hj)
    logLikj[lj] = log(det(Sj/Nj))
  }
  # exits
  ## Chain
  if (chain) {
    Chain = NULL
    Chain$Theta = thetachain
    Chain$Gamma = gamchain
    Chain$Sigma = sigmachain
    if (l != 1) {Chain$r = rchain}
  }
  ## estimates
  estimates = NULL
  estimates$Theta = thetaest
  estimates$Gamma = gamest
  estimates$Sigma = sigmaest
  data = NULL
  data$yt = t(Yt)
  data$Ut = t(Ut)
  orders = NULL
  orders$pj = list(pmax = pjmax,pj = sapply(Rest,function(x){length(x$phi)}))
  orders$qj = list(qmax = qjmax,qj = sapply(Rest,function(x){length(x$beta)}))
  orders$dj = list(dmax = djmax,dj = sapply(Rest,function(x){length(x$delta)}))
  if (l != 1) {estimates$r = rest}
  if (chain) {
    results = list(Nj = listj$Nrg,estimates = estimates,regime = Rest,Chain = Chain,logLikj = logLikj,data = data,r = rest[,2],orders = orders)
    class(results) = 'regime-model' 
  }else{
    results = list(Nj = listj$Nrg,estimates = estimates,regime = Rest,logLikj = logLikj,data = data,r = rest[,2],orders = orders)
    class(results) = 'regime-model'
  }
  return(results)
}