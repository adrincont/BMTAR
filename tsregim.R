tsregim = function(Yt, Zt = NULL, Xt = NULL, r = NULL, ...){
  list_result = list()
  if (!is.null(r)) {
    if (!is.numeric(r)) {stop('r must be a numeric vector')}
    l = length(r) + 1
    if (l > 2) {for (i in 1:{l - 2}) {
      if (r[i] >= r[i + 1]) {stop('r[i] must be smaller than r[i+1]')}}
    }
    list_result$l = l
    if (is.null(Zt)) {
      stop('Zt must be enter with threshold process')
    }
  }
  if (!is.numeric(Yt)) {stop('Yt must be a real matrix of dimension Nxk')}
  if (!is.matrix(Yt)) {Yt = as.matrix(Yt)}
  if (!is.null(Zt)) {
    if (!is.numeric(Zt)) {stop('Zt must be a real matrix of dimension Nx1')}
    if (!is.matrix(Zt)) {Zt = as.matrix(Zt)}
    if (nrow(Zt) != nrow(Yt)) {stop('Zt and Yt number of rows must match')}
  }else{l = 1}
  if (!is.null(Xt)) {
    if (!is.numeric(Xt)) {stop('Xt must be a real matrix of dimension Nx(nu+1)')}
    if (!is.matrix(Xt)) {Xt = as.matrix(Xt)}
    if (nrow(Xt) != nrow(Yt)) {stop('Xt and Yt number of rows must match')}
    nu = ncol(Xt)
    list_result$nu = nu
  }
  k = ncol(Yt)
  N = nrow(Yt)
  list_result$Yt = Yt
  list_result$Zt = Zt
  if (!is.null(Xt)) {
    list_result$Xt = Xt 
  }
  list_result$r = r
  if (!is.null(r)) {
    # Calcular regimen que pertenecen las funciones
    list_result$Ind = lists_ind(r,Zt,l)
    Table_r = data.frame('N_reg' = c(table(list_result$Ind)),
                         'Prop_reg' = 100*c(prop.table(table(list_result$Ind))))
    rownames(Table_r) = paste('Regim',1:l)
    list_result$Summary_r = Table_r
  }
  list_result$N = N
  list_result$k = k
  class(list_result) = 'tsregim'
  return(list_result)
}
print.tsregim = function(x,...){
  cat('Threshold time series:\n','N =',x$N,'\n')
  dats = x
  class(dats) = NULL
  if (!is.null(x$r)) {
    cat('======================','\n')
    cat('r = ',x$r,'\n')
    print(x$Summary_r)
    cat('======================','\n')
  }else{
    if (!is.null(x$Zt)) {
      cat('Unknown threshold values','\n')
    }
  }
  str(dats)
}
lists_ind = function(r,Zt,l,...){
  N = length(Zt)
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
  return(list(Ind = Ind))
}
# Ejemplo
Yt = datasim$Yt
Ut = Ut$Yt
Zt = t(Ut)[1,]
Xt = t(Ut)[-1,]
datos = tsregim(Yt,Zt,Xt)
