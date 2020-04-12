autoplot.mtarsim = function(object, type = 1 , ...) {
  if (!requireNamespace('ggplot2', quietly = TRUE)) {
    stop('ggplot2 is needed for this function to work. Install it via install.packages(\'ggplot2\')', call. = FALSE)
  }
  else {
    if (!inherits(object, 'mtarsim')) {
      stop('autoplot.mtarsim requires a mtarsim object, use object=object')
    }}
  dats_Yt = as.data.frame(object$Yt)
  time = seq(1,nrow(dats_Yt))
  dat = data.frame(name = 'Series.1',time = time,value = dats_Yt[,1])
  for (i in 2:ncol(object$Yt)) {
    dat = rbind(dat,data.frame(name = paste0('Series.',i),time = time,value = dats_Yt[,i]))
  }
  p = ggplot2::ggplot(data = dat,ggplot2::aes(x = time,y = value))
  p = p + ggplot2::geom_line() + ggplot2::facet_grid(name~.) + theme_bw()
  p = p +  labs(title = "Output processes")
  if (!is.null(object$Zt)) {
    dats_Zt = data.frame(time = time,value = object$Zt)
    p2 = ggplot(data = dats_Zt,aes(x = time,y = value))
    p2 = p2 + geom_line() + theme_bw()
    Nrg_plot = paste0(paste0(paste0("r",1:object$l),"="),object$Summary_r$Prop_reg,"%")
    p2 = p2 + labs(title = "Threshold processes",
                   subtitle = paste0("(",paste(Nrg_plot,collapse = ","),")"))
    if (!is.null(object$r)) {
      for (i in c(object$r)) {
        p2 = p2 + geom_hline(yintercept = i,linetype = "dashed",color = "blue")
      }
    }
  }
  if (!is.null(object$Xt)) {
    dats_Xt = as.data.frame(object$Xt)
    dat2 = data.frame(name = 'Series.1',time = time,value = dats_Xt[,1])
    if (ncol(dats_Xt) > 1) {
      for (i in 2:ncol(object$Xt)) {
        dat2 = rbind(dat2,data.frame(name = paste0('Series.',i),
                                     time = time,value = dats_Xt[,i]))
      }
    }
    p3 = ggplot2::ggplot(data = dat2,ggplot2::aes(x = time,y = value))
    p3 = p3 + ggplot2::geom_line() + ggplot2::facet_grid(name~.) + theme_bw()
    p3 = p3 +  labs(title = "Covariates processes")
  }
  if (type == 1) {
    return(p)
  }
  if (type == 2) {
    if (!is.null(object$Zt)) {
      stop("Threshold processes does not exist")}
    return(p2)
  }
  if (type == 3) {
    if (!is.null(object$Xt)) {
      stop("Covariates processes does not exist")}
    return(p3)
  }
}