#==================================================================================================#
# Date: 14/04/2020
# Description:
# Function:
#==================================================================================================#
diagnostic_mtar = function(regim_model,lagmax = NULL){
  if (!requireNamespace('ggplot2', quietly = TRUE)) {
    stop('ggplot2 is needed for this function to work. Install it via install.packages(\'ggplot2\')', call. = FALSE)
  }else {
    if (!inherits(regim_model, 'regim_model')) {
      stop('diagnostic.mtar requires a regim_model object')
    }}
  e_k = tsregim(regim_model$residuals)
  p1 = autoplot.tsregim(e_k) + ggplot2::geom_hline(yintercept = 0,color = "red") +
    ggplot2::ggtitle('Residual serie plot')
  e_data = as.data.frame(e_k$Yt)
  time = seq(1,nrow(e_data))
  dat = data.frame(label = 'Series.1',time = time,value = e_data[,1],
                   cusum = cumsum(e_data[,1])/sd(e_data[,1]),
                   cumsq = c(cumsum(e_data[,1]^2)/sum(e_data[,1]^2)))
  if (ncol(e_data) > 1) {
    for (i in 2:ncol(e_data)) {
      dat = rbind(dat,data.frame(label = paste0('Series.',i),time = time,value = e_data[,i],
                                 cusum = cumsum(e_data[,i])/sd(e_data[,i]),
                                 cumsq = c(cumsum(e_data[,i]^2)/sum(e_data[,i]^2))))
    }
  }
  p2 = ggplot2::ggplot(data = as.data.frame(dat),  ggplot2::aes(x = value, color = label )) +
    ggplot2::geom_density() + ggplot2::theme_bw()
  p2 = p2 + ggplot2::stat_function(fun = dnorm,color = "black")
  p2 = p2 + ggplot2::ggtitle("Residual density plot")

  Af = 0.948 ###Cuantil del 95% para cusum
  LS = Af*sqrt(e_k$N) + 2*Af*c(1:e_k$N)/sqrt(e_k$N)
  LI = -LS
  p3 = ggplot2::ggplot(data = dat, ggplot2::aes(x = time, y = cusum,color = label))
  p3 = p3 + ggplot2::geom_ribbon(ggplot2::aes(ymin = rep(LS,e_k$k), ymax = rep(LI,e_k$k)),
                                 fill = "gray",color = NA,alpha = 0.5)
  p3 = p3 + ggplot2::geom_line() + ggplot2::theme_bw() + ggplot2::ggtitle('CUSUM statistic for residuals')

  co = 0.13461 ####Valor del cuantil aproximado para cusumsq para T/2-1=71 alpha=0.05
  LQS = co + (1:e_k$N)/e_k$N
  LQI = -co + (1:e_k$N)/e_k$N
  p4 = ggplot2::ggplot(data = dat, ggplot2::aes(x = time, y = cumsq,color = label))
  p4 = p4 + ggplot2::geom_ribbon(ggplot2::aes(ymin = rep(LQS,e_k$k), ymax = rep(LQI,e_k$k)),
                                 fill = "gray",color = NA,alpha = 0.5)
  p4 = p4 + ggplot2::geom_line() + ggplot2::theme_bw() + ggplot2::ggtitle('CUSUMSQ statistic for residuals')

  acf_i = stats::acf(regim_model$residuals[,1],lag.max = lagmax,plot = FALSE,type = 'correlation')
  acf_Yt = data.frame(lag = acf_i$lag, value = acf_i$acf,names = 'Serie.1',type = 'ACF')
  pacf_i = stats::acf(regim_model$residuals[,1],lag.max = lagmax,plot = FALSE,type = 'partial')
  pacf_Yt = data.frame(lag = pacf_i$lag, value = pacf_i$acf,names = 'Serie.1',type = 'PACF')
  if (ncol(regim_model$residuals) > 1) {
    for (i in 2:ncol(regim_model$residuals)) {
      acf_i = stats::acf(regim_model$residuals[,i],lag.max = lagmax,plot = FALSE,type = 'correlation')
      acf_Yt = rbind(acf_Yt,data.frame(lag = acf_i$lag + 0.1*i, value = acf_i$acf,names = paste0('Series.',i),type = 'ACF'))
      pacf_i = stats::acf(regim_model$residuals[,i],lag.max = lagmax,plot = FALSE,type = 'partial')
      pacf_Yt = rbind(pacf_Yt,data.frame(lag = pacf_i$lag + 0.1*i, value = pacf_i$acf,names = paste0('Series.',i),type = 'PACF'))
    }
  }
  dat_cor = rbind(acf_Yt,pacf_Yt)
  p5 = ggplot2::ggplot(data = dat_cor[floor(dat_cor$lag) != 0,],ggplot2::aes(x = lag, y = value))
  p5 = p5 + ggplot2::geom_hline(yintercept = 0) + ggplot2::facet_grid(type~names)
  p5 = p5 + ggplot2::geom_segment(ggplot2::aes(xend = lag,yend = 0)) + ggplot2::geom_point(color = "blue",size = 0.4)
  ci = qnorm((1 + 0.95)/2)/sqrt(nrow(regim_model$residuals))
  p5 = p5 + ggplot2::geom_ribbon(ggplot2::aes(ymax = ci ,ymin = -ci),color = NA,fill = "blue",alpha = 0.2)
  p5 = p5 + ggplot2::ggtitle('ACF and PACF plots for residuals series') + ggplot2::theme_bw()
  return(list(p1,p2,p3,p4,p5))
}
