#==================================================================================================#
# Date: 14/04/2020
# Description: Display some graphics for residuals analysis
# Function:
#==================================================================================================#
diagnostic_mtar = function(regime_model,lagmax = NULL, CusumSQ = NULL){
  if (!requireNamespace('ggplot2', quietly = TRUE)) {
    stop('ggplot2 is needed for this function to work')
  }else {
    if (!inherits(regime_model, 'regime_model')) {
      stop('diagnostic.mtar requires a regime_model object')
    }}
  regime_model$residuals = as.data.frame(regime_model$residuals[-nrow(regime_model$residuals),])
  e_k = tsregime(as.matrix(regime_model$residuals))
  p1 = autoplot.tsregime(e_k) + ggplot2::geom_hline(yintercept = 0,color = "red") +
    ggplot2::ggtitle('Residual serie plot')
  e_data = as.data.frame(e_k$Yt)
  time = seq(1,nrow(e_data))
  dat = data.frame(label = 'Series.1',time = time,value = e_data[,1],
                   cusum = cumsum(e_data[,1])/stats::sd(e_data[,1]),
                   cumsq = c(cumsum(e_data[,1]^2)/sum(e_data[,1]^2)))
  if (ncol(e_data) > 1) {
    for (i in 2:ncol(e_data)) {
      dat = rbind.data.frame(dat,data.frame(label = paste0('Series.',i),time = time,value = e_data[,i],
                                 cusum = cumsum(e_data[,i])/stats::sd(e_data[,i]),
                                 cumsq = c(cumsum(e_data[,i]^2)/sum(e_data[,i]^2))))
    }
  }
  p2 = ggplot2::ggplot(ggplot2::aes_(x = ~value, color = ~label),data = dat) +
    ggplot2::geom_density() + ggplot2::theme_bw()
  p2 = p2 + ggplot2::stat_function(fun = stats::dnorm,color = "black")
  p2 = p2 + ggplot2::ggtitle("Residual density plot")

  Af = 0.948 ###Cuantil del 95% para cusum
  LS = Af*sqrt(e_k$N) + 2*Af*c(1:e_k$N)/sqrt(e_k$N)
  LI = -LS
  p3 = ggplot2::ggplot(ggplot2::aes_(x = ~time, y = ~cusum,color = ~label),data = dat)
  p3 = p3 + ggplot2::geom_ribbon(ggplot2::aes(ymin = rep(LS,e_k$k), ymax = rep(LI,e_k$k)),
                                 fill = "gray",color = NA,alpha = 0.5)
  p3 = p3 + ggplot2::geom_line() + ggplot2::theme_bw() + ggplot2::ggtitle('CUSUM statistic for residuals')
# Tabla CusumSQ
  tablasq = data.frame(n = c(61,63,65,67,69,71,73,75,77,79,81,83,85,87,89,91,93,95,97,99,100,110,120,
                             130,140,150,160,170,180,190,200,300,400,500,600,700,800,900,1000,2000,
                             3000,4000,5000,6000,7000,8000,9000,10000),
                       alpha0.1 = c(0.12522,0.12342,0.12170,0.12006,0.11848,0.11696,0.11550,
                                    0.11409,0.11273,0.11143,0.11017,0.10895,0.10777,0.10663,
                                    0.10553,0.10446,0.10342,0.10241,0.10144,0.10049,0.10002,0.09571,
                                    0.09193,0.08856,0.08555,0.08283,0.08035,0.07809,0.07601,
                                    0.07409,0.07231,0.05960,0.05190,0.04659,0.04265,0.03957,0.03707,
                                    0.03500,0.03324,0.02365,0.01936,0.01680,0.01504,0.01374,
                                    0.01273,0.01191,0.01124,0.01066),
                       alpha0.05 = c(0.14422,0.14213,0.14013,0.13821,0.13637,0.13461,0.13291,
                                     0.13128,0.12970,0.12819,0.12672,0.12531,0.12394,0.12262,
                                     0.12134,0.12010,0.11889,0.11773,0.11660,0.11550,0.11496,
                                     0.10997,0.10558,0.10169,0.09821,0.09506,0.09220,0.08959,0.08719,
                                     0.08498,0.08293,0.06828,0.05943,0.05333,0.04880,0.04526,
                                     0.04240,0.04002,0.03801,0.02702,0.02212,0.01918,0.01717,
                                     0.01569,0.01453,0.01360,0.01283,0.01217),
                       alpha0.025 = c(0.16109,0.15874,0.15649,0.15433,0.15227,0.15028,0.14838,
                                      0.14654,0.14478,0.14307,0.14143,0.13984,0.13831,0.13682,
                                      0.13538,0.13399,0.13264,0.13134,0.13007,0.12883,0.12823,0.12263,
                                      0.11772,0.11336,0.10946,0.10594,0.10274,0.09982,0.09714,0.09466,
                                      0.09237,0.07601,0.06612,0.05932,0.05427,0.05033,0.04714,
                                      0.04449,0.04225,0.03002,0.02457,0.02130,0.01907,0.01742,0.01614,
                                      0.01510,0.01424,0.01351),
                       alpha0.01 = c(0.18107,0.17841,0.17587,0.17344,0.17110,0.16886,0.16671,
                                     0.16463,0.16264,0.16071,0.15886,0.15706,0.15533,0.15366,
                                     0.15204,0.15046,0.14894,0.14747,0.14603,0.14464,0.14396,0.13765,
                                     0.13211,0.12720,0.12280,0.11884,0.11524,0.11195,0.10893,
                                     0.10614,0.10356,0.08517,0.07406,0.06642,0.06076,0.05634,
                                     0.05276,0.04980,0.04728,0.03358,0.02747,0.02382,0.02132,
                                     0.01948,0.01804,0.01688,0.01592,0.01511),
                       alpha0.005 = c(0.19486,0.19199,0.18925,0.18662,0.18410,0.18168,0.17936,
                                      0.17712,0.17497,0.17289,0.17089,0.16896,0.16709,0.16528,
                                      0.16353,0.16184,0.16020,0.15861,0.15706,0.15556,0.15483,
                                      0.14803,0.14206,0.13676,0.13202,0.12775,0.12387,0.12033,
                                      0.11708,0.11408,0.11130,0.09150,0.07955,0.07134,0.06525,0.06049,
                                      0.05665,0.05346,0.05076,0.03605,0.02949,0.02556,0.02288,
                                      0.02090,0.01936,0.01811,0.01708,0.01621))
  if (!is.null(CusumSQ)) {
    co = CusumSQ ####Valor del cuantil aproximado para cusumsq para T/2-1=71 alpha=0.05
    LQS = co + (1:e_k$N)/e_k$N
    LQI = -co + (1:e_k$N)/e_k$N
    p4 = ggplot2::ggplot(ggplot2::aes_(x = ~time, y = ~cumsq,color = ~label),data = dat)
    p4 = p4 + ggplot2::geom_ribbon(ggplot2::aes(ymin = rep(LQS,e_k$k), ymax = rep(LQI,e_k$k)),
                                   fill = "gray",color = NA,alpha = 0.5)
    p4 = p4 + ggplot2::geom_line() + ggplot2::theme_bw() + ggplot2::ggtitle('CUSUMSQ statistic for residuals')

  }
  acf_i = stats::acf(regime_model$residuals[,1],lag.max = lagmax,plot = FALSE,type = 'correlation')
  acf_Yt = data.frame(Lag = acf_i$lag, value = acf_i$acf,names = 'Serie.1',type = 'ACF')
  pacf_i = stats::acf(regime_model$residuals[,1],lag.max = lagmax,plot = FALSE,type = 'partial')
  pacf_Yt = data.frame(Lag = pacf_i$lag, value = pacf_i$acf,names = 'Serie.1',type = 'PACF')
  if (ncol(regime_model$residuals) > 1) {
    for (i in 2:ncol(regime_model$residuals)) {
      acf_i = stats::acf(regime_model$residuals[,i],lag.max = lagmax,plot = FALSE,type = 'correlation')
      acf_Yt = rbind(acf_Yt,data.frame(Lag = acf_i$lag + 0.1*i, value = acf_i$acf,names = paste0('Series.',i),type = 'ACF'))
      pacf_i = stats::acf(regime_model$residuals[,i],lag.max = lagmax,plot = FALSE,type = 'partial')
      pacf_Yt = rbind(pacf_Yt,data.frame(Lag = pacf_i$lag + 0.1*i, value = pacf_i$acf,names = paste0('Series.',i),type = 'PACF'))
    }
  }
  dat_cor = rbind.data.frame(acf_Yt,pacf_Yt)
  p5 = ggplot2::ggplot(ggplot2::aes_(x = ~Lag, y = ~value),data = dat_cor[floor(dat_cor$Lag) != 0,])
  p5 = p5 + ggplot2::geom_hline(yintercept = 0) + ggplot2::facet_grid(type~names)
  p5 = p5 + ggplot2::geom_segment(ggplot2::aes(xend = dat_cor[floor(dat_cor$Lag) != 0,]$Lag,yend = 0)) + ggplot2::geom_point(color = "blue",size = 0.4)
  ci = stats::qnorm((1 + 0.95)/2)/sqrt(nrow(regime_model$residuals))
  p5 = p5 + ggplot2::geom_ribbon(ggplot2::aes(ymax = ci ,ymin = -ci),color = NA,fill = "blue",alpha = 0.2)
  p5 = p5 + ggplot2::ggtitle('ACF and PACF plots for residuals series') + ggplot2::theme_bw()
  return(list(p1,p2,p3,p4,p5))
}
