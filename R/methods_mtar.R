autoplot = function(object, ...) UseMethod("autoplot")
autoplot.regime_model = function(object, type = 1, ...) {
  if (!requireNamespace('ggplot2', quietly = TRUE)) {
    stop('ggplot2 is needed for this function to work')
  }else {
    if (!inherits(object, 'regime_model')) {
      stop('autoplot.regime_model requires a regime_model object')
    }}
  if (!{type %in% c(1:5)}) {stop('type should take values in c (1,2,3,4)')}
  if (is.null(object$Chain)) {stop('There are no chains to graph')}
  if (type == 1) {
    if (is.null(object$Chain$r)) {stop('r unknown')}
    Chain_r = t(object$Chain$r)
    time = seq(1,nrow(Chain_r))
    dat2 = data.frame(name = 'r.1',time = time,value = Chain_r[,1])
    if (ncol(Chain_r) > 1) {
      for (i in 2:ncol(Chain_r)) {
        dat2 = rbind(dat2,data.frame(name = paste0('r.',i),
                                     time = time,value = Chain_r[,i]))
      }
    }
    p = ggplot2::ggplot(ggplot2::aes_(x = ~time,y = ~value),data = dat2)
    p = p + ggplot2::geom_line() + ggplot2::facet_grid(name~.,scales = 'free') + ggplot2::ggtitle('Threshold variable chains')
    p = p + ggplot2::theme_bw() + ggplot2::scale_y_continuous(labels = function(x) format(x, scientific = TRUE))
    return(p)
    }
  if (type == 2) {
    # Sigma Chains
    Chain_Sig = object$Chain$Sigma
    dat3l = vector(mode = 'list',length = length(Chain_Sig))
    p3 = vector(mode = 'list',length = length(Chain_Sig))
    names(p3) = names(Chain_Sig)
    names(dat3l) = names(Chain_Sig)
    for (j in names(Chain_Sig)) {
      if (!is.matrix(Chain_Sig[[j]])) {
        Chain_Sig[[j]] = t(as.matrix(Chain_Sig[[j]]))
      }
      time = seq(1,ncol(Chain_Sig[[j]]))
      dat3 = data.frame(comp = '11',time = time,value = Chain_Sig[[j]][1,])
      k = dim(object$regime[[j]]$sigma)[1]
      names_sig = paste0(1:k,1)
      for (i3 in 2:k) {names_sig = c(names_sig,paste0(1:k,i3))}
      if (nrow(Chain_Sig[[j]]) > 1) {
        ii = 1
        for (i in names_sig[-1]) {
          dat3 = rbind(dat3,data.frame(comp = i,time = time,value = Chain_Sig[[j]][ii,]))
          ii = ii + 1
        }
      }
      p3[[j]] =  ggplot2::ggplot(ggplot2::aes_(x = ~time, y = ~value),data = dat3) +
        ggplot2::geom_line() + ggplot2::facet_grid(dat3$comp~.,scales = 'free') +
        ggplot2::theme_bw() +
        ggplot2::labs(title = paste('Sigma chains',j)) +
        ggplot2::scale_y_continuous(labels = function(x) format(x, scientific = TRUE))
      }
    return(p3)
    }
  if (type == 3) {
    # Theta Chains
    Chain_Theta = object$Chain$Theta
    dat3l = vector(mode = 'list',length = length(Chain_Theta))
    p4 = vector(mode = 'list',length = length(Chain_Theta))
    names(p4) = names(Chain_Theta)
    if (!is.matrix(Chain_Theta$R1)) {
      Chain_Theta$R1 = t(as.matrix(Chain_Theta$R1))
    }
    time = seq(1,ncol(Chain_Theta$R1))
    for (j in names(Chain_Theta)) {
      dat3 = NULL
      for (i in 1:nrow(Chain_Theta[[j]])) {
        dat3 = rbind(dat3,data.frame(comp = rownames(object$estimates$Theta[[j]])[i],
                                     time = time,value = Chain_Theta[[j]][i,]))
      }
      dat3$name = sub('\\..*','',dat3$comp)
      dat3$num = sub('.*\\.','',dat3$comp)
      dat3$name2 = gsub('[0-9]+','', dat3$name)
      dat3$name2[dat3$name == 'phi0'] = 'c'
      p4[[j]] = vector(mode = 'list',length = length(unique(dat3$name2)))
      names(p4[[j]]) = unique(dat3$name2)
      for (nn in unique(dat3$name2)) {
        p4[[j]][[nn]] = ggplot2::ggplot(ggplot2::aes_(x = ~time,y = ~value),data = dat3[dat3$name2 == nn,]) +
          ggplot2::theme_bw() +
          ggplot2::geom_line() + ggplot2::facet_grid(comp~.,scales = 'free') +
          ggplot2::labs(title = paste(nn,'chains',j)) +
          ggplot2::scale_y_continuous(labels = function(x) format(x, scientific = TRUE))
       }
      }
    return(p4)
  }
  if (type == 4) {
    # Gamma Chains
    if (is.null(object$Chain$Gamma)) {
      stop('Object $Chain$Gamma does not exist (known orders)')
    }
    Chain_Gamma = object$Chain$Gamma
    dat3l = vector(mode = 'list',length = length(Chain_Gamma))
    p5 = vector(mode = 'list',length = length(Chain_Gamma))
    names(p5) = names(Chain_Gamma)
    time = seq(1,ncol(Chain_Gamma$R1))
    for (j in names(Chain_Gamma)) {
      dat3 = NULL
      for (i in 1:nrow(Chain_Gamma[[j]])) {
        dat3 = rbind(dat3,data.frame(comp = rownames(object$estimates$Gamma[[j]])[i],
                                     time = time,value = Chain_Gamma[[j]][i,]))
      }
      dat3$name = sub('\\..*','',dat3$comp)
      dat3$num = sub('.*\\.','',dat3$comp)
      dat3$name2 = gsub('[0-9]+','', dat3$name)
      dat3$name2[dat3$name == 'phi0'] = 'c'
      p5[[j]] = vector(mode = 'list',length = length(unique(dat3$name2)))
      names(p5[[j]]) = unique(dat3$name2)
      for (nn in unique(dat3$name2)) {
        p5[[j]][[nn]] = ggplot2::ggplot(ggplot2::aes_(x = ~time,y = ~value),data = dat3[dat3$name2 == nn,]) +
          ggplot2::theme_bw() +
          ggplot2::geom_area() + ggplot2::facet_grid(comp~.,scales = 'free') +
          ggplot2::labs(title = paste('Gamma chains',nn,j)) +
          ggplot2::scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
          ggplot2::theme(axis.ticks.y = ggplot2::element_blank(),axis.text.y = ggplot2::element_blank())
      }
    }
    return(p5)
  }
  if (type == 5) {
    Chain_Yt = as.data.frame(object$data$Yt)
    Chain_fit = as.data.frame(object$fitted.values)
    Chain_Yt = as.data.frame(Chain_Yt[-nrow(Chain_Yt),])
    Chain_fit = as.data.frame(Chain_fit[-nrow(Chain_fit),])
    time = seq(1,nrow(Chain_fit))
    dat1 = data.frame(type = 'obs',name = 'Series.1',time = time, value = Chain_Yt[,1])
    dat1 = rbind.data.frame(dat1,data.frame(type = 'fit',name = 'Series.1',time = time, value = Chain_fit[,1]))
    if (ncol(Chain_Yt) > 1) {
      for (i in 2:ncol(Chain_Yt)) {
        dati = data.frame(type = 'obs',name = paste0('Series.',i),time = time,value = Chain_Yt[,i])
        dati = rbind.data.frame(dati,data.frame(type = 'fit',name = paste0('Series.',i),time = time, value = Chain_fit[,i]))
        dat1 = rbind.data.frame(dat1,dati)
      }
    }
    p = ggplot2::ggplot(ggplot2::aes_(x = ~time,y = ~value, color = ~type),data = dat1)
    p = p + ggplot2::geom_line() + ggplot2::facet_grid(name~.) + ggplot2::theme_bw()
    p = p + ggplot2::labs(title = 'Output process')
    p = p + ggplot2::scale_color_manual(values = c("black","blue")) +
      ggplot2::scale_y_continuous(labels = function(x) format(x, scientific = TRUE))
    return(p)
  }
}
autoplot.regime_missing = function(object, type = 1, ...) {
    if (!requireNamespace('ggplot2', quietly = TRUE)) {
      stop('ggplot2 is needed for this function to work')
    }else {
      if (!inherits(object, 'regime_missing')) {
        stop('autoplot.regime_missing requires a regime_missing object')
      }}
    if (is.null(object$Chains$Y)) {stop('There are no chains to graph')}
    if (!{type %in% c(1:4)}) {stop('type should take values in c (1,2,3,4)')}
    if (type == 1) {
      if (is.null(object$estimates$Yt)) {stop('Yt has no missing data')}
      Chain_Yt = t(object$Chains$Yt)
      time = seq(1,nrow(Chain_Yt))
      names_yt = rownames(object$estimates$Yt)
      dat2 = data.frame(name = names_yt[1],time = time,value = Chain_Yt[,1])
      if (ncol(Chain_Yt) > 1) {
        for (i in 2:ncol(Chain_Yt)) {
          dat2 = rbind(dat2,data.frame(name = names_yt[i],time = time,value = Chain_Yt[,i]))
        }
      }
      p = ggplot2::ggplot(ggplot2::aes_(x = ~time,y = ~value),data = dat2)
      p = p + ggplot2::geom_line() + ggplot2::facet_grid(name~.,scales = 'free') + ggplot2::theme_bw()
      p = p +  ggplot2::labs(title = 'Missing data (Yt) chains')
      return(p)
    }
    if (type == 2) {
      if (is.null(object$estimates$Zt)) {stop('Zt has no missing data')}
      Chain_Zt = t(object$Chains$Zt)
      time = seq(1,nrow(Chain_Zt))
      names_Zt = rownames(object$estimates$Zt)
      dat2 = data.frame(name = names_Zt[1],time = time,value = Chain_Zt[,1])
      if (ncol(Chain_Zt) > 1) {
        for (i in 2:ncol(Chain_Zt)) {
          dat2 = rbind(dat2,data.frame(name = names_Zt[i],time = time,value = Chain_Zt[,i]))
        }
      }
      p = ggplot2::ggplot(ggplot2::aes_(x = ~time,y = ~value),data = dat2)
      p = p + ggplot2::geom_line() + ggplot2::facet_grid(name~.,scales = 'free') + ggplot2::theme_bw()
      p = p +  ggplot2::labs(title = 'Missing data (Zt) chains')
      return(p)
    }
    if (type == 3) {
      if (is.null(object$estimates$Xt)) {stop('Xt has no missing data')}
      Chain_Xt = t(object$Chains$Xt)
      time = seq(1,nrow(Chain_Xt))
      names_Xt = rownames(object$estimates$Xt)
      dat2 = data.frame(name = names_Xt[1],time = time,value = Chain_Xt[,1])
      if (ncol(Chain_Xt) > 1) {
        for (i in 2:ncol(Chain_Xt)) {
          dat2 = rbind(dat2,data.frame(name = names_Xt[i],time = time,value = Chain_Xt[,i]))
        }
      }
      p = ggplot2::ggplot(ggplot2::aes_(x = ~time,y = ~value),data = dat2)
      p = p + ggplot2::geom_line() + ggplot2::facet_grid(name~.,scales = 'free') + ggplot2::theme_bw()
      p = p +  ggplot2::labs(title = 'Missing data (Xt) chains')
      return(p)
    }
    if (type == 4) {
      list_plots_miss = list()
      time = seq(1,object$tsregime$N)
      if ('Yt' %in% names(object$estimates)) {
        dats_Yt = t(object$tsregime$Yt)
        dats_Yt[object$pos_na[[1]]] = NA
        dats_Yt = t(dats_Yt)
        dats_Yt_NA_mean = t(dats_Yt*NA)
        dats_Yt_NA_mean[object$pos_na[[1]]] = object$estimates$Yt[,2]
        dats_Yt_NA_up = t(dats_Yt*NA)
        dats_Yt_NA_up[object$pos_na[[1]]] = object$estimates$Yt[,3]
        dats_Yt_NA_low = t(dats_Yt*NA)
        dats_Yt_NA_low[object$pos_na[[1]]] = object$estimates$Yt[,1]
        dat = data.frame(name = 'Series.1',time = time,
                         value = dats_Yt[,1],mean_miss = dats_Yt_NA_mean[1,],up_miss = dats_Yt_NA_up[1,],low_miss = dats_Yt_NA_low[1,])
        if (ncol(dats_Yt) > 1) {
          for (i in 2:ncol(dats_Yt)) {
            dat = rbind(dat,data.frame(name = paste0('Series.',i),time = time,value = dats_Yt[,i],
                                       mean_miss = dats_Yt_NA_mean[i,],up_miss = dats_Yt_NA_up[i,],low_miss = dats_Yt_NA_low[i,]))
          }
        }
        p = ggplot2::ggplot(ggplot2::aes_(x = ~time,y = ~value),data = dat)
        p = p + ggplot2::geom_line() + ggplot2::theme_bw()
        p = p + ggplot2::geom_errorbar(ggplot2::aes_(ymin = ~low_miss, ymax =  ~up_miss),width = 0.2,color = 'blue')
        p = p + ggplot2::geom_point(ggplot2::aes_(x = ~time,y = ~mean_miss),color = 'blue')
        p = p + ggplot2::labs(title = 'Output process')
        p = p + ggplot2::facet_grid(name~.,scales = 'free_y')
        list_plots_miss[[1]] = p
      }
      if (!is.null(object$tsregime$Zt) & ('Zt' %in% names(object$estimates))) {
        dats_Zt = object$tsregime$Zt
        dats_Zt[object$pos_na[[2]][1,]] = NA
        dats_Zt_NA_mean = dats_Zt*NA
        dats_Zt_NA_mean[object$pos_na[[2]][1,]] = object$estimates$Zt[,2]
        dats_Zt_NA_up =  dats_Zt*NA
        dats_Zt_NA_up[object$pos_na[[2]][1,]] = object$estimates$Zt[,3]
        dats_Zt_NA_low =  dats_Zt*NA
        dats_Zt_NA_low[object$pos_na[[2]][1,]] = object$estimates$Zt[,1]
        dat = data.frame(time = time,
                         value = dats_Zt,mean_miss = dats_Zt_NA_mean,up_miss = dats_Zt_NA_up,low_miss = dats_Zt_NA_low)
        p2 = ggplot2::ggplot(ggplot2::aes_(x = ~time,y = ~value),data = dat)
        p2 = p2 + ggplot2::geom_line() + ggplot2::theme_bw()
        p2 = p2 + ggplot2::geom_errorbar(ggplot2::aes_(ymin = ~low_miss, ymax =  ~up_miss),width = 0.2,color = 'blue')
        p2 = p2 + ggplot2::geom_point(ggplot2::aes_(x = ~time,y = ~mean_miss),color = 'blue')
        if (!is.null(object$r)) {
          Nrg_plot = paste0(paste0(paste0('Reg_',1:object$l),'='),object$Summary_r$Prop_reg,'%')
          p2 = p2 + ggplot2::labs(title = 'Threshold process',subtitle = paste0('(',paste(Nrg_plot,collapse = ','),')'))
          for (i in c(object$r)) {
            p2 = p2 + ggplot2::geom_hline(yintercept = i,linetype = 'dashed',color = 'blue')
          }
        }else{
          p2 = p2 + ggplot2::labs(title = 'Threshold process')
        }
        list_plots_miss[[2]] = p2
      }
      if (!is.null(object$tsregime$Xt) & ('Xt' %in% names(object$estimates))) {
        dats_Xt = t(object$tsregime$Xt)
        dats_Xt[object$pos_na[[2]][-1,]] = NA
        dats_Xt = t(dats_Xt)
        dats_Xt_NA_mean = t(dats_Xt*NA)
        dats_Xt_NA_mean[object$pos_na[[2]][-1,]] = object$estimates$Xt[,2]
        dats_Xt_NA_up = t(dats_Xt*NA)
        dats_Xt_NA_up[object$pos_na[[2]][-1,]] = object$estimates$Xt[,3]
        dats_Xt_NA_low = t(dats_Xt*NA)
        dats_Xt_NA_low[object$pos_na[[2]][-1,]] = object$estimates$Xt[,1]
        dat = data.frame(name = 'Series.1',time = time,
                         value = dats_Xt[,1],mean_miss = dats_Xt_NA_mean[1,],up_miss = dats_Xt_NA_up[1,],low_miss = dats_Xt_NA_low[1,])
        if (ncol(dats_Xt) > 1) {
          for (i in 2:ncol(dats_Xt)) {
            dat = rbind(dat,data.frame(name = paste0('Series.',i),time = time,value = dats_Xt[,i],
                                       mean_miss = dats_Xt_NA_mean[i,],up_miss = dats_Xt_NA_up[i,],low_miss = dats_Xt_NA_low[i,]))
          }
        }
        p3 = ggplot2::ggplot(ggplot2::aes_(x = ~time,y = ~value),data = dat)
        p3 = p3 + ggplot2::geom_line() + ggplot2::theme_bw()
        p3 = p3 + ggplot2::geom_errorbar(ggplot2::aes_(ymin = ~low_miss, ymax =  ~up_miss),width = 0.2,color = 'blue')
        p3 = p3 + ggplot2::geom_point(ggplot2::aes_(x = ~time,y = ~mean_miss),color = 'blue')
        p3 = p3 + ggplot2::labs(title = 'Covariates process')
        p3 = p3 + ggplot2::facet_grid(name~.,scales = 'free_y')
        list_plots_miss[[3]] = p3
      }
      return(list_plots_miss)
    }
  }
autoplot.tsregime = function(object, type = 1, ...) {
  if (!requireNamespace('ggplot2', quietly = TRUE)) {
    stop('ggplot2 is needed for this function to work')
  }else {
    if (!inherits(object, 'tsregime')) {
      stop('autoplot.tsregime requires a tsregime object')
    }}
  if (!{type %in% c(1:3)}) {stop('type should take values in c (1,2,3)')}
  dats_Yt = as.data.frame(object$Yt)
  time = seq(1,nrow(dats_Yt))
  dat = data.frame(name = 'Series.1',time = time,value = dats_Yt[,1])
  if (ncol(dats_Yt) > 1) {
    for (i in 2:ncol(object$Yt)) {
      dat = rbind(dat,data.frame(name = paste0('Series.',i),time = time,value = dats_Yt[,i]))
    }
  }
  dat_NA = c()
  N = length(dats_Yt[,1])
  for (i in 1:object$k) {
    xl = c(1:N)[is.na(dats_Yt[,i])]
    dat_NA = rbind(dat_NA,data.frame(name = rep(paste0('Series.',i),length(xl)),xl = xl))
  }
  p = ggplot2::ggplot(ggplot2::aes_(x = ~time,y = ~value),data = dat)
  p = p + ggplot2::geom_line() + ggplot2::theme_bw()
  p = p + ggplot2::labs(title = 'Output process')
  p = p + ggplot2::geom_vline(ggplot2::aes(xintercept = xl),color = "red",linetype = 'dashed',data = dat_NA)
  p = p + ggplot2::facet_grid(name~.,scales = 'free_y')
  if (!is.null(object$Zt)) {
    dats_Zt = data.frame(time = time,value = object$Zt)
    p2 = ggplot2::ggplot(ggplot2::aes_(x = ~time,y = ~value),data = dats_Zt)
    p2 = p2 + ggplot2::geom_line() + ggplot2::theme_bw()
    p2 = p2 + ggplot2::geom_vline(xintercept = dats_Zt$time[is.na(dats_Zt$value)],color = "red",linetype = 'dashed')
    if (!is.null(object$r)) {
      Nrg_plot = paste0(paste0(paste0('Reg_',1:object$l),'='),object$Summary_r$Prop_reg,'%')
      p2 = p2 + ggplot2::labs(title = 'Threshold process',subtitle = paste0('(',paste(Nrg_plot,collapse = ','),')'))
      for (i in c(object$r)) {
        p2 = p2 + ggplot2::geom_hline(yintercept = i,linetype = 'dashed',color = 'blue')
      }
    }else{
      p2 = p2 + ggplot2::labs(title = 'Threshold process')
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
    dat_NA = c()
    for (i in 1:object$nu) {
      xl = c(1:N)[is.na(dats_Xt[,i])]
      dat_NA = rbind(dat_NA,data.frame(name = rep(paste0('Series.',i),length(xl)),xl = xl))
    }
    p3 = ggplot2::ggplot(ggplot2::aes_(x = ~time,y = ~value),data = dat2)
    p3 = p3 + ggplot2::geom_line() + ggplot2::theme_bw()
    p3 = p3 + ggplot2::labs(title = 'Covariates process')
    p3 = p3 + ggplot2::geom_vline(ggplot2::aes(xintercept = xl),color = "red",linetype = 'dashed',data = dat_NA)
    p3 = p3 + ggplot2::facet_grid(name~.,scales = 'free_y')
  }
  if (type == 1) {
    return(p)
  }
  if (type == 2) {
    if (is.null(object$Zt)) {
      stop('Threshold process does not exist')}
    return(p2)
  }
  if (type == 3) {
    if (is.null(object$Xt)) {
      stop('Covariates process does not exist')}
    return(p3)
  }
}
autoplot.regime_number = function(object, type = 1, ...){
  if (!requireNamespace('ggplot2', quietly = TRUE)) {
    stop('ggplot2 is needed for this function to work')
  }else {
    if (!inherits(object,'regime_number')) {
      stop('autoplot.regime_number requires a regime_number object')
    }}
  if (!{type %in% c(1:2)}) {stop('type should take values in c (1,2)')}
  if (type == 1){
    plot1_list = list()
    for (m_i in names(object$list_m)) {
      tex_aux = paste(round(rev(rev(object$list_m[[m_i]]$par$r)[-1]),3),collapse=" | ")
      plot1_list[[m_i]] = autoplot.regime_model(object$list_m[[m_i]]$par,1) +
        ggplot2::ggtitle(paste0('Threshold variable chains (r = ',tex_aux,')')) +
        ggplot2::geom_smooth()
    }
    return(plot1_list)
  }
  if (type == 2){
    plot2_list = list()
    for (m_i in names(object$list_m)) {
      tex_aux = paste(round(rev(rev(object$list_m[[m_i]]$par$r)[-1]),3),collapse=" | ")
      data_i = object$list_m[[m_i]]
      r_i = rev(rev(data_i$par$r)[-1])
      tsregime_i = tsregime(Yt = data_i$par$data$Yt,Zt = data_i$par$data$Zt,Xt = data_i$par$data$Xt,r = r_i)
      plot2_list[[m_i]] = autoplot.tsregime(tsregime_i,2) +
        ggplot2::ggtitle(paste0('Threshold variable chains (r = ',tex_aux,')'))
      cat(paste(m_i,'===========================|','\n'))
      print(tsregime_i$Summary_r)
    }
    return(plot2_list)
  }
}
autoplot.regime_forecast = function(object, type = 1, ...){
  if (!requireNamespace('ggplot2', quietly = TRUE)) {
    stop('ggplot2 is needed for this function to work')
  }else {
    if (!inherits(object, 'regime_forecast')) {
      stop('autoplot.regime_forecast requires a regime_forecast object')
    }}
  if (!{type %in% c(1:2)}) {stop('type should take values in c (1,2)')}
  if (is.null(object$forecast$chain)) {stop('There are no chains to graph')}
  if (type == 1) {
    plot1_list = list()
    for (name_chain in names(object$forecast$estim)) {
      dat3 = NULL
      if (name_chain %in% c('Xt','Zt')) {
        if ('Xt' %in% c('Xt','Zt')){
          if (name_chain == 'Xt') {
            id_names = !{rownames(object$forecast$estim_Ut) %in% rownames(object$forecast$estim$Zt)}
          }else{
            id_names = rownames(object$forecast$estim_Ut) %in% rownames(object$forecast$estim$Zt)
          }
          Chain_t = t(object$forecast$chain$Ut)[id_names,]
          time = seq(1,ncol(Chain_t))
          for (i in 1:nrow(Chain_t)) {
            dat3 = rbind(dat3,data.frame(comp = rownames(object$forecast$estim[[name_chain]])[i],
                                         time = time,value = Chain_t[i,]))
          }
        }else{
          Chain_t = t(object$forecast$chain$Ut)
          time = seq(1,ncol(Chain_t))
          for (i in 1:nrow(Chain_t)) {
            dat3 = rbind(dat3,data.frame(comp = rownames(object$forecast$estim[[name_chain]])[i],
                                         time = time,value = Chain_t[i,]))
          }
        }
      }
      if (name_chain == 'Yt') {
        Chain_t = t(object$forecast$chain[[name_chain]])
        time = seq(1,ncol(Chain_t))
        for (i in 1:nrow(Chain_t)) {
          dat3 = rbind(dat3,data.frame(comp = rownames(object$forecast$estim[[name_chain]])[i],
                                       time = time,value = Chain_t[i,]))
        }
      }
      plot1_list[[name_chain]] = ggplot2::ggplot(ggplot2::aes_(x = ~time,y = ~value),data = dat3) +
        ggplot2::theme_bw() +
        ggplot2::geom_line() + ggplot2::facet_grid(comp~.,scales = 'free') +
        ggplot2::labs(title = paste0('Chains forecast ',name_chain))
    }
    return(plot1_list)
  }
  if (type == 2) {
    plot2_list = list()
    dats_Yt = as.data.frame(object$tsregime$Yt)
    aux_vec = dats_Yt*NA
    forecast_mean = aux_vec
    forecast_mean[as.numeric(colnames(object$forecast$Yth)),] =  t(object$forecast$Yth)
    forecast_low = aux_vec
    forecast_low[as.numeric(colnames(object$forecast$Yth)),] = t(ks::invvec(object$forecast$estim$Yt[,1],ncol = length(colnames(object$forecast$Yth)),nrow = object$tsregime$k))
    forecast_Up = aux_vec
    forecast_Up[as.numeric(colnames(object$forecast$Yth)),] = t(ks::invvec(object$forecast$estim$Yt[,3],ncol = length(colnames(object$forecast$Yth)),nrow = object$tsregime$k))
    time = seq(1,nrow(dats_Yt))
    dat = NULL
    for (i in 1:ncol(dats_Yt)) {
      dat = rbind(dat,data.frame(name = paste0('Series.',i),time = time,
                                 value = dats_Yt[,i],low = forecast_low[,i],up = forecast_Up[,i],mean = forecast_mean[,i]))
    }
    p1 = ggplot2::ggplot(ggplot2::aes_(x = ~time,y = ~value),data = dat)
    p1 = p1 + ggplot2::geom_line() + ggplot2::theme_bw()
    p1 = p1 + ggplot2::geom_ribbon(ggplot2::aes_(ymin = ~low,ymax = ~up),fill = 'blue',alpha = 0.4)
    p1 = p1 + ggplot2::geom_line(ggplot2::aes_(x = ~time,y = ~mean),color = 'blue')
    p1 = p1 + ggplot2::labs(title = 'Output process and forecast')
    p1 = p1 + ggplot2::facet_grid(name~.,scales = 'free_y')
    plot2_list[['Yt']] = p1
    if (!is.null(object$tsregime$Zt)){
      dats_Zt = as.data.frame(object$tsregime$Zt)
      aux_vec = dats_Zt*NA
      forecast_mean = aux_vec
      forecast_mean[as.numeric(colnames(object$forecast$Zth)),] =  c(object$forecast$Zth)
      forecast_low = aux_vec
      forecast_low[as.numeric(colnames(object$forecast$Zth)),] = c(object$forecast$estim$Zt[,1])
      forecast_Up = aux_vec
      forecast_Up[as.numeric(colnames(object$forecast$Zth)),] = c(object$forecast$estim$Zt[,3])
      time = seq(1,nrow(dats_Zt))
      dat = data.frame(time = c(time),value = dats_Zt[,1],low = forecast_low[,1],up = forecast_Up[,1],mean = forecast_mean[,1])
      p2 = ggplot2::ggplot(ggplot2::aes_(x = ~time,y = ~value),data = dat)
      p2 = p2 + ggplot2::geom_line() + ggplot2::theme_bw()
      p2 = p2 + ggplot2::geom_ribbon(ggplot2::aes_(ymin = ~low,ymax = ~up),fill = 'blue',alpha = 0.4)
      p2 = p2 + ggplot2::geom_line(ggplot2::aes_(x = ~time,y = ~mean),color = 'blue')
      Nrg_plot = paste0(paste0(paste0('Reg_',1:object$tsregime$l),'='),object$tsregime$Summary_r$Prop_reg,'%')
      p2 = p2 + ggplot2::labs(title = 'Threshold process and forecast',
                              subtitle = paste0('(',paste(Nrg_plot,collapse = ','),')'))
      for (i in c(object$tsregime$r)) {
        p2 = p2 + ggplot2::geom_hline(yintercept = i,linetype = 'dashed',color = 'blue')
      }
      plot2_list[['Zt']] = p2
    }
    if (!is.null(object$tsregime$Xt)) {
      if (object$tsregime$nu == 1) {
        names_xt = as.numeric(names(object$forecast$Xth))
      }else{
        names_xt = as.numeric(colnames(object$forecast$Xth))
      }
      dats_Xt = as.data.frame(object$tsregime$Xt)
      aux_vec = dats_Xt*NA
      forecast_mean = aux_vec
      forecast_mean[names_xt,] = t(ks::invvec(object$forecast$estim$Xt[,2],ncol = length(names_xt),nrow = object$tsregime$nu))
      forecast_low = aux_vec
      forecast_low[names_xt,] = t(ks::invvec(object$forecast$estim$Xt[,1],ncol = length(names_xt),nrow = object$tsregime$nu))
      forecast_Up = aux_vec
      forecast_Up[names_xt,] = t(ks::invvec(object$forecast$estim$Xt[,3],ncol = length(names_xt),nrow = object$tsregime$nu))
      time = seq(1,nrow(dats_Xt))
      dat = NULL
      for (i in 1:ncol(dats_Xt)) {
        dat = rbind(dat,data.frame(name = paste0('Series.',i),time = time,
                                   value = dats_Xt[,i],low = forecast_low[,i],up = forecast_Up[,i],mean = forecast_mean[,i]))
      }
      p3 = ggplot2::ggplot(ggplot2::aes_(x = ~time,y = ~value),data = dat)
      p3 = p3 + ggplot2::geom_line() + ggplot2::theme_bw()
      p3 = p3 + ggplot2::geom_ribbon(ggplot2::aes_(ymin = ~low,ymax = ~up),fill = 'blue',alpha = 0.4)
      p3 = p3 + ggplot2::geom_line(ggplot2::aes_(x = ~time,y = ~mean),color = 'blue')
      p3 = p3 + ggplot2::labs(title = 'Covariates process and forecast')
      p3 = p3 + ggplot2::facet_grid(name~.,scales = 'free_y')
      plot2_list[['Xt']] = p3
    }
    return(plot2_list)
  }
}
# print = function(object, ...) UseMethod('print')
print.tsregime = function(object, ...){
  cat('Threshold time series:\n',
      'N =',object$N,
      'k =',ifelse(is.null(object$k),0,object$k),
      'nu = ',ifelse(is.null(object$nu),0,object$nu),'\n')
  dats = object
  class(dats) = NULL
  cat('======================','\n')
  if (!is.null(object$Yt)){
    if (is.na(sum(object$Yt))) {
      cat('Missing data in Yt:','\n')
      nas_Yt = apply(is.na(object$Yt),2,sum)
      names(nas_Yt) = paste0('Yt.',1:dim(object$Yt)[2])
      print(nas_Yt)
    }
  }
  if (!is.null(object$Zt)){
    if (is.na(sum(object$Zt))) {
      cat('Missing data in Zt:','\n')
      nas_Zt = apply(is.na(object$Zt),2,sum)
      names(nas_Zt) = paste0('Zt.',1:dim(object$Zt)[2])
      print(nas_Zt)
    }
  }
  if (!is.null(object$Xt)){
    if (is.na(sum(object$Xt))) {
      cat('Missing data in Xt:','\n')
      nas_Xt = apply(is.na(object$Xt),2,sum)
      names(nas_Xt) = paste0('Xt.',1:dim(object$Xt)[2])
      print(nas_Xt)
    }
  }
  if (!is.null(object$r)) {
    cat('======================','\n')
    cat('r = ',object$r,'\n')
    print(object$Summary_r)
  }else{
    if (!is.null(object$Zt)) {
      cat('Unknown threshold values','\n')
    }else{
      cat('Non-existent threshold variable','\n')
    }
  }
  cat('======================','\n')
  utils::str(dats)
}
print.regime_number = function(object, ...) {
  print(object$estimates)
  for (i in 1:length(object$list_m)) {
    cat(names(object$list_m)[[i]],'===============================|','\n')
    m_j = length(object$list_m[[i]]$orders$pj)
    nai_mj = c(object$NAIC[[i]]$AICj,object$NAIC[[i]]$NAIC)
    names(nai_mj) = c(paste0('AIC(R',1:m_j,')'),'NAIC')
    print(nai_mj,3)
    cat('\n')
    print(object$list_m[[i]]$par$estimates$r)
  }
}
print.regime_model = function(object, ...) {
  print(object$estimates)
}
print.regime_missing = function(object, ...) {
  print(object$estimates)
}
print.regime_forecast = function(object, ...) {
  print(object$forecast$estim)
}
