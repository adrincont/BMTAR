
autoplot.mtarsim = function(object, ci=0.95, ...) {
  if (!requireNamespace('ggplot2', quietly = TRUE)) {
    stop('ggplot2 is needed for this function to work. Install it via install.packages(\'ggplot2\')', call. = FALSE)
  }
  else {
    if (!inherits(object, 'mtarsim')) {
      stop('autoplot.mtarsim requires a mtarsim object, use object=object')
    }}
  dats = as.data.frame(object$Yt)
  time = seq(1,nrow(dats))
  dat = data.frame(name = 'Series.1',time = time,value = dats[,1])
  for (i in 2:ncol(object$Yt)) {
    dat = rbind(dat,data.frame(name = paste0('Series.',i),time = time,value = dats[,i]))
  }
  p = ggplot2::ggplot(data = dat,ggplot2::aes(x = time,y = value))
  p + ggplot2::geom_line() + ggplot2::facet_grid(name~.)
}