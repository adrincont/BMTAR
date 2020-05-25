MTAR <img src="man/figures/logoMTAR.png" align="right" />
======================
Bayesian Analysis of Multivariate Threshold Autoregressive Models with Missing Data


![coverage](https://img.shields.io/badge/coverage-90%25-yellowgreen)
![version](https://img.shields.io/badge/version-0.1.0-blue)

The R package *MTAR* implements parameter estimation using a Bayesian approach for MTAR models with missing data using Markov Chain Monte Carlo methods. This package performs the simulation of MTAR processes ("mtarsim"), estimation of matrix parameters and the threshold values ("mtarns"), identification of the autoregressive orders using Bayesian variable selection ("mtarstr"), identification of the number of regimes using Metropolised Carlin and Chib ("mtarnumreg") and estimate missing data, coefficients and covariance matrices conditional on the autoregressive orders, the threshold values and the number of regimes  together ("mtarmissing").

## Installation
You can install the **development** version from [Github](https://github.com/adrincont/libreria-MTAR).
```s
install.packages("devtools")
library(devtools)
devtools::install_github("adrincont/libreria-MTAR")
```
## Example of use
```s
library(MTAR)
library(ggplot2)

data(datasim_miss)

data = tsregim(datasim_miss$Yt,datasim_miss$Zt,ddatasim_miss$Xt)
autoplot.tsregim(data,1)
autoplot.tsregim(data,2)
autoplot.tsregim(data,3)

Y_temp = t(datasim_miss$Yt)
meanY = apply(Y_temp,1,mean,na.rm = T)
Y_temp[apply(Y_temp,2,is.na)] = meanY
Y_temp = t(Y_temp)
X_temp = datasim_miss$Xt
meanX = mean(X_temp,na.rm = T)
X_temp[apply(X_temp,2,is.na)] = meanX
Z_temp = datasim_miss$Zt
meanZ = mean(Z_temp,na.rm = T)
Z_temp[apply(Z_temp,2,is.na)] = meanZ

data_temp = tsregim(Y_temp,Z_temp,X_temp)
initial = mtarinipars(tsregim_obj = data_temp,list_model = list(l0_max = 3),method = 'KUO')
estim_nr = mtarnumreg(ini_obj = initial,iterprev = 500,niter_m = 500,burn_m = 500, list_m = TRUE,
ordersprev = list(maxpj = 2,maxqj = 2,maxdj = 2))
print(estim_nr)

nitial = mtarinipars(tsregim_obj = data_temp,method = 'KUO',
list_model = list(pars = list(l = estim_nr$final_m),orders = list(pj = c(2,2))))
estruc = mtarstr(ini_obj = initial,niter = 500,chain = TRUE)
autoplot.regim_model(estruc,1)
autoplot.regim_model(estruc,2)
autoplot.regim_model(estruc,3)
autoplot.regim_model(estruc,4)
autoplot.regim_model(estruc,5)


initial = mtarinipars(tsregim_obj = data_temp,
list_model = list(pars = list(l = estim_nr$final_m,r = estruc$regim$r, orders = estruc$regim$orders)
missingest = mtarmissing(ini_obj = initial,chain = TRUE, niter = 500,burn = 500)
print(missingest)
autoplot.regim_missing(missingest,1)
data_c = missingest$tsregim

initial = mtarinipars(tsregim_obj = data_c,list_model = list(l0_max = 3),method = 'KUO')
estim_nr = mtarnumreg(ini_obj = initial,iterprev = 500,niter_m = 500,burn_m = 500, list_m = TRUE,
ordersprev = list(maxpj = 2,maxqj = 2,maxdj = 2))
print(estim_nr)

nitial = mtarinipars(tsregim_obj = data_c,method = 'KUO',
list_model = list(pars = list(l = estim_nr$final_m),orders = list(pj = c(2,2))))
estruc = mtarstr(ini_obj = initial,niter = 500,chain = TRUE)
autoplot.regim_model(estruc,1)
autoplot.regim_model(estruc,2)
autoplot.regim_model(estruc,3)
autoplot.regim_model(estruc,4)
autoplot.regim_model(estruc,5)

```
## For more information
You will find the theoretical basis of the method in the documents:
  * https://www.tandfonline.com/doi/abs/10.1080/03610926.2014.990758
  * https://core.ac.uk/download/pdf/77274943.pdf
## License
This package is free and open source software, licensed under GPL-3.

## References
 * Calderón Villanueva, S. A. (2014). Bayesian Analysis of Multivariate Threshold Autoregressive Models with Missing Data (Doctoral dissertation, Universidad Nacional de Colombia).

