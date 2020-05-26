MTAR <img src="man/figures/logoMTAR.png" align="right" />
======================
Bayesian Analysis of Multivariate Threshold Autoregressive Models with Missing Data


![coverage](https://img.shields.io/badge/coverage-90%25-yellowgreen)
![version](https://img.shields.io/badge/version-0.1.0-blue)

The R package *MTAR* implements parameter estimation using a Bayesian approach for MTAR models with missing data using Markov Chain Monte Carlo methods. This package performs the simulation of MTAR process (`mtarsim`). Estimation of matrix parameters and the threshold values conditional on the autoregressive orders and number of regimes (`mtarns`). Identification of the autoregressive orders using Bayesian variable selection, together with coefficients and covariance matrices and the threshold values conditional on the number of regimes (`mtarstr`). Identification of the number of regimes using Metropolised Carlin and Chib or via NAIC criteria (`mtarnumreg`), to calculate NAIC of any estimated model (`mtarNAIC`). Estimate missing values together with matrix parameters conditional to threshold values, autoregressive orders and numbers of regimes (`mtarmissing`). The diagnostic of the residuals in any estimated model can be done (`diagnostic_mtar`). The package manage several class objects for autoplot and print, functions like (`tsregime`),(`mtaregime`) and (`mtarinipars`) make its construction. Finally, (`auto_mtar`) its an automatic function that perfoms all above.

## Installation
You can install the **development** version from [Github](https://github.com/adrincont/MTAR).
```s
install.packages("devtools")
devtools::install_github("adrincont/mtar")
```

## Overview
As mention in the first paragraph lets introduce the objects class and usage in the different functions.

- `tsregime` return an object class 'tsregime' which is how the package manage data.
- `mtaregime` return an object of class 'regime' use for simulation purposes and as standard presentation of the final estimations.
- `mtarsim` return an object of class 'mtarsim' use in autoplot methods. Its practical to conditionate some functions for different known parameters.
- `mtarinipars` return an object of class 'regime_inipars' that itself contains an object of class 'tsregime', it is the main object that save known parameters and parameters of the prior distributions for each parameter in a MTAR model. This object needs to be provided in every estimation function.
- `mtarns` and `mtarstr` return an object of class 'regime_model' use in print and autoplot methods, its an stardard presentation for estimations done in this functions. It is the object to introduce in `mtarNAIC`. 
- `mtarmissing` return an object of class 'regime_missing' for print and autoplot methods.
- `mtarnumreg` return an object of class 'regime_number'

## Example of use
```s
library(mtar)
library(ggplot2)

data(datasim_miss)

data = tsregime(datasim_miss$Yt,datasim_miss$Zt,datasim_miss$Xt)
autoplot.tsregime(data,1)
autoplot.tsregime(data,2)
autoplot.tsregime(data,3)

# Fill in the missing data with the component average
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

# Estimate the number of regimens with the completed series
data_temp = tsregime(Y_temp,Z_temp,X_temp)
initial = mtarinipars(tsregime_obj = data_temp,list_model = list(l0_max = 3),method = 'KUO')
estim_nr = mtarnumreg(ini_obj = initial,iterprev = 500,niter_m = 500,burn_m = 500, list_m = TRUE,ordersprev = list(maxpj = 2,maxqj = 2,maxdj = 2),parallel = TRUE)
print(estim_nr)

# Estimate the structural and non-structural parameters 
# for the series once we know the number of regimes and some idea of its orders
initial = mtarinipars(tsregime_obj = data_temp,method = 'KUO',
                      list_model = list(pars = list(l = estim_nr$final_m),
                      orders = list(pj = c(2,2))))
estruc = mtarstr(ini_obj = initial,niter = 500,chain = TRUE)
autoplot.regime_model(estruc,1)
autoplot.regime_model(estruc,2)
autoplot.regime_model(estruc,3)
autoplot.regime_model(estruc,4)
autoplot.regime_model(estruc,5)
cc = 1/2*(nrow(data$Yt)-2)-1
## in the table of critical values for CusumSQ we have a value of 
## 0.05333 for n = 500 and alpha = 0.05
diagnostic_mtar(estruc,CusumSQ = 0.05333)

# With the known structural parameters we estimate the missing data
list_model = list(pars = list(l = estim_nr$final_m,r = estruc$estimates$r[,2],orders = estruc$orders))
initial = mtarinipars(tsregime_obj = data_temp,list_model = list_model)
missingest = mtarmissing(ini_obj = initial,chain = TRUE, niter = 500,burn = 500)
print(missingest)
autoplot.regime_missing(missingest,1)
data_c = missingest$tsregim
# ============================================================================================#
# Once the missing data has been estimated, we make the estimates again for all the structural 
# and non-structural parameters.
# ============================================================================================#
initial = mtarinipars(tsregime_obj = data_c,list_model = list(l0_max = 3),method = 'KUO')
estim_nr = mtarnumreg(ini_obj = initial,iterprev = 500,niter_m = 500,burn_m = 500, list_m = TRUE,ordersprev = list(maxpj = 2,maxqj = 2,maxdj = 2))
print(estim_nr)

initial = mtarinipars(tsregime_obj = data_c,method = 'KUO',
list_model = list(pars = list(l = estim_nr$final_m),orders = list(pj = c(2,2))))
estruc = mtarstr(ini_obj = initial,niter = 500,chain = TRUE)
autoplot.regime_model(estruc,1)
autoplot.regime_model(estruc,2)
autoplot.regime_model(estruc,3)
autoplot.regime_model(estruc,4)
autoplot.regime_model(estruc,5)
diagnostic_mtar(estruc,CusumSQ = 0.05333)
```

## Other useful examples
MTAR is a general model were it is possible to specificated other kind of models we are familiar with, like:

- Basic auto-regressive model AR(p): 
- Vector auto-regressive model VAR(p): 
- Threshold auto-regressive model TAR(l,p,q): 

| spec/Model | AR | VAR | TAR |
| :---  | :---: | :---: | :---: |
| k | 1 | >= 1 | 1 |
| Regimes | 1 | 1 | > 1 |
| Threshold process   | x | x | ✓ |

This can be useful when you have missing data in one of this types of models and use **mtar package** for its estimation based on a bayesian approach.

- AR
```s
library(mtar)
library(ggplot2)
library(forecast)
# AR = MTAR k = 1, l = 1, Zt = NO
R1 = mtaregime(orders = list(p = 2),Phi = list(phi1 = 0.4,phi2 = 0.3),Sigma = 2)
data = mtarsim(100,list(R1))
ardata = arima.sim(list(ar = c(0.4,0.3),sd = 2),100)
ggpubr::ggarrange(
autoplot(tsregime(ardata)) + ggplot2::labs(title = 'base package'),
autoplot(data$Sim) + ggplot2::labs(title = 'mtar package'),ncol = 2)
arima1 = arima(ts(data$Sim$Yt),c(2,0,0))
parameters = list(l = 1,orders = list(pj = 2))
initial = mtarinipars(tsregime_obj = data$Sim,list_model = list(pars = parameters))
estim1 = mtarns(ini_obj = initial,niter = 1000,chain = TRUE)
print.regime_model(estim1)
ggpubr::ggarrange(
autoplot(estim1,5) + theme(legend.position = 'none') + 
labs(title = 'mtar package'),
ggplot(data = NULL,aes(x = 1:100,y = data$Sim$Yt)) + 
geom_line(col = 'black') + geom_line(data = NULL,
aes(x = 1:100,y = fitted(arima1)),col = "blue") + theme_bw() + 
labs(title = 'forecast package'),ncol = 2)
diagnostic_mtar(estim1)
```

- VAR
```s
library(mtar)
library(ggplot2)
# VAR = MTAR k > 1, l = 1, Zt = NO
library(vars)
library(BVAR)
R1 = mtaregime(orders = list(p = 1,q = 0,d = 0),
              Phi = list(phi1 = matrix(c(0.3,0.2,0.1,0.4),2,2)),Sigma = matrix(c(1,0.5,0.5,1),2,2))
data = mtarsim(100,list(R1))
autoplot(data$Sim)
var1 = vars::VAR(y = data$Sim$Yt,p = 1)
BVAR::bvar(data = data$Sim$Yt,lags = 1)
parameters = list(l = 1,orders = list(pj = 1))
initial = mtarinipars(tsregim_obj = data$Sim,list_model = list(pars = parameters))
estim1 = mtarns(ini_obj = initial,niter = 1000,chain = TRUE)
estim1$regime
print.regime_model(estim1)
autoplot(estim1,5)
diagnostic_mtar(estim1)
```

- TAR
```s
# Example 1, TAR model with 2 regimes
Z = arima.sim(n = 500,list(ar = c(0.5)))
l = 2;r = 0;K = c(2,1)
theta = matrix(c(1,-0.5,0.5,-0.7,-0.3,NA), nrow = l)
H = c(1, 1.5)
X = simu.tar.norm(Z,l,r,K,theta,H)
Yt = tsregim(Yt = X,Zt = Z,r = r)
R1 = mtaregim(orders = list(p = 2),cs = 1,Phi = list(phi1 = -0.5,phi2 = 0.5),
              Sigma = 1)
R2 = mtaregim(orders = list(p = 1),cs = -0.7,Phi = list(phi1 = -0.3),
              Sigma = sqrt(1.5))
YtSim = mtarsim(500,list(R1,R2),r,Zt = Z)
ggpubr::ggarrange(
autoplot(Yt) + ggplot2::labs(title = 'TAR package'),
autoplot(YtSim$Sim) + ggplot2::labs(title = 'mtar package'),ncol = 2)
# number of regimes
res = reg.thr.norm(Z,X)
res$L.est
res$L.prob
res$R.est
res$R.CI
initial = mtarinipars(Yt,list_model = list(l0_min = 2,l0_max = 3),method = 'KUO')
resmtar = mtarnumreg(initial)
# structural parameters
res2 = ARorder.norm(Z,X,l,r)
res2$K.est
res2$K.prob
initial = mtarinipars(Yt,list_model = list(pars = list(l = 2),
orders = list(pj = c(2,2),dj = c(1,1))),method = 'KUO')
res2mtar = mtarstr(initial)
res2mtar$orders
# non-structural parameters
res3 = Param.norm(Z,X,l,r,K) #gibbs
res4 = LS.norm(Z,X,l,r,c(0,0)) #least square
initial = mtarinipars(Yt,list(pars = list(l = 2,orders = list(pj = c(1,1)))))
res3mtar = mtarns(initial)
```
## For more information
You will find the theoretical basis of the method in the documents:

  - https://www.tandfonline.com/doi/abs/10.1080/03610926.2014.990758
  - https://core.ac.uk/download/pdf/77274943.pdf

## License
This package is free and open source software, licensed under GPL-3.

## References
 * Calderón Villanueva, S. A. (2014). Bayesian Analysis of Multivariate Threshold Autoregressive Models with Missing Data (Doctoral dissertation, Universidad Nacional de Colombia).

