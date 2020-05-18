![coverage](https://img.shields.io/badge/coverage-60%25-yellowgreen)
![version](https://img.shields.io/badge/version-0-blue)
[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/ggplot2)](https://cran.r-project.org/package=ggplot2)
# libreria-MTAR

The R package *MTAR* implements parameter estimation using a Bayesian approach for MTAR models with missing data using Markov Chain Monte Carlo methods. This package performs the simulation of MTAR processes ("mtarsim"), estimation of matrix parameters and the threshold values ("mtarns"), identification of the autoregressive orders using Bayesian variable selection ("mtarstr"), identification of the number of regimes using Metropolised Carlin and Chib ("mtarnumreg") and estimate missing data, coefficients and covariance matrices conditional on the autoregressive orders, the threshold values and the number of regimes  together ("mtarmissing").

## Installation
You can install the **development** version from [Github](https://github.com/adrincont/libreria-MTAR).
```s
devtools::install_github("adrincont/libreria-MTAR")
```
## Example of use
```s
library(MTAR)
library(ggplot2)  
```
## For more information
You will find the theoretical basis of the method in the documents:
  * https://www.tandfonline.com/doi/abs/10.1080/03610926.2014.990758
  * https://core.ac.uk/download/pdf/77274943.pdf
## License
This package is free and open source software, licensed under GPL-3.

