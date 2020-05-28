# CRAN comments
## Test environments
* - using platform: x86_64-pc-linux-gnu (64-bit)
  - using R version 3.6.3 (2020-02-29)
  - using session charset: UTF-8
* - using platform: x86_64-w64-mingw32 (64-bit) 
  - using R version 3.6.0 (2019-04-26)
  - using session charset: ISO8859-1
* Rcloud
  - using platform: x86_64-pc-linux-gnu (64-bit)
  - using R version 3.6.0 (2019-04-26) 
  - using session charset: UTF-8
* - using platform: x86_64-apple-darwin15.6.0 (64-bit)
  - using R version 3.6.1 (2019-07-05)
  - using session charset: UTF-8

## R CMD check results
0 errors ✓ | 0 warnings ✓ | 0 notes ✓

## Time spending management
Because of the nature of MCMC methods, its time may depend on pc cores and its processers. As an example we show the duration of all examples were check() was done.

* x86_64-w64-mingw32 (64-bit)

v  checking examples (2h 34m 10.8s)
   Examples with CPU or elapsed time > 5s
   
   |                           |user |system| elapsed|
   | :---  | :---: | :---: | :---: |
   |auto_mtar               |3773.01 |200.39| 4195.75|
   |mtarnumreg              |2732.78 |150.72| 3026.69|
   |mtarmissing              |893.22  | 2.15| 1007.02|
   |mtarstr                  |683.47  |13.41 | 752.04|
   |mtarns                   |178.95|   6.59 | 222.50|
   |autoplot.regime_missing  | 2.11  | 0.44  | 10.59|
* x86_64-pc-linux-gnu (64-bit)

✓  checking examples (1h 17m 24.5s)
   Examples with CPU or elapsed time > 5s
   
   |                |user| system|  elapsed|
   | :---  | :---: | :---: | :---: |
   |auto_mtar   |2016.272  |5.551 |2025.287|
   |mtarnumreg  |1647.642 |37.849 |1688.025|
   |mtarmissing  |429.186  |0.383  |430.249|
   |mtarstr      |404.442  |9.430  |414.639|
   |mtarns        |73.629  |0.508   |74.271|
* Rcloud

✓  checking examples (58m 23.6s)
   Examples with CPU or elapsed time > 5s
   
   |                |user |system  |elapsed|
   | :---  | :---: | :---: | :---: |                
   |auto_mtar   |1504.099 |28.190 |1568.435|
   |mtarnumreg  |1160.001 |40.333 |1229.247|
   |mtarmissing  |318.512  |0.037  |327.572|
   |mtarstr      |291.368 |12.277  |310.413|
   |mtarns        |50.160  |1.196   |57.664|
   
* x86_64-apple-darwin15.6.0 (64-bit)

  checking examples ... OK
  Examples with CPU or elapsed time > 5s
  
   |               |user  |system  |elapsed|
   | :---  | :---: | :---: | :---: |               
   |auto_mtar   |2268.366 |326.850 |2673.075|
   |mtarnumreg  |1674.550 |311.636 |2033.274|
   |mtarstr      |415.585  |84.863  |511.333|
   |mtarmissing  |382.359   |1.117  |392.103|
   |mtarns        |69.469   |3.057   |75.231|

**Taking into account CRAN policies, the examples of auto_mtar,mtarnumreg,mtarstr and mtarmissing were commented, with this final change approximate final time was 3minutes.**


* This is a new release.
