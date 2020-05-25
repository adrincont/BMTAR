## Test environments
* ─  using platform: x86_64-pc-linux-gnu (64-bit)
  ─  using R version 3.6.3 (2020-02-29)
  ─  using session charset: UTF-8
* windows x64. R.3.6.0
* Rcloud
  - using R version 3.5.3 (2019-03-11)
  - using platform: x86_64-pc-linux-gnu (64-bit)
  - using session charset: UTF-8
* local OS X install, R 3.6.3

## R CMD check results
0 errors ✓ | 0 warnings ✓ | 0 notes ✓

## Time spending 
Because of the nature of MCMC methods, to run all examples when checking its done on Rcloud was: Duration: 1h 35s.

## Downstream dependencies
I have also run R CMD check on downstream dependencies of httr 
(https://github.com/wch/checkresults/blob/master/httr/r-release). 
All packages that I could install passed except:

* Ecoengine: this appears to be a failure related to config on 
  that machine. I couldn't reproduce it locally, and it doesn't 
  seem to be related to changes in httr (the same problem exists 
  with httr 0.4).
