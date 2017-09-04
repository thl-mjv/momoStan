###
### misc testing etc
###
### 1. setup the environment accoring to https://cran.r-project.org/web/packages/rstantools/vignettes/developer-guidelines.html
library("rstantools")
getwd()
rstan_package_skeleton("momoStan",code_files=list.files("R/",full=TRUE),stan_files=list.files(patt="*.stan"))
