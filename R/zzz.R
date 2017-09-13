#' @useDynLib momoStan, .registration = TRUE
#' 
.onLoad <- function(libname, pkgname) {
  modules <- paste0("stan_fit4", names(stanmodels), "_mod")
  if(length(modules)>0) {
    for (m in modules) {
    # try(loadModule(m, what = TRUE))
      packageStartupMessage(m)
    }
    packageStartupMessage("SUCCESS")
  } else {
    packageStartupMessage("FAIL")
  }
}

