.onLoad <- function(libname, pkgname) {
  modules <- paste0("stan_fit4", names(stanmodels), "_mod")
  if(length(modules)>0) {
    for (m in modules) {
      packageStartupMessage(m)
      #loadModule(m, what = TRUE)
    }
  } else {
    packageStartupMessage("FAIL")
  }
}
