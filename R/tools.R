#' Extract interesting information from a fit
#'
#' @param obj fit object
#' @param temp either numeric vector or a regular expression to identify temperature effects
#' @param infl either numeric vector or a regular expression to identify influenxa effects
#' @return a list with components
#' @export
momoEffects<-function(obj,temp="TEMP",infl="InfA") {
    if(is.null(obj$standata$z)) {
        Qp<-obj$standata$Qp
        Qf<-obj$standata$Qf
        Z<-with(obj$standata,array(NA,dim=c(N,M,Qp+Qf),
                                   dimnames=list(dimnames(zp)[[1]],dimnames(zp)[[2]],
                                                 c(dimnames(zf)[[3]],dimnames(zp)[[3]]))))
        Z[,,   1:Qf]<-obj$standata$zf
        Z[,,Qf+1:Qp]<-obj$standata$zp
        print(dim(alphap<-extract(obj$fit,pars="alpha_covariates_plus")[[1]]))
        print(dim(alphaf<-extract(obj$fit,pars="alpha_covariates_free")[[1]]))
        alpha<-array(NA,dim=c(dim(alphap)[1],Qp+Qf,dim(alphap)[3]))
        alpha[,   1:Qf,]<-alphaf
        alpha[,Qf+1:Qp,]<-alphap
    } else {
        Z<-obj$standata$z
        alpha<-extract(obj$fit,pars="alpha_covariates")[[1]]
    }
    M<-dim(Z)[2] # number of subgroups
    print(dim(Z))
    print(dim(alpha))
    baseline<-extract(obj$fit,pars="pred_baseline")[[1]]
    full    <-extract(obj$fit,pars="pred_full"    )[[1]]
    if(is.character(temp)) temp<-grep(temp,covar)
    if(is.character(infl)) infl<-grep(infl,covar)
    print(temp)
    print(infl)
    tempeff<-lapply(1:M,function(m) alpha[,temp,m]%*%t(Z[,m,temp]))
    infleff<-lapply(1:M,function(m) alpha[,infl,m]%*%t(Z[,m,infl]))
    temptotal<-lapply(1:M,function(m) {
        apply(apply(baseline[,,m]*(exp(tempeff[[m]])-1),1,tapply,momoSeason(obj$date),sum),1,quantile,c(0.5,0.025,0.975))
    })
    infltotal<-lapply(1:M,function(m) {
        apply(apply(baseline[,,m]*(exp(infleff[[m]])-1),1,tapply,momoSeason(obj$date),sum),1,quantile,c(0.5,0.025,0.975))
    })
    list(Z=Z,alpha=alpha,baseline=baseline,full=full,tempeff=tempeff,infleff=infleff,
         temptotal=temptotal,infltotal=infltotal)
}
