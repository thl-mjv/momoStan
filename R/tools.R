#' Extract interesting information from a fit
#'
#' @param obj fit object
#' @param temp either numeric vector or a regular expression to identify temperature effects
#' @param infl either numeric vector or a regular expression to identify influenxa effects
#' @param quants which quantiles to calculate for the effects
#' @return a list with components ...
#' @export
momoEffects<-function(obj,temp="TEMP",infl="InfA",quants=c(0.5,0.025,0.975)) {
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
    tempeff<-lapply(1:M,function(m) baseline[,,m]*(exp(alpha[,temp,m]%*%t(Z[,m,temp]))-1))
    infleff<-lapply(1:M,function(m) baseline[,,m]*(exp(alpha[,infl,m]%*%t(Z[,m,infl]))-1))
    temptotal<-lapply(1:M,function(m) {
        apply(apply(tempeff[[m]],1,tapply,momoSeason(obj$date),sum),1,quantile,quants)
    })
    infltotal<-lapply(1:M,function(m) {
        apply(apply(infleff[[m]],1,tapply,momoSeason(obj$date),sum),1,quantile,quants)
    })
    res<-list(
      date=obj$date,
      y=obj$standata$y,
      Z=Z,
      quants=quants,
      alpha   =apply(alpha   ,2:3,quantile,quants),
      baseline=apply(baseline,2:3,quantile,quants),
      full    =apply(full    ,2:3,quantile,quants),
      tempeff=lapply(tempeff,apply,2,quantile, quants),
      infleff=lapply(infleff,apply,2,quantile, quants),
      temptotal=temptotal,infltotal=infltotal)
    class(res)<-"momoStanEffs"
    return(res)
}
#' Plot methods for momoStanEffs
#'
#' @param obj an object to plot
#' @param type type of the plot
#' @param colour colour of the ribbon
#' @param alpha transparency of the ribbon
#' @param group which group to be plotted
#' @param plt a plot to be added to
#' @return a data or a ggplot
#' @export
momoStanPlotData<-function(obj,type="baseline",colour="red",group=1,alpha=.2) {
  if(!type%in%names(obj)) stop("Unknown type")
  if(type%in%c("baseline","full"))
    data<-with(obj,data.frame(x=date,
                              y=y[,group],
                              Y=get(type)[1,,group],
                              Ymin=get(type)[2,,group],
                              Ymax=get(type)[3,,group]))
  if(type%in%c("tempeff","infleff"))
    data<-with(obj,data.frame(x=date,y=0,
                              Y=get(type)[[group]][1,],Ymin=get(type)[[group]][2,],Ymax=get(type)[[group]][3,]))
  return(data)
}
#' @describeIn momoStanPlotData additional features to a plot
#' @export
momoStanFeature<-function(plt,obj,type="baseline",colour="red",group=1,alpha=.2) {
  data<-momoStanPlotData(obj,type,group=group)
  plt+geom_line(data=data,aes(y=Y),size=2,colour=colour)+
    geom_ribbon(data=data,aes(ymin=Ymin,ymax=Ymax),colour=NA,fill=colour,alpha=alpha)
}
#' @describeIn momoStanPlotData plot data from momoStanEffs
#' @export
plot.momoStanEffs<-function(obj,type="baseline",colour="red",group=1,alpha=.2) {
  data<-momoStanPlotData(obj,type,group=group)
  ggplot(data,aes(x=x))+geom_line(aes(y=y))+geom_line(aes(y=Y),size=2,colour=colour)+
    geom_ribbon(aes(ymin=Ymin,ymax=Ymax),colour=NA,fill=colour,alpha=alpha)
}
