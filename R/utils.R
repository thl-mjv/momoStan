#' Date of Monday
#'
#' @param date A Date object to convert
#' @return a Date of the last Monday
#' @export
momoMonday<-function(date) {
    w<-format(date,"%w")
    date-ifelse(w==0,6,w-1)
}
#' A calendar of weeks
#'
#' Create a dataset of weeks with additional metadata. Useful for extending the observed data to incorporate future dates for forecasting
#' @param start first week to consider (a Date rounded to a Monday)
#' @param end last week to consider (a Date rounded to a Monday)
#' @return a data frame with fields date ...
#' @export
momoCalendar<-function(start,end) {
    ## missing
}
#' Sinusoidal covariates
#'
#' Create n x 2 covariates with varying wave lengths
#' @param date a Date variable to base the waves on
#' @param n number of wavelengths
#' @param freq base frequency of the longest wave in years. Default 1 means longest wave is yearly
#' @param use which frequencies to use. Default is 1 to n
#' @return a matrix of components
#' @export
momoSin<-function(date,n=1,freq=NA,use=NULL) {
  if(is.na(freq)) freq<-1 
  x<-as.POSIXlt(date)
  a<-x$year+x$yday/365.25
  if(is.null(use))
       use<-1:n
  n<-length(use)
  res<-matrix(NA,nrow=length(x),ncol=2*n)
  dimnames(res)<-list(as.character(x),paste(rep(c("Sin","Cos"),c(n,n)),rep(use,2)))
  for(i in use) res[,paste("Sin",i)]<-sin(2*pi*i*a/(freq))
  for(i in use) res[,paste("Cos",i)]<-cos(2*pi*i*a/(freq))
  res
}
#' A trend effect from a date
#'
#' Turns a Date variable into a trend suitable for modeling. Currently only linear trend is available
#'
#' @param date a Date variable
#' @param origin origin date, defaults to the midpoint
#' @return a numeric vector with values between -1 and 1
#' @export
momoTrend<-function(date,origin=NULL) {
    dmin<-as.numeric(min(date,na.rm=TRUE))
    dmax<-as.numeric(max(date,na.rm=TRUE))
    if(is.null(origin)) origin<-(dmax+dmin)/2
    as.numeric(date-origin)/(dmax-dmin)
}
#' Replace missing values with zeroes
#'
#' @param a variable to be transformed
#' @param b value that used to replace missing values, defaults 0. Recycled
#' @return transformed variable
na.0<-function(a,b=0) ifelse(is.na(a),b,a)
#' Calculate the season from a date
#'
#' Define season as a year that starts given week. If this week is 1, season corresponds to a normal year
#'
#' @param date sate to use
#' @param start the week when year starts
#' @param fact if TRUE, return a factor, otherwise return a number
#' @return number giving the start year of the season or a factor
#' @export
momoSeason<-function(date,start=26,fact=FALSE) {
    require("ISOweek")
    yw<-ISOweek(date)
    year<-as.numeric(substring(yw,1,4))
    week<-as.numeric(substring(yw,7))
    season<-ifelse(week<start,year-1,year)
    if(fact) {
        yrs<-seq(min(season),max(season))
        lab<-ifelse(start>1,paste(yrs,yrs+1,sep="-"),paste(yrs))
        return(factor(season),levels=yrs,labels=lab)
    } else {
        return(season)
    }
}
