#' Tool to fix the data in analyses
#'
#' @param data a dataset to be transformed
#' @param keep subgroups to keep in the variable keepvar
#' @param keepvar which var to use in keep
#' @param seas which seasons to use
#' @param datevar which variable contains date
#' @param tempvar variable with the temperature indicator
#' @param tempsin number of fourier components
#' @return data frame suitable to analysis
#' @export
fixMomoData<-function(data=NULL,keep=NULL,keepvar="age",seas=NULL,datevar="date",tempvar="TEMP",tempsin=3) {
  if(is.null(data)) {
    data("fimomodata")
    data<-fimomodata
  }
  if(is.null(keep)) keep<-unique(data[[keepvar]]) # fixme to include any other variable
  if(is.null(seas)) seas<-unique(momoSeason(data[[datevar]]))
  form<-paste(tempvar,"~momoSin(",datevar,",",tempsin,")")
  tmp.pred<-predict(lm(formula(form),data=data))
  if(!".TEMP"%in%names(data)) data$.TEMP<-data[[tempvar]]
  data$.TEMP.pred<-tmp.pred
  data$.TEMP.hotter<-with(data,pmax(.TEMP,.TEMP.pred)-.TEMP.pred) # above prediction
  data$.TEMP.hot   <-with(data,pmax(.TEMP,max(.TEMP.pred))-max(.TEMP.pred)) # above hot
  data$.TEMP.hot2  <-with(data,ifelse(.TEMP>max(.TEMP.pred),max(.TEMP.pred)-pmax(.TEMP.pred,mean(.TEMP.pred)),
                                                 ifelse(.TEMP<mean(.TEMP.pred),0,ifelse(.TEMP>.TEMP.pred,.TEMP-pmax(mean(.TEMP.pred),.TEMP.pred),0))))
  data$.TEMP.hot3  <-with(data,ifelse(.TEMP>mean(.TEMP.pred),mean(.TEMP.pred)-pmin(.TEMP.pred,mean(.TEMP.pred)),
                                                 ifelse(.TEMP>.TEMP.pred,.TEMP-.TEMP.pred,0)))
  data$.TEMP.hots  <-with(data,.TEMP.hot+.TEMP.hot2+.TEMP.hot3)
  data$.TEMP.colder<-with(data,-pmin(.TEMP,.TEMP.pred)+.TEMP.pred)
  data$.TEMP.cold  <-with(data,min(.TEMP.pred)-pmin(.TEMP,min(.TEMP.pred)))
  data$.TEMP.cold2 <-with(data,ifelse(.TEMP<min(.TEMP.pred),-min(.TEMP.pred)+pmin(.TEMP.pred,mean(.TEMP.pred)),
                                                 ifelse(.TEMP>mean(.TEMP.pred),0,ifelse(.TEMP<.TEMP.pred,pmin(mean(.TEMP.pred),.TEMP.pred)-.TEMP,0))))
  data$.TEMP.cold3 <-with(data,ifelse(.TEMP<mean(.TEMP.pred),-mean(.TEMP.pred)+pmax(.TEMP.pred,mean(.TEMP.pred)),
                                                 ifelse(.TEMP<.TEMP.pred,.TEMP.pred-.TEMP,0)))
  data$.TEMP.colds <-with(data,.TEMP.cold+.TEMP.cold2+.TEMP.cold3)

  res<-subset(data,get(keepvar)%in%keep & momoSeason(get(datevar))%in%seas)
}
