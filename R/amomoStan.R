#' A bayesian version of the A-MOMO algorithm
#'
#' Fits a Serfling -type model to a mortality time series by omiting summer and winter seasons from the estimation, and making predictions for the whole time span of the data. This version does not include delay adjustments
#'
#' @param data a data frame with all necessary variables
#' @param spring numbers of the weeks defining Spring (use 1:53 to use all weeks)
#' @param autumn numbers of the weeks defining Autumn (use 1:53 to use all weeks)
#' @param datevar name of the variable containing the date of the observation (preferrably the date of the Monday)
#' @param mortvar name of the variable containing the number of deaths
#' @param popvar (optionally) the name of the variable containing the population denominator
#' @param byvar (optionally) the name of the variable by which groups the analysis is done
#' @param iter number of iterations
#' @param penalties control the shrinkage between
#' @param data.only if TRUE just collect all the data together and don't fit
#' @return stan fit object
#' @export
amomoStan<-function(data,spring=16:25,autumn=31:47,datevar="date",mortvar="n",popvar=NA,byvar=NA,iter=2000,penalties=NULL,
                    data.only=FALSE) {
    ## Check the parameter
    if(!datevar%in%names(data)) stop("datevar not found")
    if(!mortvar%in%names(data)) stop("mortvar not found")
    if(is.na(popvar)) popvar<-"_const_"
    if(is.na(byvar)) byvar<-"_const_"
    if(!"_const_"%in%names(data)) data[[".const."]]<-rep(1,nrow(data))
    if(!popvar%in%names(data)) stop("popvar not found")
    if(!byvar%in%names(data)) stop("byvar not found")
    ## calculate the week and define the active period
    require(ISOweek)
    data[[".yw."]]  <-ISOweek::ISOweek(data[[datevar]])
    data[[".week."]]<-as.numeric(substring(data[[".yw."]],7))
    data[[".OK."]]  <-with(data,(.week.%in%spring |.week.%in%autumn))
    ## Input for Stan
    ## create response matrix, replace missing values with zeroes
    y   <-with(data,na.0(tapply(get(mortvar),list(get(datevar),get(byvar)),sum)))
    ## create population matrix, replace missing values with zeroes
    pop <-with(data,na.0(tapply(log(get(popvar)),list(get(datevar),get(byvar)),sum)))
    pop<-(pop-mean(pop))/pmax(1,diff(range(pop)))
    ## create indicator vector for inclusions
    OK <-with(data,na.0(tapply(.OK.,get(datevar),max))) # using max means we are inclusive
    ## create date vector based on the names
    date<-as.Date(dimnames(y)[[1]])
    ## create covariate matrix, shared by all by groups
    X<-cbind(1,momoTrend(date),momoSin(date)) # only one harmonic
    ## create new index to transform the data back to original order
    ndate<-c(date[OK==1],date[OK!=1])
    ## Penalties
    if(is.null(penalties)) penalties<-c(0.0,0.0,1.0)
    ## penalties<-pmax(penalties,0.00001)
    ## Create the data object for Stan
    standata<-list()
    ## Fit data
    standata$y<-y[OK==1,,drop=FALSE]
    standata$pop<-rbind(pop[OK==1,,drop=FALSE],pop[OK!=1,,drop=FALSE])
    standata$x<-rbind(X[OK==1,,drop=FALSE],X[OK!=1,,drop=FALSE])
    standata$N<-sum(OK==1)
    standata$M<-ncol(y)
    standata$NP<-sum(OK!=1)
    standata$NT<-nrow(pop)
    standata$P<-ncol(X)
    standata$penalty<-penalties
    ## Run the jewels
    if(!data.only) {
        require("rstan")
        ##print(names(stanmodels))
        ## eventually we will use this:
        ## print(system.time(fit <- try(rstan::sampling(stanmodels$amomo, data=standata,iter=1000, chains=4,verbose=TRUE))))
        ## a kludge
        print(system.time(fit <- try(rstan::stan(model_code=stanmodels$amomo@model_code,
                                                 data=standata,iter=iter, chains=4,verbose=FALSE))))
    } else {
        fit<-NULL
    }
    ## cleanup
    ## Reorganize the predictions so that you can compare with date
    ## MISSING
    ## Done
    return(list(data=data,standata=standata,fit=fit,OK=OK,date=date,ndate=ndate,
                vars=list(datevar=datevar,mortvar=mortvar,popvar=popvar,byvar=byvar)))
}


