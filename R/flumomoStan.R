#' A bayesian version of the FluMOMO algorithm
#'
#' Fits a regression model to a mortality time series using any number of covariates with effects constrained to be positive.
#'
#' @param data a data frame with all necessary variables
#' @param datevar name of the variable containing the date of the observation (preferrably the date of the Monday)
#' @param mortvar name of the variable containing the number of deaths
#' @param covar character vector of covariates included in the model as such
#' @param seasvar character vector of covariates included in the model by season
#' @param popvar (optionally) the name of the variable containing the population denominator
#' @param byvar (optionally) the name of the variable by which groups the analysis is done
#' @param penalties control the shrinkage between
#' @param data.only if TRUE just collect all the data together and don't fit
#' @return stan fit object
#' @export
flumomoStan<-function(data,spring=16:25,autumn=31:47,datevar="date",mortvar="n",
                      covar=NULL,seasvar=NULL,
                      popvar=NA,byvar=NA,penalties=NULL,
                      data.only=FALSE) {
    res<-amomoStan(data,spring=1:53,autumn=1:53,datevar=datevar,mortvar=mortvar,
                   popvar=popvar,byvar=byvar,penalties=c(0,0,1),
                   data.only=TRUE)
    ##    print(str(res))
    if(is.null(covar) & is.null(seasvar)) stop("Must have at least one variable")
    res$data[[".season."]]<-momoSeason(res$data[[datevar]])
    if(!is.null(seasvar)) {
        seasons<-sort(unique(res$data[[".season."]]))
        for(i in seasons) {
            cat(i,":")
            for(j in seasvar) {
                cat(j,"...")
                if(j%in%names(res$data)) {
                    res$data[[paste0(j,i)]]<-with(res$data,ifelse(.season.==i,get(j),0))
                    covar<-c(covar,paste0(j,i))
                } else {
                    cat("! ")
                }
            }
            cat("\n")
        }
    }
    covar<-covar[covar%in%names(res$data)]
    Z<-with(res$standata,array(0,dim=c(N,M,length(covar)),dimnames=list(1:N,1:M,covar)))
    print(dim(Z))
    for(i in covar) {
        cat(i)
        z<-with(res$data,na.0(tapply(get(i),list(get(res$vars$datevar),get(res$vars$byvar)),sum)))
        z<-z/pmax(1,max(abs(z))) # scale!
        cat(dim(z))
        Z[,,i]<-z
        cat("\n")
    }
    ## print(str(res))
    res$standata$z<-Z
    res$standata$Q<-length(covar)
    res$standata$NT<-res$standata$penalty<-res$standata$NP<-NULL
    ## Run the jewels
    if(!data.only) {
        require("rstan")
        ##print(names(stanmodels))
        ## eventually we will use this:
        ## print(system.time(fit <- try(rstan::sampling(stanmodels$amomo, data=standata,iter=1000, chains=4,verbose=TRUE))))
        ## a kludge
        print(system.time(fit <- try(rstan::stan(model_code=stanmodels$flumomo@model_code,
                                                 data=res$standata,iter=2000, chains=4,verbose=FALSE))))
    } else {
        fit<-NULL
    }
    res$fit<-fit
    res$vars<-c(res$vars,covar=covar,seasvar=seasvar)
    res
}



