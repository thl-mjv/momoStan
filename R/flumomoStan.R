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
#' @param iter number of iterations
#' @param positive which covariates should be restricted to positives, either all, none or some
#' @param positives if positive=="some", which. Treated as list of regular expressions, anything that matches is treated as a hit
#' @param data.only if TRUE just collect all the data together and don't fit
#' @return stan fit object
#' @export
flumomoStan<-function(data,spring=16:25,autumn=31:47,datevar="date",mortvar="n",
                      covar=NULL,seasvar=NULL,
                      popvar=NA,byvar=NA,iter=2000,
                      positive=c("all","none","some"),positives=NULL,
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
    Zok<-rep(TRUE,length(covar))
    names(Zok)<-covar
    print(dim(Z))
    for(i in covar) {
        cat(i," ")
        z<-with(res$data,na.0(tapply(get(i),list(get(res$vars$datevar),get(res$vars$byvar)),sum)))
        if(any(z<0)) cat("covariates should be non-negative ")
        if(!any(z>0)) {
            cat("no positive values ")
            Zok[i]<-FALSE
        } else {
            z<-z/pmax(1,max(abs(z))) # scale!
        }
        cat(dim(z),range(z))
        Z[,,i]<-z
        cat("\n")
    }
    covar<-covar[Zok]
    ## print(str(res))
    res$standata$z<-Z[,,Zok,drop=FALSE]
    res$standata$Q<-sum(Zok)
    res$standata$NT<-res$standata$penalty<-res$standata$NP<-NULL
    ## Run the jewels
    if(!data.only) {
        require("rstan")
        positive<-match.arg(positive)
        if(positive=="some") {
            cat("Positives=some:")
            model<-stanmodels$flumomomixed
            npos<-0
            if(length(positives)>0) {
                posZ<-apply(sapply(positives,function(a) grepl(a,covar)),1,any)
                npos<-sum(posZ)
                cat(npos,"/",sum(!posZ)," of ",length(posZ),"...",sep="")
            }
            if(npos==0) {
                cat(" No covariate set as positive!")
                positive<-"none"
            } else {
                if(npos==length(covar)) {
                    cat(" All covariates set as positive!")
                    positive<-"all"
                } else {
                    cat(npos,"covariates set as positive!")
                    Zplus<-res$standata$z[,, posZ,drop=FALSE]
                    Zfree<-res$standata$z[,,!posZ,drop=FALSE]
                    res$standata$z<-NULL
                    res$standata$zp<-Zplus
                    res$standata$zf<-Zfree
                    res$standata$Qp<-sum( posZ)
                    res$standata$Qf<-sum(!posZ)
                    cat(dim(Zplus),dim(Zfree))
                    positives<-covar[posZ]
                    freepars<-covar[!posZ]
                }
            }
        }
        if(positive=="all" ) {
            model<-stanmodels$flumomopositive
            positives<-covar
            freepars<-NULL
        }
        if(positive=="none") {
            model<-stanmodels$flumomofree
            positives<-NULL
            freepars<-covar
        }
        ##print(names(stanmodels))
        ## eventually we will use this:
        print(system.time(fit <- try(rstan::sampling(stanmodels$amomo, data=res$standata,
                                                     iter=iter, chains=4,verbose=TRUE))))
        ## a kludge
        if(inherits(fit,"try-error")) {
            cat("Recalculating\n")
            print(system.time(fit <- try(rstan::stan(model_code=model@model_code,
                                                     data=res$standata,iter=iter, chains=4,verbose=FALSE))))
        }
    } else {
        fit<-NULL
    }
    res$fit<-fit
    res$vars<-c(res$vars,list(covar=covar,seasvar=seasvar))
    res$positives<-list(positives=positives,free=freepars)
    cat("Done\n")
    res
}



