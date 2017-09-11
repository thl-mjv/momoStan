###
### misc testing etc
###
### 1. setup the environment accoring to https://cran.r-project.org/web/packages/rstantools/vignettes/developer-guidelines.html
library("rstantools")
getwd()
rstan_package_skeleton("momoStan",code_files=list.files("R/",full=TRUE),stan_files=list.files(patt="*.stan"))
###
library("devtools")
document()
traceback()
build()
check()
###
load("data/fimomodata.RData")
summary(fimomodata)
library("ggplot2")
ggplot(fimomodata,aes(x=date,y=n/pop,group=age,colour=age))+geom_line()
ggplot(fimomodata,aes(x=date,y=TEMP,group=age,colour=age))+geom_line()
ggplot(fimomodata,aes(x=date,y=InfA,group=age,colour=age))+geom_line()
ggplot(fimomodata,aes(y=n/pop,x=TEMP,group=age,colour=age))+geom_point()+facet_wrap(~age,scales="free_y",ncol=1)
ggplot(fimomodata,aes(y=n/pop,x=InfA,group=age,colour=age))+geom_point()+facet_wrap(~age,scales="free_y",ncol=1)

### rehash the models
source("R/stanmodels.R")
source("R/utils.R")
source("R/amomoStan.R")
source("R/flumomoStan.R")
#loadModule("src/momoStan.so")
tmp<-amomoStan(fimomodata,byvar="age",popvar="pop",penalties=c(0,0,1))
tmp$standata$pop
print(tmp$fit,pars="shrinkage")
traceplot(tmp$fit,pars="shrinkage")
traceplot(tmp$fit,pars="alpha_param")
traceplot(tmp$fit,pars="alpha_group")
traceplot(tmp$fit,pars="alpha_real")
traceplot(tmp$fit,pars="alpha_eps")
pairs(tmp$fit,pars="shrinkage")

(foomatch<-with(tmp,match(date,ndate)))
(foo<-aperm(apply(extract(tmp$fit,pars="y_pred")[[1]],3:2,quantile,c(.5,0.025,.975)),c(3,1,2))[foomatch,,])
dim(foo)
par(mfcol=c(2,3))
for(i in 1:6) matplot(foo[,,i],type="l")
par(mfcol=c(1,1))
with(subset(tmp$data,age=="All"),plot(date,n))
matlines(tmp$date,foo[,,6])
with(tmp,plot(date,ndate))
with(tmp$data,plot(date,n,type="l"))

apply(apply(extract(tmp$fit,pars="shrinkage")[[1]],1,function(a) a[1:2]/sum(a)),1,quantile,c(0.5,0.025,0.975))

tmp$fit3<-sampling(stanmodels$amomo,tmp$standata,verbose=TRUE)
stanmodels$amomo@mk_cppmodule(stanmodels$amomo) # but why?
rstan:::grab_cxxfun(stanmodels$amomo@dso)

str(stanmodels$flumomo@model_code)
table(momoSeason(fimomodata$date))

ftmp<-flumomoStan(fimomodata,byvar="age",popvar="pop",penalties=c(0,0,1),covar=c("TEMP"),seasvar="InfA",data.only=TRUE)
ftmp$fit2<-stan(model_code=stanmodels$flumomo@model_code,data=ftmp$standata)

print(ftmp$fit2,pars="alpha_covariates")
traceplot(ftmp$fit2,pars="alpha_covariates")

(ffoo<-aperm(apply(extract(ftmp$fit,pars="pred_baseline")[[1]],3:2,quantile,c(.5,0.025,.975)),c(3,1,2)))
dim(ffoo)
par(mfcol=c(2,3))
for(i in 1:6) {
    matplot(ffoo[,,i],type="l")
}

ftmp$standata$z
traceback()
