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

tmp.pred<-predict(lm(TEMP~momoSin(date,3),data=fimomodata))
fimomodata$TEMP.pred<-tmp.pred
fimomodata$TEMP.hotter<-with(fimomodata,pmax(TEMP,TEMP.pred)-TEMP.pred) # above prediction
fimomodata$TEMP.hot   <-with(fimomodata,pmax(TEMP,max(TEMP.pred))-max(TEMP.pred)) # above hot
fimomodata$TEMP.hot2  <-with(fimomodata,ifelse(TEMP>max(TEMP.pred),max(TEMP.pred)-pmax(TEMP.pred,mean(TEMP.pred)),
                                        ifelse(TEMP<mean(TEMP.pred),0,ifelse(TEMP>TEMP.pred,TEMP-pmax(mean(TEMP.pred),TEMP.pred),0))))
fimomodata$TEMP.hot3  <-with(fimomodata,ifelse(TEMP>mean(TEMP.pred),mean(TEMP.pred)-pmin(TEMP.pred,mean(TEMP.pred)),
                                        ifelse(TEMP>TEMP.pred,TEMP-TEMP.pred,0)))
fimomodata$TEMP.hots  <-with(fimomodata,TEMP.hot+TEMP.hot2+TEMP.hot3)
fimomodata$TEMP.colder<-with(fimomodata,-pmin(TEMP,TEMP.pred)+TEMP.pred)
fimomodata$TEMP.cold  <-with(fimomodata,min(TEMP.pred)-pmin(TEMP,min(TEMP.pred)))
fimomodata$TEMP.cold2 <-with(fimomodata,ifelse(TEMP<min(TEMP.pred),-min(TEMP.pred)+pmin(TEMP.pred,mean(TEMP.pred)),
                                        ifelse(TEMP>mean(TEMP.pred),0,ifelse(TEMP<TEMP.pred,pmin(mean(TEMP.pred),TEMP.pred)-TEMP,0))))
fimomodata$TEMP.cold3 <-with(fimomodata,ifelse(TEMP<mean(TEMP.pred),-mean(TEMP.pred)+pmax(TEMP.pred,mean(TEMP.pred)),
                                        ifelse(TEMP<TEMP.pred,TEMP.pred-TEMP,0)))
fimomodata$TEMP.colds <-with(fimomodata,TEMP.cold+TEMP.cold2+TEMP.cold3)



par(mfcol=c(2,2))
for(i in 0:1) {
    with(subset(fimomodata,age=="All"),matplot(date,cbind(TEMP,TEMP.pred,mean(TEMP.pred),max(TEMP.pred),
                                                          TEMP.hot+i*max(TEMP.pred),
                                                          TEMP.hot2+i*pmax(TEMP.pred,mean(TEMP.pred)),
                                                          TEMP.hot3+i*TEMP.pred),type="l",
                                               lty=c(1,1,3,3,1,1,1),col=c(1,2,1,1,1,2,3),lwd=c(2,2,1,1,4,4,4)))
    with(subset(fimomodata,age=="All"),matplot(date,cbind(TEMP,TEMP.pred,mean(TEMP.pred),max(TEMP.pred),
                                                          -TEMP.cold+i*min(TEMP.pred),
                                                          -TEMP.cold2+i*pmin(TEMP.pred,mean(TEMP.pred)),
                                                          -TEMP.cold3+i*TEMP.pred),type="l",
                                               lty=c(1,1,3,3,1,1,1),col=c(1,2,1,1,1,2,3),lwd=c(2,2,1,1,4,4,4)))
}

with(subset(fimomodata,age=="All" & abs(TEMP.hots-TEMP.hotter)>0),cbind(TEMP,TEMP.pred,mean(tmp.pred),max(tmp.pred),TEMP.hotter,TEMP.hot,TEMP.hot2,TEMP.hot3,TEMP.hotter-TEMP.hots))

5.5979-c(3.7940, 4.668 )

source("R/stanmodels.R")
source("R/utils.R")
source("R/amomoStan.R")
source("R/flumomoStan.R")
ftmp<-flumomoStan(fimomodata,byvar="age",popvar="pop",covar=c("TEMP.hot","TEMP.hot2","TEMP.hot3","TEMP.cold","TEMP.cold3","TEMP.cold3"),
                  seasvar="InfA",data.only=FALSE)
ptmp<-flumomoStan(fimomodata,byvar="age",popvar="pop",covar=c("TEMP.hot","TEMP.hot2","TEMP.hot3","TEMP.cold","TEMP.cold3","TEMP.cold3"),
                  seasvar="InfA",data.only=FALSE,positive=FALSE)
str(ftmp)
if(inherts(ftmp$fit,"try-error"))
    ftmp$fit<-stan(model_code=stanmodels$flumomo@model_code,data=ftmp$standata)

print(ftmp$fit2,pars="alpha_covariates")
traceplot(ftmp$fit2,pars="alpha_covariates")

dim(bfoo<-aperm(apply(extract(ftmp$fit,pars="pred_full")[[1]],3:2,quantile,c(.5,0.025,.975)),c(3,1,2)))
dim(nfoo<-aperm(apply(extract(ftmp$fit,pars="pred_baseline_null")[[1]],3:2,quantile,c(.5,0.025,.975)),c(3,1,2)))
dim(ffoo<-aperm(apply(extract(ftmp$fit,pars="pred_baseline")[[1]],3:2,quantile,c(.5,0.025,.975)),c(3,1,2)))
dim(ffoo)
par(mfcol=c(2,3))
for(i in 1:5) {
    with(ftmp,plot(date,standata$y[,i],main=i,type="l"))
    matlines(ftmp$date,bfoo[,,i],type="l",col=4,lty=c(1,2,2),lwd=4)
    matlines(ftmp$date,nfoo[,,i],type="l",col=3,lty=c(1,2,2),lwd=3)
    matlines(ftmp$date,ffoo[,,i],type="l",col=2,lty=c(1,2,2),lwd=2)
    matlines(ftmp$date, foo[,,i],type="l",col=1,lty=c(1,2,2),lwd=1)
}

ftmp$standata$z
traceback()

head(b1<-apply(extract(ftmp$fit,pars="pred_full")[[1]],2:3,mean))
head(b2<-apply(extract(ftmp$fit,pars="pred_baseline")[[1]],2:3,mean))
head(b3<-apply(extract(ftmp$fit,pars="pred_baseline_null")[[1]],2:3,mean))
#dim(b4<-apply(extract(ftmp$fit,pars="pred_effects")[[1]],2:4,mean))
head(b5<-apply(extract( tmp$fit,pars="y_pred")[[1]],2:3,mean)[foomatch,])

par(mfcol=c(1,1))
matplot(cbind(ftmp$standata$y[,5],b1[,5],b2[,5],b3[,5],b5[,5]),type="l",lwd=3,lty=1)
matplot(b1,type="l")
matlines(b2,type="l")
matlines(b3,type="l")
matlines(b5,type="l")

matplot(t(apply(b6,1,cumsum)),type="l",lty=1,col=1)
lines(b1[,5]-b2[,5],col="red")
dimnames(ftmp$standata$z)[[3]]
stan_dens(ftmp$fit,pars=paste0("alpha_covariates[",1:13,",5]"))
traceplot(ftmp$fit,pars=paste0("alpha_covariates[",1:14,",5]"))
pairs(ftmp$fit,pars=paste0("alpha_covariates[",1:14,",5]"))
pairs(ptmp$fit,pars=paste0("alpha_covariates[",1:13,",5]"))
pairs(ftmp$fit,pars=paste0("alpha_baseline[",1:4,",5]"))
matplot(c4<-apply(extract(ftmp$fit,pars="alpha_covariates")[[1]][,,5],2,sort),type="l")

addcol<-function(a,b=1) a+b*col(a)
par(mfcol=c(1,1))
matplot(addcol(ftmp$standata$z[,5,]),type="l")
apply(100*ftmp$standata$z[,5,],2,max)
apply(100*ftmp$standata$z[,5,],2,min)

dim(a4<-extract(ftmp$fit,pars="pred_baseline")[[1]][,,5])
dim(b4<-extract(ptmp$fit,pars="pred_baseline")[[1]][,,5])
dim(c4<-extract(ftmp$fit,pars="alpha_covariates")[[1]][,,5])
dim(d4<-extract(ptmp$fit,pars="alpha_covariates")[[1]][,,5])
dim(e4<-extract(ftmp$fit,pars="pred_full")[[1]][,,5])
dim(f4<-extract(ptmp$fit,pars="pred_full")[[1]][,,5])
dim(z4<-t(ftmp$standata$z[,5,]))
dim(c4%*%z4)
dimnames(z4)

par(mfcol=c(2,2))
matplot(t(apply(a4*(-1+exp(c4[, 1:6]%*%z4[ 1:6,])),2,quantile,c(.5,.025,.975))),type="l")
matplot(t(apply(a4*(-1+exp(c4[,7:13]%*%z4[7:13,])),2,quantile,c(.5,.025,.975))),type="l")
matplot(t(apply(b4*(-1+exp(d4[, 1:6]%*%z4[ 1:6,])),2,quantile,c(.5,.025,.975))),type="l")
matplot(t(apply(b4*(-1+exp(d4[,7:13]%*%z4[7:13,])),2,quantile,c(.5,.025,.975))),type="l")

matplot(cbind(apply(a4*(exp(c4%*%z4))-e4,2,mean),apply(b4*(exp(d4%*%z4))-f4,2,mean)),type="l",lty=1)
