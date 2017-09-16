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
matlines(tmp$date,foo[,,5])
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

dim(smallmomodata<-subset(fimomodata,age=="All" & momoSeason(date)%in%2011:2016))
with(smallmomodata,table(momoSeason(date)))

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
summary(smallmomodata)
5.5979-c(3.7940, 4.668 )

source("R/stanmodels.R")
source("R/utils.R")
source("R/amomoStan.R")
source("R/flumomoStan.R")
ptmp<-flumomoStan(smallmomodata,byvar="age",popvar="pop",covar=c("TEMP.hot","TEMP.hot2","TEMP.hot3","TEMP.cold","TEMP.cold3","TEMP.cold3"),
                  seasvar="InfA",data.only=FALSE,positive="all" ,iter=1000) # 560s
ftmp<-flumomoStan(smallmomodata,byvar="age",popvar="pop",covar=c("TEMP.hot","TEMP.hot2","TEMP.hot3","TEMP.cold","TEMP.cold3","TEMP.cold3"),
                  seasvar="InfA",data.only=FALSE,positive="none",iter=1000) #351s 615
mtmp<-flumomoStan(smallmomodata,byvar="age",popvar="pop",covar=c("TEMP.hot","TEMP.hot2","TEMP.hot3","TEMP.cold","TEMP.cold3","TEMP.cold3"),
                  seasvar="InfA",data.only=FALSE,positive="some",iter=1000,positives="InfA") # 1193
Mtmp<-flumomoStan(smallmomodata,byvar="age",popvar="pop",covar=c("TEMP.hot","TEMP.hot2","TEMP.hot3","TEMP.cold","TEMP.cold3","TEMP.cold3"),
                  seasvar="InfA",data.only=FALSE,positive="some",iter=1000,positives="TEMP") # 1047
system.time(save(ptmp,ftmp,mtmp,Mtmp,file="tmp.rda")) #464s
if(inherts(ftmp$fit,"try-error"))
    ftmp$fit<-stan(model_code=stanmodels$flumomo@model_code,data=ftmp$standata)

str(ptmp)

print(ftmp$fit,pars="alpha_covariates")
pairs(ptmp$fit,pars="alpha_baseline")
stan_dens(ftmp$fit,pars="alpha_baseline")
stan_dens(mtmp$fit,pars=c("alpha_covariates_free","alpha_covariates_plus"),separate_chains = TRUE)+xlim(c(-.4,.4))+
  geom_vline(aes(xintercept=0))
stan_dens(Mtmp$fit,pars=c("alpha_covariates_plus","alpha_covariates_free"),separate_chains = TRUE)+xlim(c(-.4,.4))+
  geom_vline(aes(xintercept=0))

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

source("R/tools.R")
system.time(feff<-momoEffects(ftmp,temp=1:6,infl=7:12)) # 20s
system.time(peff<-momoEffects(ptmp,temp=1:6,infl=7:12)) # 21
system.time(meff<-momoEffects(mtmp,temp=1:6,infl=7:12)) # 10s
system.time(Meff<-momoEffects(Mtmp,temp=8:12,infl=1:7)) # 10s

p1<-plot(peff,type="infleff")
momoStanFeature(p1,feff,type="infleff",colour="blue")


par(mfrow=c(2,2))
matplot (2011:2016,t(feff$temptotal[[1]]),type="l",lty=c(1,2,2),lwd=c(4,1,1),col=1,ylim=c(-1000,2000))
matlines(2011:2016,t(peff$temptotal[[1]]),type="l",lty=c(1,2,2),lwd=c(4,1,1),col=2)
matlines(2011:2016,t(meff$temptotal[[1]]),type="l",lty=c(1,2,2),lwd=c(4,1,1),col=3)
matlines(2011:2016,t(Meff$temptotal[[1]]),type="l",lty=c(1,2,2),lwd=c(4,1,1),col=4)
abline(h=0)
matplot (2011:2016,t(feff$infltotal[[1]]),type="l",lty=c(1,2,2),lwd=c(4,1,1),col=1,ylim=c(-1000,2000))
matlines(2011:2016,t(peff$infltotal[[1]]),type="l",lty=c(1,2,2),lwd=c(4,1,1),col=2)
matlines(2011:2016,t(meff$infltotal[[1]]),type="l",lty=c(1,2,2),lwd=c(4,1,1),col=3)
matlines(2011:2016,t(Meff$infltotal[[1]]),type="l",lty=c(1,2,2),lwd=c(4,1,1),col=4)
abline(h=0)
plot(ftmp$standata$y[,1],type="l")
matlines(t(apply(feff$baseline[,,1],2,quantile,c(0.5,0.025,0.975))),lty=c(1,2,2),lwd=c(4,1,1),col=1)
matlines(t(apply(peff$baseline[,,1],2,quantile,c(0.5,0.025,0.975))),lty=c(1,2,2),lwd=c(4,1,1),col=2)
matlines(t(apply(meff$baseline[,,1],2,quantile,c(0.5,0.025,0.975))),lty=c(1,2,2),lwd=c(4,1,1),col=3)
matlines(t(apply(Meff$baseline[,,1],2,quantile,c(0.5,0.025,0.975))),lty=c(1,2,2),lwd=c(4,1,1),col=4)
plot(ftmp$standata$y[,1],type="l")
matlines(t(apply(feff$full    [,,1],2,quantile,c(0.5,0.025,0.975))),lty=c(1,2,2),lwd=c(4,1,1),col=1)
matlines(t(apply(peff$full    [,,1],2,quantile,c(0.5,0.025,0.975))),lty=c(1,2,2),lwd=c(4,1,1),col=2)
matlines(t(apply(meff$full    [,,1],2,quantile,c(0.5,0.025,0.975))),lty=c(1,2,2),lwd=c(4,1,1),col=3)
matlines(t(apply(Meff$full    [,,1],2,quantile,c(0.5,0.025,0.975))),lty=c(1,2,2),lwd=c(4,1,1),col=4)

100*apply(feff$alpha>0,2:3,mean)
100*apply(peff$alpha>0,2:3,mean)
100*apply(meff$alpha>0,2:3,mean)
100*apply(Meff$alpha>0,2:3,mean)

qnts<-c(.5,0.025,0.975,.25,.75)
par(mfcol=c(2,5))
for(i in 1:5) {
matplot (t(apply((ftmp$standata$y[,i])-t(feff$full[,,i]),1,function(a) quantile((a-mean(a))/sd(a),qnts))),type="l",lty=c(1,2,2),col=1,main=i)
matlines(t(apply((ptmp$standata$y[,i])-t(peff$full[,,i]),1,function(a) quantile((a-mean(a))/sd(a),qnts))),type="l",lty=c(1,2,2),col=2)
matlines(t(apply((mtmp$standata$y[,i])-t(meff$full[,,i]),1,function(a) quantile((a-mean(a))/sd(a),qnts))),type="l",lty=c(1,2,2),col=3)
matlines(t(apply((Mtmp$standata$y[,i])-t(Meff$full[,,i]),1,function(a) quantile((a-mean(a))/sd(a),qnts))),type="l",lty=c(1,2,2),col=4)
matplot (t(apply((ftmp$standata$y[,i])-t(feff$baseline[,,i]),1,function(a) quantile((a-mean(a))/sd(a),qnts))),type="l",lty=c(1,2,2),col=1)
matlines(t(apply((ptmp$standata$y[,i])-t(peff$baseline[,,i]),1,function(a) quantile((a-mean(a))/sd(a),qnts))),type="l",lty=c(1,2,2),col=2)
matlines(t(apply((mtmp$standata$y[,i])-t(meff$baseline[,,i]),1,function(a) quantile((a-mean(a))/sd(a),qnts))),type="l",lty=c(1,2,2),col=3)
matlines(t(apply((Mtmp$standata$y[,i])-t(Meff$baseline[,,i]),1,function(a) quantile((a-mean(a))/sd(a),qnts))),type="l",lty=c(1,2,2),col=4)
}
meff$temptotal
str(Mtmp$standata)

par(mfcol=c(1,1))
plot((apply((ftmp$standata$y[,i])-t(feff$full[,,i]),1,function(a) mean((a-mean(a))/sd(a)>2))))
plot(apply((ftmp$standata$y[,i])-t(feff$baseline[,,i]),1,function(a) mean(a>0)))
abline(h=c(0.05,0.95))

