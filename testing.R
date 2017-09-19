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
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

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

apply(apply(extract(tmp$fit,pars="shrinkage")[[1]],1,function(a) a[1:2]/sum(a)),1,quantile,c(0.5,0.025,0.975))

### Fix the data for further analyses
smalldata<-fixMomoData()
summary(smalldata)


### Larger model
source("R/stanmodels.R")
source("R/utils.R")
source("R/amomoStan.R")
source("R/flumomoStan.R")
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

atmp<-amomoStan(smallmomodata,byvar="age",popvar="pop",penalties=c(0,0,1))

ptmp<-flumomoStan(smalldata,byvar="age",popvar="pop",covar=c("TEMP.hot","TEMP.hot2","TEMP.hot3","TEMP.cold","TEMP.cold3","TEMP.cold3"),
                  seasvar="InfA",data.only=FALSE,positive="all" ,iter=1000) # 560s
ftmp<-flumomoStan(smalldata,byvar="age",popvar="pop",covar=c("TEMP.hot","TEMP.hot2","TEMP.hot3","TEMP.cold","TEMP.cold3","TEMP.cold3"),
                  seasvar="InfA",data.only=FALSE,positive="none",iter=1000) #351s 615
mtmp<-flumomoStan(smalldata,byvar="age",popvar="pop",covar=c("TEMP.hot","TEMP.hot2","TEMP.hot3","TEMP.cold","TEMP.cold3","TEMP.cold3"),
                  seasvar="InfA",data.only=FALSE,positive="some",iter=1000,positives="InfA") # 1193
Mtmp<-flumomoStan(smalldata,byvar="age",popvar="pop",covar=c("TEMP.hot","TEMP.hot2","TEMP.hot3","TEMP.cold","TEMP.cold3","TEMP.cold3"),
                  seasvar="InfA",data.only=FALSE,positive="some",iter=1000,positives="TEMP") # 1047
system.time(save(tmp,ptmp,ftmp,mtmp,Mtmp,file="tmp.rda")) #464s

names(tmp$positives)
names(ftmp$positives)

print(ftmp$fit,pars="alpha_covariates")
pairs(ptmp$fit,pars="alpha_baseline")
stan_dens(ftmp$fit,pars="alpha_baseline")
stan_dens(mtmp$fit,pars=c("alpha_covariates_free","alpha_covariates_plus"),separate_chains = TRUE)+xlim(c(-.4,.4))+
  geom_vline(aes(xintercept=0))
stan_dens(Mtmp$fit,pars=c("alpha_covariates_plus","alpha_covariates_free"),separate_chains = TRUE)+xlim(c(-.4,.4))+
  geom_vline(aes(xintercept=0))

source("R/tools.R")
system.time(aeff<-momoEffects(atmp,temp=NULL,infl=NULL))
system.time(feff<-momoEffects(ftmp,temp=1:6,infl=7:12)) # 20s
system.time(peff<-momoEffects(ptmp,temp=1:6,infl=7:12)) # 21
system.time(meff<-momoEffects(mtmp,temp=1:6,infl=7:12)) # 10s
system.time(Meff<-momoEffects(Mtmp,temp=7:12,infl=1:6)) # 10s

infleff<-rbind(momoStanPlotData(peff,type="infleff",modelname="restricted"),
               momoStanPlotData(feff,type="infleff",modelname="free"))#,
#               momoStanPlotData(meff,type="infleff",modelname="temp free, infl pos"),
#               momoStanPlotData(Meff,type="infleff",modelname="temp pos, infl free"))

tempeff<-rbind(momoStanPlotData(peff,type="tempeff",modelname="restricted"),
               momoStanPlotData(feff,type="tempeff",modelname="free"))#,
#               momoStanPlotData(meff,type="tempeff",modelname="temp free, infl pos"),
#               momoStanPlotData(Meff,type="tempeff",modelname="temp pos, infl free"))

baselines<-rbind(momoStanPlotData(aeff,type="baseline",modelname="a-momo(like)"),
                 momoStanPlotData(peff,type="baseline",modelname="restricted"),
                 momoStanPlotData(feff,type="baseline",modelname="free"))#,
#                 momoStanPlotData(meff,type="baseline",modelname="temp free, infl pos"),
#                 momoStanPlotData(Meff,type="baseline",modelname="temp pos, infl free"))
baselines$yy<-with(subset(baselines,modelname=="free"),y[match(baselines$x,x)])
summary(baselines)
pdf("restrict.pdf",paper="a4r",width=10,height=8)
ggplot(infleff,aes(x=x))+geom_line(aes(y=Y),alpha=.5)+geom_ribbon(aes(ymin=Ymin,ymax=Ymax),alpha=.25)+facet_wrap(~modelname,ncol=1)+theme_minimal()+
  labs(y="Excees",x="Time",title="Influenza")
ggplot(tempeff,aes(x=x))+geom_line(aes(y=Y),alpha=.5)+geom_ribbon(aes(ymin=Ymin,ymax=Ymax),alpha=.25)+facet_wrap(~modelname,ncol=1)+theme_minimal()+
  labs(y="Excess",x="Time",title="Temperature")
ggplot(baselines,aes(x=x))+geom_line(aes(y=yy),colour="blue")+
  geom_line(aes(y=Y),alpha=.75)+geom_ribbon(aes(ymin=Ymin,ymax=Ymax),alpha=.25,colour=NA)+facet_wrap(~modelname,ncol=1)+theme_minimal()+
  labs(y="Deaths",x="Time",title="Baselines")
ggplot(baselines,aes(x=x,group=modelname,colour=modelname))+geom_line(aes(y=yy),colour="black",alpha=.25)+
  geom_line(aes(y=Y),size=2,alpha=.75)+theme_minimal()+
  labs(y="Deaths",x="Time",title="Baselines compared")
dev.off()
