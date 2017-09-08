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
get_cppcode(stanmodels$amomo)
loadModule("src/momoStan.so")
tmp<-amomoStan(fimomodata,byvar="age",popvar="pop")
tmp$fit2<-stan(file="exec/amomo.stan",data=tmp$standata)
str(tmp)
table(tmp$OK)
tmp$standata$pop
print(tmp$fit2,pars="alpha")

tmp$fit3<-sampling(stanmodels$amomo,tmp$standata,verbose=TRUE)
stanmodels$amomo@mk_cppmodule(stanmodels$amomo)
rstan:::grab_cxxfun(stanmodels$amomo@dso)
    
str(stanmodels$amomo)
