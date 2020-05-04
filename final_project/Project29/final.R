rm(list = ls())

require(doParallel)
cl <- makeCluster(25)# The number of cores on this machine 
registerDoParallel(cl)
mcopts <- list(set.seed=TRUE)



## ----setup,echo=F,results=F,cache=F-------------------------------------------
myround<- function (x, digits = 1) {
  # taken from the broman package
  if (digits < 1) 
    stop("This is intended for the case digits >= 1.")
  if (length(digits) > 1) {
    digits <- digits[1]
    warning("Using only digits[1]")
  }
  tmp <- sprintf(paste("%.", digits, "f", sep = ""), x)
  zero <- paste0("0.", paste(rep("0", digits), collapse = ""))
  tmp[tmp == paste0("-", zero)] <- zero
  tmp
}

set.seed(2050320976)

options(
  keep.source=TRUE,
  encoding="UTF-8"
)



## ----prelims,echo=F,cache=F---------------------------------------------------
set.seed(594709947L)
library(ggplot2)
theme_set(theme_bw())
library(tidyverse)
library(plyr)
library(reshape2)
library(foreach)
#library(doMC)
library(pomp)
stopifnot(packageVersion("pomp")>="2.0")


## ----aapl_plot_code,echo=T,eval=F--------------------------------------------
## dat <- read.table("aapl.csv",sep=",",header=TRUE)
## plot(as.Date(dat$Date),dat$Close,
##   xlab="date",ylab="S&P 500",type="l")
## plot(as.Date(dat$Date),dat$Close, log="y",
##   xlab="date",ylab="S&P 500",type="l")


## ----aapl_plot,out.width="10cm",echo=F---------------------------------------
par(mfrow=c(1,2),mai=c(0.8,0.8,0.1,0.3))
dat <- read.table("aapl.csv",sep=",",header=TRUE)
plot(as.Date(dat$Date),dat$Close,
     xlab="date",ylab="S&P 500",type="l")
plot(as.Date(dat$Date),dat$Close, log="y",
     xlab="date",ylab="S&P 500",type="l")


## ----data_rda_code,echo=T,eval=F----------------------------------------------
load(file="aapl-2002-2012.rda")
plot(aapl.ret.demeaned, type="l",xlab="business day (2002-2012)",ylab="demeaned S&P 500 return")


## ----data_rda_plot,out.width="8cm",echo=F-------------------------------------
par(mai=c(0.8,0.8,0.1,0.1))
load(file="aapl-2002-2012.rda")
plot(aapl.ret.demeaned, type="l",
     xlab="business day (2002-2012)",ylab="demeaned S&P 500 return")


## ----garch--------------------------------------------------------------------
require(tseries)
fit.garch <- garch(aapl.ret.demeaned,grad = "numerical",
                   trace = FALSE)
L.garch <- tseries:::logLik.garch(fit.garch)


## ----names--------------------------------------------------------------------
aapl_statenames <- c("H","G","Y_state")
aapl_rp_names <- c("sigma_nu","mu_h","phi","sigma_eta")
aapl_ivp_names <- c("G_0","H_0")
aapl_paramnames <- c(aapl_rp_names,aapl_ivp_names)


## ----rproc--------------------------------------------------------------------
rproc1 <- "
double beta,omega,nu;
omega = rnorm(0,sigma_eta * sqrt( 1- phi*phi ) * 
sqrt(1-tanh(G)*tanh(G)));
nu = rnorm(0, sigma_nu);
G += nu;
beta = Y_state * sigma_eta * sqrt( 1- phi*phi );
H = mu_h*(1 - phi) + phi*H + beta * tanh( G ) 
* exp(-H/2) + omega;
"
rproc2.sim <- "
Y_state = rnorm( 0,exp(H/2) );
"

rproc2.filt <- "
Y_state = covaryt;
"
aapl_rproc.sim <- paste(rproc1,rproc2.sim)
aapl_rproc.filt <- paste(rproc1,rproc2.filt)


## ----rinit--------------------------------------------------------------------
aapl_rinit <- "
G = G_0;
H = H_0;
Y_state = rnorm( 0,exp(H/2) );
"


## ----measure------------------------------------------------------------------
aapl_rmeasure <- "
y=Y_state;
"

aapl_dmeasure <- "
lik=dnorm(y,0,exp(H/2),give_log);
"


## ----transforms---------------------------------------------------------------
aapl_partrans <- parameter_trans(
  log=c("sigma_eta","sigma_nu"),
  logit="phi"
)


## ----sp_pomp------------------------------------------------------------------
aapl.filt <- pomp(data=data.frame(
  y=aapl.ret.demeaned,time=1:length(aapl.ret.demeaned)),
  statenames=aapl_statenames,
  paramnames=aapl_paramnames,
  times="time",
  t0=0,
  covar=covariate_table(
    time=0:length(aapl.ret.demeaned),
    covaryt=c(0,aapl.ret.demeaned),
    times="time"),
  rmeasure=Csnippet(aapl_rmeasure),
  dmeasure=Csnippet(aapl_dmeasure),
  rprocess=discrete_time(step.fun=Csnippet(aapl_rproc.filt),
                         delta.t=1),
  rinit=Csnippet(aapl_rinit),
  partrans=aapl_partrans
)


## ----sim_pomp-----------------------------------------------------------------
params_test <- c(
  sigma_nu = exp(-4.5),  
  mu_h = -0.25,  	 
  phi = expit(4),	 
  sigma_eta = exp(-0.07),
  G_0 = 0,
  H_0=0
)

sim1.sim <- pomp(aapl.filt, 
                 statenames=aapl_statenames,
                 paramnames=aapl_paramnames,
                 rprocess=discrete_time(
                   step.fun=Csnippet(aapl_rproc.sim),delta.t=1)
)

sim1.sim <- simulate(sim1.sim,seed=1,params=params_test)


## ----build_sim1.filt----------------------------------------------------------

sim1.filt <- pomp(sim1.sim, 
                  covar=covariate_table(
                    time=c(timezero(sim1.sim),time(sim1.sim)),
                    covaryt=c(obs(sim1.sim),NA),
                    times="time"),
                  statenames=aapl_statenames,
                  paramnames=aapl_paramnames,
                  rprocess=discrete_time(
                    step.fun=Csnippet(aapl_rproc.filt),delta.t=1)
)




## ----run_level----------------------------------------------------------------
run_level <- 3
aapl_Np <-           switch(run_level, 100, 1e3, 2e3)
aapl_Nmif <-         switch(run_level,  10, 100, 200)
aapl_Nreps_eval <-   switch(run_level,   4,  10,  20)
aapl_Nreps_local <-  switch(run_level,  10,  20,  20)
aapl_Nreps_global <- switch(run_level,  10,  20, 100)


## ----parallel-setup,cache=FALSE-----------------------------------------------
library(doParallel)
registerDoParallel()
library(doRNG)
registerDoRNG(34118892)

## ----pf1----------------------------------------------------------------------
stew(file=sprintf("pf1-%d.rda",run_level),{
  t.pf1 <- system.time(
    pf1 <- foreach(i=1:aapl_Nreps_eval,
                   .packages='pomp') %dopar% pfilter(sim1.filt,Np=aapl_Np))
},seed=493536993,kind="L'Ecuyer")
(L.pf1 <- logmeanexp(sapply(pf1,logLik),se=TRUE))


## ----mif_setup----------------------------------------------------------------
aapl_rw.sd_rp <- 0.02
aapl_rw.sd_ivp <- 0.1
aapl_cooling.fraction.50 <- 0.5
aapl_rw.sd <- rw.sd(
  sigma_nu  = aapl_rw.sd_rp,
  mu_h      = aapl_rw.sd_rp,
  phi       = aapl_rw.sd_rp,
  sigma_eta = aapl_rw.sd_rp,
  G_0       = ivp(aapl_rw.sd_ivp),
  H_0       = ivp(aapl_rw.sd_ivp)
)	 


## ----mif----------------------------------------------------------------------
stew(file=sprintf("mif1-%d.rda",run_level),{
  t.if1 <- system.time({
    if1 <- foreach(i=1:aapl_Nreps_local,
                   .packages='pomp', .combine=c) %dopar% mif2(aapl.filt,
                                                              params=params_test,
                                                              Np=aapl_Np,
                                                              Nmif=aapl_Nmif,
                                                              cooling.fraction.50=aapl_cooling.fraction.50,
                                                              rw.sd = aapl_rw.sd)
    L.if1 <- foreach(i=1:aapl_Nreps_local,
                     .packages='pomp', .combine=rbind) %dopar% logmeanexp(
                       replicate(aapl_Nreps_eval, logLik(pfilter(aapl.filt,
                                                                  params=coef(if1[[i]]),Np=aapl_Np))), se=TRUE)
  })
},seed=318817883,kind="L'Ecuyer")

r.if1 <- data.frame(logLik=L.if1[,1],logLik_se=L.if1[,2],
                    t(sapply(if1,coef)))
if (run_level>1) write.table(r.if1,file="aapl_params.csv",
                             append=TRUE,col.names=FALSE,row.names=FALSE)


## ----pairs_code,echo=T,eval=F-------------------------------------------------
## pairs(~logLik+sigma_nu+mu_h+phi+sigma_eta,
##   data=subset(r.if1,logLik>max(logLik)-20))


## ----pairs_plot,echo=F,eval=T,out.width="11cm"--------------------------------
pairs(~logLik+sigma_nu+mu_h+phi+sigma_eta,
      data=subset(r.if1,logLik>max(logLik)-20))


## ----box----------------------------------------------------------------------
aapl_box <- rbind(
  sigma_nu=c(0.005,0.05),
  mu_h    =c(-1,0),
  phi = c(0.95,0.99),
  sigma_eta = c(0.5,1),
  G_0 = c(-2,2),
  H_0 = c(-1,1)
)


## ----box_eval-----------------------------------------------------------------
stew(file=sprintf("box_eval-%d.rda",run_level),{
  t.box <- system.time({
    if.box <- foreach(i=1:aapl_Nreps_global,
                      .packages='pomp',.combine=c) %dopar% mif2(if1[[1]],
                                                                params=apply(aapl_box,1,function(x)runif(1,x)))
    L.box <- foreach(i=1:aapl_Nreps_global,
                     .packages='pomp',.combine=rbind) %dopar% {
                       logmeanexp(replicate(aapl_Nreps_eval, logLik(pfilter(
                         aapl.filt,params=coef(if.box[[i]]),Np=aapl_Np))), 
                         se=TRUE)
                     }
  })
},seed=290860873,kind="L'Ecuyer")

r.box <- data.frame(logLik=L.box[,1],logLik_se=L.box[,2],
                    t(sapply(if.box,coef)))
if(run_level>1) write.table(r.box,file="aapl_params.csv",
                            append=TRUE,col.names=FALSE,row.names=FALSE)
summary(r.box$logLik,digits=5)


## ----pairs_global_code,eval=F,echo=T------------------------------------------
## pairs(~logLik+log(sigma_nu)+mu_h+phi+sigma_eta+H_0,
##   data=subset(r.box,logLik>max(logLik)-10))


## ----pairs_global_plot,eval=T,echo=F,out.width="11cm"-------------------------
pairs(~logLik+log(sigma_nu)+mu_h+phi+sigma_eta+H_0,
      data=subset(r.box,logLik>max(logLik)-10))


## ----garch_benchmark,eval=F,echo=F--------------------------------------------
## library(tseries)
## fit.garch.benchmark <- garch(aapl.ret.demeaned,grad = "numerical", trace = FALSE)
## L.garch.benchmark <- tseries:::logLik.garch(fit.garch.benchmark)

