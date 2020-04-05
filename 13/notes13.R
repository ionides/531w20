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


## ----data---------------------------------------------------------------------
polio_data <- read.csv("polio_wisconsin.csv",comment="#")
head(polio_data,4)


## ----data_plot,echo=F,out.width="12cm"----------------------------------------
polio_data %>%
  gather(variable,value,-time) %>%
  ggplot(aes(x=time,y=value))+
  geom_line()+
  facet_wrap(~variable,ncol=1,scales='free_y',strip.position = "left")+
  theme_bw()+
  theme(
    strip.background=element_rect(fill=NA,color=NA),
    strip.placement="outside"
  )+
  labs(x="",y="")


## ----statenames---------------------------------------------------------------
polio_statenames <- c("SB1","SB2","SB3","SB4","SB5","SB6",
  "IB","SO","IO")
polio_obsnames <- "cases"
polio_t0 <- 1932+4/12


## ----covariates---------------------------------------------------------------
polio_K <- 6
polio_covar <- covariate_table(
  t=polio_data$time,
  B=polio_data$births,
  P=predict(smooth.spline(x=1931:1954,
    y=polio_data$pop[seq(12,24*12,by=12)]))$y,
  periodic.bspline.basis(t,nbasis=polio_K,
    degree=3,period=1,names="xi%d"),
  times="t"
)


## ----rp_names-----------------------------------------------------------------
polio_rp_names <- c("b1","b2","b3","b4","b5","b6",
  "psi","rho","tau","sigma_dem","sigma_env")


## ----ivp_names----------------------------------------------------------------
polio_ivp_names <- c("SO_0","IO_0")
polio_paramnames <- c(polio_rp_names,polio_ivp_names)


## ----fixed_names--------------------------------------------------------------
polio_fp_names <- c("delta","K",
  "SB1_0","SB2_0","SB3_0","SB4_0","SB5_0","SB6_0")
polio_paramnames <- c(polio_rp_names,
  polio_ivp_names,polio_fp_names)
covar_index_t0 <- which(abs(polio_covar@times-polio_t0)<0.01)
polio_initial_births <- polio_covar@table["B",covar_index_t0-0:5]
names(polio_initial_births) <- c("SB1_0","SB2_0",
  "SB3_0","SB4_0","SB5_0","SB6_0") 
polio_fixed_params <- c(delta=1/60,K=polio_K,
  polio_initial_births)


## ----polio_read_mle-----------------------------------------------------------
polio_params_guess <- c(b1=3,b2=0,b3=1.5,b4=6,b5=5,b6=3,
  psi=0.002,rho=0.01,tau=0.001,sigma_dem=0.04,sigma_env=0.5,
  SO_0=0.12,IO_0=0.001,polio_fixed_params)


## ----rprocess-----------------------------------------------------------------
polio_rprocess <- Csnippet("
  double beta = exp(dot_product( (int) K, &xi1, &b1));
  double lambda = (beta * (IO+IB) / P + psi);
  double var_epsilon = pow(sigma_dem,2)/ lambda +  
    pow(sigma_env,2);
  lambda *= (var_epsilon < 1.0e-6) ? 1 : 
    rgamma(1/var_epsilon,var_epsilon);
  double p = exp(- (delta+lambda)/12);
  double q = (1-p)*lambda/(delta+lambda);
  SB1 = B;
  SB2= SB1*p;
  SB3=SB2*p;
  SB4=SB3*p;
  SB5=SB4*p;
  SB6=SB5*p;
  SO= (SB6+SO)*p;
  IB=(SB1+SB2+SB3+SB4+SB5+SB6)*q;
  IO=SO*q;
")


## ----measure------------------------------------------------------------------
polio_dmeasure <- Csnippet("
  double tol = 1.0e-25;
  double mean_cases = rho*IO;
  double sd_cases = sqrt(pow(tau*IO,2) + mean_cases);
  if(cases > 0.0){
    lik = pnorm(cases+0.5,mean_cases,sd_cases,1,0)
      - pnorm(cases-0.5,mean_cases,sd_cases,1,0) + tol; 
  } else{
    lik = pnorm(cases+0.5,mean_cases,sd_cases,1,0) + tol;
  }
  if (give_log) lik = log(lik);
")

polio_rmeasure <- Csnippet("
  cases = rnorm(rho*IO, sqrt( pow(tau*IO,2) + rho*IO ) );
  if (cases > 0.0) {
    cases = nearbyint(cases);
  } else {
    cases = 0.0;
  }
")


## ----initializer--------------------------------------------------------------
polio_rinit <- Csnippet("
  SB1 = SB1_0;
  SB2 = SB2_0;
  SB3 = SB3_0;
  SB4 = SB4_0;
  SB5 = SB5_0;
  SB6 = SB6_0;
  IB = 0;
  IO = IO_0 * P;
  SO = SO_0 * P;
")


## ----trans--------------------------------------------------------------------
polio_partrans <- parameter_trans(
  log=c("psi","rho","tau","sigma_dem","sigma_env"),
  logit=c("SO_0","IO_0")
)


## ----pomp---------------------------------------------------------------------
polio <- pomp(
  data=subset(polio_data, 
    (time > polio_t0 + 0.01) & (time < 1953+1/12+0.01),
    select=c("cases","time")),
  times="time",
  t0=polio_t0,
  params=polio_params_guess,
  rprocess = euler(step.fun = polio_rprocess, delta.t=1/12),
  rmeasure= polio_rmeasure,
  dmeasure = polio_dmeasure,
  covar=polio_covar,
  obsnames = polio_obsnames,
  statenames = polio_statenames,
  paramnames = polio_paramnames,
  rinit=polio_rinit,
  partrans=polio_partrans
)
plot(polio)


## ----run_level----------------------------------------------------------------
run_level=2
polio_Np <-          switch(run_level,100, 1e3, 1e4)
polio_Nmif <-        switch(run_level, 10, 100, 400)
polio_Nreps_eval <-  switch(run_level,  2,  10,  20)
polio_Nreps_local <- switch(run_level, 10,  20,  40)
polio_Nreps_global <-switch(run_level, 10,  20, 100)
polio_Nsim <-        switch(run_level, 50, 100, 500) 


## ----parallel-setup,cache=FALSE-----------------------------------------------
require(doParallel)
registerDoParallel()


## ----pf1----------------------------------------------------------------------
stew(file=sprintf("pf1-%d.rda",run_level),{
  t1 <- system.time(
    pf1 <- foreach(i=1:20,.packages='pomp') %dopar% pfilter(
      polio,Np=polio_Np)
  )
},seed=493536993,kind="L'Ecuyer")
L1 <- logmeanexp(sapply(pf1,logLik),se=TRUE)


## ----persistence-sim----------------------------------------------------------
stew(sprintf("persistence-%d.rda",run_level),{
  t_sim <- system.time(
    sim <- foreach(i=1:polio_Nsim,
      .packages='pomp') %dopar% simulate(polio)
  )
},seed=493536993,kind="L'Ecuyer")


## ----persistence-analysis-----------------------------------------------------
no_cases_data <- sum(obs(polio)==0)
no_cases_sim <- sum(sapply(sim,obs)==0)/length(sim)
fadeout1_sim <- sum(sapply(sim,function(po)states(po)["IB",]
  +states(po)["IO",]<1))/length(sim)
fadeout100_sim <- sum(sapply(sim,function(po)states(po)["IB",]
  +states(po)["IO",]<100))/length(sim)
imports_sim <- coef(polio)["psi"]*mean(sapply(sim,
  function(po) mean(states(po)["SO",]+states(po)["SB1",]
    +states(po)["SB2",]+states(po)["SB3",]+states(po)["SB4",]
    +states(po)["SB5",]+states(po)["SB6",])))/12


## ----plot_simulated-----------------------------------------------------------
mle_simulation <- simulate(polio,seed=127)
plot(mle_simulation)


## ----rwsd---------------------------------------------------------------------
polio_rw.sd_rp <- 0.02
polio_rw.sd_ivp <- 0.2
polio_cooling.fraction.50 <- 0.5
polio_rw.sd <- rw.sd(
  b1=polio_rw.sd_rp, b2=polio_rw.sd_rp,
  b3=polio_rw.sd_rp, b4=polio_rw.sd_rp,
  b5=polio_rw.sd_rp, b6=polio_rw.sd_rp,
  psi=polio_rw.sd_rp, rho=polio_rw.sd_rp,
  tau=polio_rw.sd_rp, sigma_dem=polio_rw.sd_rp,
  sigma_env=polio_rw.sd_rp,
  IO_0=ivp(polio_rw.sd_ivp), SO_0=ivp(polio_rw.sd_ivp)
)


## ----mif----------------------------------------------------------------------
stew(sprintf("mif-%d.rda",run_level),{
  t2 <- system.time({
    m2 <- foreach(i=1:polio_Nreps_local,
      .packages='pomp', .combine=c) %dopar% mif2(polio,
           Np=polio_Np,
           Nmif=polio_Nmif,
           cooling.fraction.50=polio_cooling.fraction.50,
           rw.sd=polio_rw.sd)	
    lik_m2 <- foreach(i=1:polio_Nreps_local,
      .packages='pomp', .combine=rbind) %dopar% logmeanexp(
        replicate(polio_Nreps_eval,logLik(
	  pfilter(polio,params=coef(m2[[i]]),Np=polio_Np))),
        se=TRUE)
  })
},seed=318817883,kind="L'Ecuyer")


## ----search-save--------------------------------------------------------------
r2 <- data.frame(logLik=lik_m2[,1],logLik_se=lik_m2[,2],
  t(sapply(m2,coef)))
if (run_level>1) 
  write.table(r2,file="polio_params.csv",append=TRUE,
  col.names=FALSE,row.names=FALSE)
summary(r2$logLik,digits=5)


## ----pairs--------------------------------------------------------------------
pairs(~logLik+psi+rho+tau+sigma_dem+sigma_env,
  data=subset(r2,logLik>max(logLik)-20))


## ----box----------------------------------------------------------------------
polio_box <- rbind(
  b1=c(-2,8), b2=c(-2,8),
  b3=c(-2,8), b4=c(-2,8),
  b5=c(-2,8), b6=c(-2,8),
  psi=c(0,0.1), rho=c(0,0.1), tau=c(0,0.1),
  sigma_dem=c(0,0.5), sigma_env=c(0,1),
  SO_0=c(0,1), IO_0=c(0,0.01)
)


## ----box_eval-----------------------------------------------------------------
stew(file=sprintf("box_eval-%d.rda",run_level),{
  t3 <- system.time({
    m3 <- foreach(i=1:polio_Nreps_global,.packages='pomp',
      .combine=c) %dopar% mif2(
        m2[[1]],
        params=c(apply(polio_box,1,function(x)runif(1,x[1],x[2])),
	  polio_fixed_params)
      )
    lik_m3 <- foreach(i=1:polio_Nreps_global,.packages='pomp',
      .combine=rbind) %dopar% logmeanexp(
        replicate(polio_Nreps_eval,
        logLik(pfilter(polio,
	  params=coef(m3[[i]]),Np=polio_Np))), 
        se=TRUE)
  })
},seed=290860873,kind="L'Ecuyer")


## ----global-results-----------------------------------------------------------
r3 <- data.frame(logLik=lik_m3[,1],logLik_se=lik_m3[,2],t(sapply(m3,coef)))
if(run_level>1) write.table(r3,file="polio_params.csv",append=TRUE,col.names=FALSE,row.names=FALSE)
summary(r3$logLik,digits=5)


## ----pairs_global-------------------------------------------------------------
pairs(~logLik+psi+rho+tau+sigma_dem+sigma_env,
  data=subset(r3,logLik>max(logLik)-20))


## ----nbinom-------------------------------------------------------------------
nb_lik <- function(theta) -sum(dnbinom(as.vector(obs(polio)),
  size=exp(theta[1]),prob=exp(theta[2]),log=TRUE))
nb_mle <- optim(c(0,-5),nb_lik)
-nb_mle$value


## ----arma---------------------------------------------------------------------
log_y <- log(as.vector(obs(polio))+1)
arma_fit <- arima(log_y,order=c(2,0,2),seasonal=list(order=c(1,0,1),period=12))
arma_fit$loglik-sum(log_y)


## ----param_file,out.width="10cm"----------------------------------------------
polio_params <- read.table("polio_params.csv",row.names=NULL,
  header=TRUE)
pairs(~logLik+psi+rho+tau+sigma_dem+sigma_env,
  data=subset(polio_params,logLik>max(logLik)-20))


## ----global_rho_code,eval=F,echo=T--------------------------------------------
## plot(logLik~rho,data=subset(r3,logLik>max(r3$logLik)-10),log="x")


## ----global_rho_plot,eval=T,echo=F,out.width="7cm",fig.width=5,fig.height=3----
par(mai=c(0.8,0.8,0.1,0.1))
plot(logLik~rho,data=subset(r3,logLik>max(r3$logLik)-10),log="x")


## ----mif_diagnostics_code,echo=T,eval=F---------------------------------------
## plot(m3[r3$logLik>max(r3$logLik)-10])


## ----mif_diagnostics_plot,out.width="7.5cm",echo=F----------------------------
plot(m3[r3$logLik>max(r3$logLik)-10])


## ----likelihood_convergence,out.width="10cm"----------------------------------
loglik_convergence <- do.call(cbind,
  traces(m3[r3$logLik>max(r3$logLik)-10],"loglik"))
matplot(loglik_convergence,type="l",lty=1,
  ylim=max(loglik_convergence,na.rm=T)+c(-10,0))

