options(
  keep.source=TRUE,
  encoding="UTF-8"
)

set.seed(594709947L)
library(ggplot2)
library(dplyr)
library(plyr)
library(reshape2)
library(foreach)
library(pomp)
library(doParallel)

cores = 5  # The number of cores on this machine 
cl = makeCluster(cores)
registerDoParallel(cl)
stopifnot(packageVersion("pomp")>="2.0")
theme_set(theme_bw())


pertussis = read.csv("pertussis.csv", header = TRUE) 

pertussis = pertussis %>% 
  filter(state == "MI") %>% 
  filter(week > 200000) %>%
  mutate(time = 2000 + as.integer((week-200000)/100) + ((week-200000)-1)%%100/52,
         year = 2000 + as.integer((week-200000)/100),
         week = ((week-200000))%%100) %>%
  select(time, year, week, cases, incidence_per_capita)
pertussis = pertussis %>% filter(year >= 2006)

demographics = data.frame(year = 2000:2011,
                          N = c(9918000, 10010000, 10040000, 10070000, 10090000, # population
                                10090000, 10080000, 10050000, 9999000, 9955000,
                                9931000, 9882000),
                          mu_D = c(0.0087, 0.0086, 0.0087, 0.0086, 0.0085, # death rate
                                   0.0086, 0.0086, 0.0087, 0.0089, 0.0087,
                                   0.0089, 0.0091),
                          mu_B = c(0.0137, 0.0133, 0.0129, 0.013, 0.0128, # birth rate
                                   0.0126, 0.0126, 0.0124, 0.0121, 0.0118, 
                                   0.0116, 0.0116),
                          v = 0.832)
pertussis = pertussis %>% left_join(demographics, by = "year")

pertussis_covar = covariate_table(time = pertussis$time,
                                  N = pertussis$N,
                                  mu_D = pertussis$mu_D,
                                  mu_B = pertussis$mu_B,
                                  v = pertussis$v,
                                  times = "time")



pertussis_rprocess = "
double dN_VSv = rbinom(V, 1-exp(-mu_VSv*dt));

double dN_SvI = rbinom(Sv, 1-exp(-Beta*I*dt));
double dN_SnvI = rbinom(Snv, 1-exp(-Beta*I*dt));
double dN_IR = rbinom(I, 1-exp(-mu_IR*dt));

double dN_VD = rbinom(V, 1-exp(-mu_D*dt));
double dN_SvD = rbinom(Sv, 1-exp(-mu_D*dt));
double dN_SnvD = rbinom(Snv, 1-exp(-mu_D*dt));
double dN_ID = rbinom(I, 1-exp(-mu_D*dt));
double dN_RD = rbinom(R, 1-exp(-mu_D*dt));

int B_Snv = nearbyint(N * mu_B * (1-v) * dt);
int B_V = nearbyint(N * mu_B * v * dt);

V += B_V - dN_VSv - dN_VD;
Sv += dN_VSv - dN_SvI - dN_SvD;
if(Sv<0) Sv = 0;

Snv += B_Snv - dN_SnvI - dN_SnvD;
if(Snv<0) Snv = 0;
I += dN_SvI + dN_SnvI- dN_IR - dN_ID;
if(I<0) I = 0;
R += dN_IR - dN_RD;
if(R<0) R = 0;
H += dN_IR;
if(H<0) H = 0;
"

pertussis_dmeasure = "
lik = dpois(cases,rho*H+1e-10,give_log);
"
pertussis_rmeasure = "
cases = rpois(rho*H+1e-10);
"
# pertussis_dmeasure = "
# lik =  dbinom(cases, H, rho+1e-10, give_log);
# "
# pertussis_rmeasure = "
# cases = rbinom(H, rho+1e-10);
# "
pertussis_rinit = "
V=6000000;
Sv=1000000;
Snv=3000000;
I=10;
R=0;
H=0;
"
pertussis_statenames = c("V","Snv","Sv","I","R","H")
pertussis_paramnames = c("Beta","mu_VSv","mu_IR","rho")

pertussis_obsnames = "cases"
pertussis_params_guess = c(Beta = 0.05, mu_VSv=0.01, 
                           mu_IR = 0.01, rho = 0.1)


pertussis.pomp = pomp(data = pertussis %>% select(cases, time),
                      times = "time", 
                      t0 = pertussis$time[1]-1/52,
                      params = pertussis_params_guess,
                      rprocess = euler(step.fun=Csnippet(pertussis_rprocess),
                                       delta.t=1/365),
                      rmeasure = Csnippet(pertussis_rmeasure),
                      dmeasure = Csnippet(pertussis_dmeasure),
                      covar = pertussis_covar,
                      obsnames = pertussis_obsnames,
                      statenames = pertussis_statenames,
                      paramnames=pertussis_paramnames,
                      rinit=Csnippet(pertussis_rinit),
                      partrans=parameter_trans(log=c("Beta","mu_IR","mu_VSv"),
                                               logit="rho"),
                      accumvars="H")


run_level = 1
switch(run_level, {
  pertussis_Np=100; pertussis_Nmif=10; pertussis_Neval=10;
  pertussis_Nglobal=10; pertussis_Nlocal=10
},{
  pertussis_Np=20000; pertussis_Nmif=100; pertussis_Neval=10;
  pertussis_Nglobal=10; pertussis_Nlocal=10
}
)

gs_box2 = rbind(Beta = c(0.1, 5),
                mu_VSv = c(0.01, 0.5),
                mu_IR =c(0.00001, 0.01),
                rho = c(0.4, 0.99))

pertussis_rw.sd = 0.02
pertussis_cooling.fraction.50 = 0.5                

stew(file=sprintf("gs-sirv-%d.rda",run_level),{
  t_local = system.time({
    mifs_global = foreach(i=1:pertussis_Nglobal,
                          .packages='pomp', .combine=c,
                          .export = c("pertussis.pomp", "pertussis_Np", 
                                      "pertussis_Nmif", "pertussis_rw.sd",
                                      "pertussis_cooling.fraction.50", 
                                      "gs_box2")) %dopar%  {
                                        mif2(pertussis.pomp,
                                             Np=pertussis_Np,
                                             Nmif=pertussis_Nmif,
                                             cooling.fraction.50=pertussis_cooling.fraction.50,
                                             rw.sd=rw.sd(
                                               Beta=pertussis_rw.sd,
                                               mu_VSv=pertussis_rw.sd,
                                               mu_IR=pertussis_rw.sd,
                                               rho=pertussis_rw.sd),
                                             params = apply(gs_box2, 1, function(x)runif(1,x[1],x[2]))
                                        )
                                      }
  })
},seed=12345,kind="L'Ecuyer")


## ----lik_local_eval,cache=FALSE------------------------------------------
stew(file=sprintf("lik_gs-sirv-%d.rda",run_level),{
  t_global_eval <- system.time({
    liks_global <- foreach(i=1:pertussis_Nglobal,
                           .combine=rbind,.packages='pomp',
                           .export = c("pertussis_Neval", "pertussis.pomp",
                                       "mifs_global", "pertussis_Np"))%dopar% {
                                         evals <- replicate(pertussis_Neval, logLik(
                                           pfilter(pertussis.pomp, params=coef(mifs_global[[i]]), Np=pertussis_Np)))
                                         logmeanexp(evals, se=TRUE)
                                       }
  })
},seed=12345,kind="L'Ecuyer")

results_global <- data.frame(logLik=liks_global[,1],
                             logLik_se=liks_global[,2],t(sapply(mifs_global,coef)))


## ----lik_local_summary---------------------------------------------------
summary(results_global$logLik,digits=5)


## ----pairs_local_code,eval=FALSE,echo=T----------------------------------
## pairs(~logLik+Beta+mu_IR+rho,
##   data=subset(results_local,logLik>max(logLik)-50))


## ----pairs_local_plot,eval=TRUE,echo=FALSE,out.width="12cm"--------------
pairs(~logLik+Beta+mu_IR+rho+mu_VSv,
      data=subset(results_global))


plot(mifs_global)
