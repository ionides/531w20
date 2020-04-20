---
title: "Solution to Homework 8"
author: "STATS 531, Winter 2020"
output: html_document
---

First, we set up a parallel computing environment which uses multiple cores on either a laptop or a slurm cluster such as Great Lakes.


```r
library(doParallel)
```

```
## Loading required package: foreach
```

```
## Loading required package: iterators
```

```
## Loading required package: parallel
```

```r
cluster_cores <- as.numeric(Sys.getenv('SLURM_NTASKS_PER_NODE'))
CLUSTER <- !is.na(cluster_cores)
if(CLUSTER) registerDoParallel(cluster_cores) else registerDoParallel()

library(doRNG)
```

```
## Loading required package: rngtools
```

```r
registerDoRNG(1929384)
```

We now loading the data and reproduce the results in the notes12.


```r
options(
  keep.source=TRUE,
  stringsAsFactors=FALSE,
  encoding="UTF-8"
  )
library(ggplot2)
theme_set(theme_bw())
library(plyr)
library(reshape2)
library(magrittr)
library(pomp)
# bsflu_data <- read.table("http://ionides.github.io/531w20/12/bsflu_data.txt")
# read from a local copy when the internet is unreliable!
bsflu_data <- read.table("bsflu_data.txt")
bsflu_statenames <- c("S","I","R1","R2")
bsflu_paramnames <- c("Beta","mu_IR","rho","mu_R1","mu_R2")
bsflu_obsnames <- "B"

bsflu_dmeasure <- "
  lik = dpois(B,rho*R1+1e-6,give_log);
"

bsflu_rmeasure <- "
  B = rpois(rho*R1+1e-6);
"

bsflu_rprocess <- "
  double dN_SI = rbinom(S,1-exp(-Beta*I*dt));
  double dN_IR1 = rbinom(I,1-exp(-dt*mu_IR));
  double dN_R1R2 = rbinom(R1,1-exp(-dt*mu_R1));
  double dN_R2R3 = rbinom(R2,1-exp(-dt*mu_R2));
  S -= dN_SI;
  I += dN_SI - dN_IR1;
  R1 += dN_IR1 - dN_R1R2;
  R2 += dN_R1R2 - dN_R2R3;
"
bsflu_rinit <- "
 S=762;
 I=1;
 R1=0;
 R2=0;
"

bsflu2 <- pomp(
  data=subset(bsflu_data,select=c(day,B)),
  times="day",
  t0=0,
  rprocess=euler(
    step.fun=Csnippet(bsflu_rprocess),
    delta.t=1/12
  ),
  rmeasure=Csnippet(bsflu_rmeasure),
  dmeasure=Csnippet(bsflu_dmeasure),
  partrans=parameter_trans(
    log=c("Beta","mu_IR","mu_R1","mu_R2"),
    logit="rho"),
  obsnames = bsflu_obsnames,
  statenames=bsflu_statenames,
  paramnames=bsflu_paramnames,
  rinit=Csnippet(bsflu_rinit)
)

run_level <- 3
bsflu_Np <-      switch(run_level, 100, 20000, 60000)
bsflu_Nmif <-    switch(run_level,  10,   100,   300)
bsflu_Neval <-   switch(run_level,  10,    10,    10)
bsflu_Nglobal <- switch(run_level,  10,    10,   100)
bsflu_Nlocal <-  switch(run_level,  10,    10,    20)

# bsflu_params <- data.matrix(read.table("http://ionides.github.io/531w20/12/mif_bsflu_params.csv",row.names=NULL,header=TRUE))
# read from a local copy when the internet is unreliable!
bsflu_params <- data.matrix(read.table("mif_bsflu_params.csv",row.names=NULL,header=TRUE))
bsflu_mle <- bsflu_params[which.max(bsflu_params[,"logLik"]),][bsflu_paramnames]
bsflu_fixed_params <- c(mu_R1=1/(sum(bsflu_data$B)/512),mu_R2=1/(sum(bsflu_data$C)/512))


bsflu_box <- rbind(
  Beta=c(0.001,0.01),
  mu_IR=c(0.5,2),
  rho = c(0.5,1)
)

bsflu_rw.sd <- 0.02
bsflu_cooling.fraction.50 <- 0.5
stew(file=sprintf("box_eval-%d.rda",run_level),{
  t_global <- system.time({
    mifs_global <- foreach(i=1:bsflu_Nglobal,
      .packages='pomp', .combine=c) %dopar% mif2(     
        bsflu2,
        params=c(apply(bsflu_box,1,function(x)runif(1,x[1],x[2])),
	  bsflu_fixed_params),
        Np=bsflu_Np,
        Nmif=bsflu_Nmif,
        cooling.fraction.50=bsflu_cooling.fraction.50,
        rw.sd=rw.sd(
          Beta=bsflu_rw.sd,
          mu_IR=bsflu_rw.sd,
          rho=bsflu_rw.sd
        )
      )
    })
})

stew(file=sprintf("lik_global_eval-%d.rda",run_level),{
  t_global_eval <- system.time({
    liks_global <- foreach(i=1:bsflu_Nglobal,
      .packages='pomp',.combine=rbind) %dopar% {
        evals <- replicate(bsflu_Neval,
	  logLik(pfilter(bsflu2,params=coef(mifs_global[[i]]),Np=bsflu_Np)))
          logmeanexp(evals, se=TRUE)
    }
  })
})

results_global <- data.frame(logLik=liks_global[,1],
  logLik_se=liks_global[,2],t(sapply(mifs_global,coef)))
summary(results_global$logLik,digits=5)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  -75.77  -74.59  -74.25  -74.26  -73.97  -73.00
```

```r
t_global
```

```
##       user     system    elapsed 
## 118781.432    326.499   5109.939
```

```r
plot(mifs_global)
```

![plot of chunk load_notes12](figure/load_notes12-1.png)![plot of chunk load_notes12](figure/load_notes12-2.png)

----------

**<big>Question 8.1</big>**. Assessing and improving algorithmic parameters.

From the diagnostic plots above, at run level 3, we can see the effective sample sizes are basically large (except the time point with less than 500 effective sample size) and all the parameters are converged after 200 MIF iteration.
Since this run level takes considerable computational resources, we investigate whether we can get comparable results for half the effort.

In the following code, we test a new set of algorithmic parameters.


```r
run_level <- 4
bsflu_Np <-      switch(run_level, 100, 20000, 60000, 50000)
bsflu_Nmif <-    switch(run_level,  10,   100,   300,   250)
bsflu_Neval <-   switch(run_level,  10,    10,    10,    10)
bsflu_Nglobal <- switch(run_level,  10,    10,   100,    75)
bsflu_Nlocal <-  switch(run_level,  10,    10,    20,    20)
bsflu_cooling.fraction.50 <- 0.6
```

We will change the number of particles (`Np`) to 5 &times; 10<sup>4</sup> the number of the MIF iteration (`Nmif`) to 250, and the number of global searches to 75.
Moreover, we will change `bsflu_cooling.fraction.50` to 0.6.
Our hypothesis is that these changes can reduce the computation time and give us satisfactory convergence with acceptable effective sample sizes.



```r
stew(file=sprintf("Mif-8.1-%d.rda",run_level),{
  
  t_global.2 <- system.time({
    mifs_global.2 <- foreach(i=1:bsflu_Nglobal,
      .packages='pomp', .combine=c) %dopar% 
      mif2(
        mifs_global[[1]],
        start=c(apply(bsflu_box,1,function(x)runif(1,x[1],x[2])),bsflu_fixed_params),
        Np=bsflu_Np,
        Nmif=bsflu_Nmif,
        cooling.fraction.50=bsflu_cooling.fraction.50        
      )
  })
})
##----eval llh--------

stew(file=sprintf("lik-8.1-%d.rda",run_level),{
  t_global_eval.2 <- system.time({
    liks_global.2 <- foreach(i=1:bsflu_Nglobal,
      .packages='pomp',.combine=rbind) %dopar% {
      evals <- replicate(bsflu_Neval,
        logLik(pfilter(bsflu2,params=coef(mifs_global.2[[i]]),Np=bsflu_Np)))
      logmeanexp(evals, se=TRUE)
    }
  })
})
t_global.2
```

```
##      user    system   elapsed 
## 61596.856   291.778  2620.883
```

```r
results_global.2 <- data.frame(logLik=liks_global.2[,1],
  logLik_se=liks_global.2[,2],t(sapply(mifs_global.2,coef)))
summary(results_global.2$logLik,digits=5)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  -75.33  -74.48  -74.24  -74.12  -73.73  -72.53
```

```r
plot(mifs_global.2)
```

![plot of chunk 8.1](figure/8.1-1.png)![plot of chunk 8.1](figure/8.1-2.png)

* We can see, this change reduce the time since less computations are needed. However, it found a better maximized likelihood
-72.53.
Moreover, comparing the median and the mean of the likelihood, these modifications generally also give us better MLEs.

* The effective sample sizes are large and basically the same as before, so it is still acceptable. 

* $\beta$ and $\rho$ are converged in most of the particles. There are only few particles in which $\mu_{IR}$ is not converged.

-----------

**<big>Question 8.2</big>**.  Finding sharp peaks in the likelihood surface.

The drawing method with scale invariant properties is:
First,taking logarithm of the orginal domain. Secondly, drawing sample uniformly from it. Thirdly, taking exponetiation of the sample.

In `R`, we can implement with command `exp(runif(1,log(a),log(b)))`, where (a,b) is the orginal domain.

Following is the code the implements this drawing method (only for $\beta$).


```r
stew(file=sprintf("box_eval_log new-%d.rda",run_level),{
    t_global.3 <- system.time({
    mifs_global.3 <- foreach(i=1:bsflu_Nglobal,
      .packages='pomp', .combine=c) %dopar% 
      mif2(
        mifs_global[[1]],
        params=c(Beta=exp(runif(1,log(bsflu_box[1,1]),log(bsflu_box[1,2])))
           ,apply(bsflu_box[2:3,],1,function(x)
	     runif(1,x[1],x[2])),bsflu_fixed_params),
        Np=bsflu_Np,
        Nmif=bsflu_Nmif
      )
  })
})

stew(file=sprintf("lik_global_eval_log new-%d.rda",run_level),{
  t_global_eval.3 <- system.time({
    liks_global.3 <- foreach(i=1:bsflu_Nglobal,
      .packages='pomp',.combine=rbind) %dopar% {
      evals <- replicate(bsflu_Neval, logLik(
        pfilter(bsflu2,params=coef(mifs_global.3[[i]]),Np=bsflu_Np)))
      logmeanexp(evals, se=TRUE)
    }
  })
})

results_global.3 <- data.frame(logLik=liks_global.3[,1],
  logLik_se=liks_global.3[,2],t(sapply(mifs_global.3,coef)))
summary(results_global.3$logLik,digits=5)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  -75.50  -74.75  -74.44  -74.37  -74.10  -72.33
```

Comparing this result with the ones above, we can see this method actually does not help much in finding sharp peaks, although it has scale invariant property.

To illustrate this, we can see the following histograms of 5000 samples with this two drawing methods.


```r
par(mfrow=c(1,2))
hist(runif(5000,0.001,0.1),ylim=c(0,2600),xlab=expression(beta),
     main="Uniformly drawing \n from orginal scale")
hist(exp(runif(5000,log(0.001),log(0.1))),xlab=expression(beta),
     ylim=c(0,2600),
     main="Uniformly drawing \n from log scale")
```

![plot of chunk hist](figure/hist-1.png)

As we can see the drawing around the peak ($\beta \approx 0.004$) does not increase.

----------------

**<big>Question 8.3</big>**.  Construct a profile likelihood.

First, we construct a parameter box for computing the profile likelihood of $\beta$. 


```r
profile_run_level <- 3
bsflu_profile_pts <- switch(profile_run_level, 5,10,30)
bsflu_profile_reps <- switch(profile_run_level, 3,10,20)
bsflu_profile_Np_mif <- switch(profile_run_level, 100,1000,5000)
bsflu_profile_Nmif <- switch(profile_run_level, 2,20,100)
bsflu_profile_Neval <- switch(profile_run_level, 10,10,10)
bsflu_profile_Np_eval <- switch(profile_run_level, 100,1000,10000)

profile.starts <- profileDesign(  
  Beta=exp(seq(log(0.001),log(0.01),length.out=bsflu_profile_pts)),
  lower=c(mu_IR=0.5,rho=0.5),
  upper=c(mu_IR=2,rho=1),
  nprof=bsflu_profile_reps
)
```

Then, from each start point, we use `mif2` to find the maximized likelihood and the correponding MLE. Since we need to find the profile likelihood of $\beta$, we need to fix it during the iterated filtering. Therefore, only $\mu_{IR}$ and $\rho$ have random walk stand deviation (`rw.sd`).


```r
stew(file=sprintf("profile_beta-%d.rda",profile_run_level),{ 
  t_global.4 <- system.time({
      prof.llh<- foreach(start=iter(profile.starts,"row"),
        .combine=rbind) %dopar%{
        # Find MLE
        mif2(
          mifs_global[[1]],
          params=c(start,bsflu_fixed_params),
          Np=bsflu_profile_Np_mif,
	  Nmif=bsflu_profile_Nmif,
          rw.sd=rw.sd(
            mu_IR=bsflu_rw.sd,
            rho=bsflu_rw.sd
          )
        )->mifs_global.4
        # evaluate llh
        evals = replicate(bsflu_profile_Neval,logLik(pfilter(mifs_global.4,
	  Np=bsflu_profile_Np_eval)))
        ll=logmeanexp(evals, se=TRUE)        
        
        data.frame(as.list(coef(mifs_global.4)),
                   loglik = ll[1],
                   loglik.se = ll[2])
      }
  })
})
```

At each value of $\beta$, we pick the MLEs which gives us the maximized likelihood evaluation, and we do the iterated filtering again.

This step is optional. It may be a good idea to carry out these additional local seaches to refine the results of the global search, or it may be unnecessary. 


```r
## filiter again at the maxima

prof.llh %>% 
  ddply(~Beta,subset,rank(-loglik)<=10) %>%
  subset(select=bsflu_paramnames) -> pars


## mif2 again
stew(file=sprintf("profile_beta-2-%d.rda",profile_run_level),{  
  t_global.5 <- system.time({
    prof.llh<- foreach(i=1:(nrow(pars)),
      .packages='pomp', .combine=rbind) %dopar%{
      # Find MLE
      mif2(
        mifs_global[[1]],
        start=unlist(pars[i,]),
        Np=bsflu_profile_Np_mif,
	Nmif=bsflu_profile_Nmif,
        rw.sd=rw.sd(
          mu_IR=bsflu_rw.sd,
          rho=bsflu_rw.sd
        )
      )->mifs_global.5
      # evaluate llh 
      pf= replicate(bsflu_profile_Neval,
        pfilter(mifs_global.5,Np=bsflu_profile_Np_eval))
      evals=sapply(pf,logLik)
      ll=logmeanexp(evals, se=TRUE)  
      nfail=sapply(pf,getElement,"nfail")
      
      data.frame(as.list(coef(mifs_global.5)),
                 loglik = ll[1],
                 loglik.se = ll[2],
                 nfail.max=max(nfail))
    }
  })
})
```

Finally, for each value of $\beta$, we use the MLE of the largest likelihood to construct the profile likelihood.
The `loess` method is used to fit an estimate of the profile through these points.
An approximate 95% confidence interval of $\beta$ is $\{\beta: \max[\ell^{profile}(\beta)]-\ell^{profile}(\beta)<1.92\}$



```r
prof.llh %<>%
  subset(nfail.max==0) %>%
  mutate(Beta=exp(signif(log(Beta),5))) %>%
  ddply(~Beta,subset,rank(-loglik)<=1)

aa=max(prof.llh$loglik)
bb=aa-1.92
CI=which(prof.llh$loglik>=bb)
cc=prof.llh$Beta[min(CI)]
dd=prof.llh$Beta[max(CI)]


prof.llh %>%
  ggplot(aes(x=Beta,y=loglik))+
  geom_point()+
  geom_smooth(method="loess")+
  geom_hline(aes(yintercept=aa),linetype="dashed")+
  geom_hline(aes(yintercept=bb),linetype="dashed")+
  geom_vline(aes(xintercept=cc),linetype="dashed")+
  geom_vline(aes(xintercept=dd),linetype="dashed")
```

```
## `geom_smooth()` using formula 'y ~ x'
```

![plot of chunk plot_profile](figure/plot_profile-1.png)

```r
c(lower=cc,upper=dd)
```

```
##       lower       upper 
## 0.004189391 0.004189391
```

As we can see, the 95% confidence interval of $\beta$ is 
$[0.00419, 0.00419]$, 
which shows that $\beta$ is statistically significant.

----------------------

**<big>Question 8.4</big>**.

All responses were given full credit as long as they were consistent with the solutions presented.


------------

### Acknowledgements

Parts of this solution are adapted from a previous homework submission by Rui Zhang.

---------------------
 
