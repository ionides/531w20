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
library(plyr)
library(reshape2)
library(foreach)
#library(doMC)
library(pomp)
stopifnot(packageVersion("pomp")>="2.0")




## ----sir-sim1,out.width="10cm"------------------------------------------------
sims <- simulate(sir,params=c(Beta=1.8,mu_IR=1,rho=0.9,N=2600),
  nsim=20,format="data.frame",include=TRUE)
ggplot(sims,mapping=aes(x=day,y=B,group=.id,color=.id=="data"))+
  geom_line()+guides(color=FALSE)


## ----sir-pfilter-1,results='markup',cache=T-----------------------------------
pf <- pfilter(sir,Np=5000,params=c(Beta=2,mu_IR=1,rho=0.8,N=2600))
logLik(pf)


## ----sir-pfilter-2,results='markup',cache=T-----------------------------------
pf <- replicate(10,
  pfilter(sir,Np=5000,params=c(Beta=2,mu_IR=1,rho=0.8,N=2600))
)
print(ll <- sapply(pf,logLik))


## ----logmeanexp---------------------------------------------------------------
logmeanexp(ll,se=TRUE)


## ----parallel-setup,cache=FALSE-----------------------------------------------
library(doParallel)
registerDoParallel()
library(doRNG)
registerDoRNG(3899882)


## ----sir-like-slice,cache=TRUE,results='hide'---------------------------------
p <- sliceDesign(
  c(Beta=2,mu_IR=1,rho=0.8,N=2600),
  Beta=rep(seq(from=0.5,to=4,length=40),each=3),
  mu_IR=rep(seq(from=0.5,to=2,length=40),each=3)) 

foreach (theta=iter(p,"row"),
  .combine=rbind,.inorder=FALSE) %dopar% {
    pfilter(sir,params=unlist(theta),Np=5000) -> pf
    theta$loglik <- logLik(pf)
    theta
  } -> p


## ----sir-like-slice-plot,cache=TRUE,results="hide",echo=F,out.width="8cm"-----
foreach (v=c("Beta","mu_IR")) %do% 
{
  x <- subset(p,slice==v)
  plot(x[[v]],x$loglik,xlab=v,ylab="loglik")
}



## ----sir-grid1,cache=TRUE-----------------------------------------------------
expand.grid(Beta=seq(from=1,to=4,length=50),
            mu_IR=seq(from=0.7,to=3,length=50),
            rho=0.8,
            N=2600) -> p

foreach (theta=iter(p,"row"),.combine=rbind,
  .inorder=FALSE) %dopar% {
     pfilter(sir,params=unlist(theta),Np=5000) -> pf
     theta$loglik <- logLik(pf)
     theta
  } -> p


## ----sir-grid1-plot,cache=TRUE,out.width="11cm"-------------------------------
pp <- mutate(p,loglik=ifelse(loglik>max(loglik)-100,loglik,NA))
ggplot(data=pp,mapping=aes(x=Beta,y=mu_IR,z=loglik,fill=loglik))+
  geom_tile(color=NA)+
  scale_fill_gradient()+
  geom_contour(color='black',binwidth=3)+
  labs(x=expression(beta),y=expression(mu[IR]))


## ----bsflu_rprocess-----------------------------------------------------------
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


## ----bsflu_measure------------------------------------------------------------
bsflu_dmeasure <- "
  lik = dpois(B,rho*R1+1e-10,give_log);
"

bsflu_rmeasure <- "
  B = rpois(rho*R1+1e-10);
"


## ----bsflu_rinit--------------------------------------------------------------
bsflu_rinit <- "
 S=762;
 I=1;
 R1=0;
 R2=0;
"


## ----Csnippets_bsflu1---------------------------------------------------------
bsflu_statenames <- c("S","I","R1","R2")
bsflu_paramnames <- c("Beta","mu_IR","rho","mu_R1","mu_R2")


## ----Csnippets_bsflu2---------------------------------------------------------
bsflu_data <- read.table("bsflu_data.txt")

bsflu2 <- pomp(
  data=subset(bsflu_data,select=c(day,B)),
  times="day",
  t0=0,
  rprocess=euler(
    step.fun=Csnippet(bsflu_rprocess),
    delta.t=1/12),
  rmeasure=Csnippet(bsflu_rmeasure),
  dmeasure=Csnippet(bsflu_dmeasure),
  partrans=parameter_trans(
    log=c("Beta","mu_IR","mu_R1","mu_R2"),
    logit="rho"),
  statenames=bsflu_statenames,
  paramnames=bsflu_paramnames,
  rinit=Csnippet(bsflu_rinit)
)


## ----run_level----------------------------------------------------------------
run_level <- 2
switch(run_level, {
  bsflu_Np=100; bsflu_Nmif=10; bsflu_Neval=10;
  bsflu_Nglobal=10; bsflu_Nlocal=10
  },{
  bsflu_Np=20000; bsflu_Nmif=100; bsflu_Neval=10;
  bsflu_Nglobal=10; bsflu_Nlocal=10
  },{
  bsflu_Np=60000; bsflu_Nmif=300; bsflu_Neval=10;
  bsflu_Nglobal=100; bsflu_Nlocal=20}
)


## ----bsflu_params-------------------------------------------------------------
bsflu_params <- data.matrix(
  read.table("mif_bsflu_params.csv",
  row.names=NULL,header=TRUE))
which_mle <- which.max(bsflu_params[,"logLik"])
bsflu_mle <- bsflu_params[which_mle,][bsflu_paramnames]


## ----fixed_params-------------------------------------------------------------
bsflu_fixed_params <- c(mu_R1=1/(sum(bsflu_data$B)/512),
  mu_R2=1/(sum(bsflu_data$C)/512) )


## ----pf,cache=FALSE-----------------------------------------------------------
stew(file=sprintf("pf-%d.rda",run_level),{
  t_pf <- system.time(
    pf <- foreach(i=1:20,.packages='pomp') %dopar% try(
                    pfilter(bsflu2,params=bsflu_mle,Np=bsflu_Np)
                  )
  )
  
},seed=1320290398,kind="L'Ecuyer")

(L_pf <- logmeanexp(sapply(pf,logLik),se=TRUE))


## ----set_cache, eval=FALSE----------------------------------------------------
## opts_chunk$set(
##   cache=TRUE,
##   )


## ----box_search_local,cache=FALSE---------------------------------------------
bsflu_rw.sd <- 0.02; bsflu_cooling.fraction.50 <- 0.5
stew(file=sprintf("local_search-%d.rda",run_level),{
  t_local <- system.time({
  mifs_local <- foreach(i=1:bsflu_Nlocal,
    .packages='pomp', .combine=c) %dopar%  {
      mif2(bsflu2,
        params=bsflu_mle,
        Np=bsflu_Np,
        Nmif=bsflu_Nmif,
        cooling.fraction.50=bsflu_cooling.fraction.50,
        rw.sd=rw.sd(
          Beta=bsflu_rw.sd,
          mu_IR=bsflu_rw.sd,
          rho=bsflu_rw.sd)
      )
    }
  })
},seed=900242057,kind="L'Ecuyer")


## ----lik_local_eval,cache=FALSE-----------------------------------------------
stew(file=sprintf("lik_local-%d.rda",run_level),{
  t_local_eval <- system.time({
  liks_local <- foreach(i=1:bsflu_Nlocal,.combine=rbind)%dopar% {
    evals <- replicate(bsflu_Neval, logLik(
      pfilter(bsflu2,params=coef(mifs_local[[i]]),Np=bsflu_Np)))
    logmeanexp(evals, se=TRUE)
    }
  })
},seed=900242057,kind="L'Ecuyer")

results_local <- data.frame(logLik=liks_local[,1],
  logLik_se=liks_local[,2],t(sapply(mifs_local,coef)))


## ----lik_local_summary--------------------------------------------------------
summary(results_local$logLik,digits=5)


## ----pairs_local_code,eval=FALSE,echo=T---------------------------------------
## pairs(~logLik+Beta+mu_IR+rho,
##   data=subset(results_local,logLik>max(logLik)-50))


## ----pairs_local_plot,eval=TRUE,echo=FALSE,out.width="12cm"-------------------
pairs(~logLik+Beta+mu_IR+rho,
  data=subset(results_local,logLik>max(logLik)-50))


## ----box----------------------------------------------------------------------
bsflu_box <- rbind(
  Beta=c(0.001,0.01),
  mu_IR=c(0.5,2),
  rho = c(0.5,1)
)


## ----box_eval,cache=FALSE-----------------------------------------------------
stew(file=sprintf("box_eval-%d.rda",run_level),{
  t_global <- system.time({
    mifs_global <- foreach(i=1:bsflu_Nglobal,.combine=c) %dopar% {
      mif2(
        mifs_local[[1]],
        params=c(
          apply(bsflu_box,1,function(x)runif(1,x[1],x[2])),
	  bsflu_fixed_params)
      )}
  })
},seed=1270401374,kind="L'Ecuyer")


## ----lik_global_eval,cache=FALSE----------------------------------------------
stew(file=sprintf("lik_global_eval-%d.rda",run_level),{
  t_global_eval <- system.time({
    liks_global <- foreach(i=1:bsflu_Nglobal,
      .combine=rbind) %dopar% {
        evals <- replicate(bsflu_Neval,
          logLik(pfilter(bsflu2,
	    params=coef(mifs_global[[i]]),Np=bsflu_Np)))
        logmeanexp(evals, se=TRUE)
    }
  })
},seed=442141592,kind="L'Ecuyer")

results_global <- data.frame(
  logLik=liks_global[,1],
  logLik_se=liks_global[,2],t(sapply(mifs_global,coef)))
summary(results_global$logLik,digits=5)


## ----save_params,eval=FALSE---------------------------------------------------
## if (run_level>2)
##   write.table(rbind(results_local,results_global),
##     file="mif_bsflu_params.csv",
##       append=TRUE,col.names=FALSE,row.names=FALSE)


## ----pairs_global_code,echo=TRUE,eval=FALSE-----------------------------------
## pairs(~logLik+Beta+mu_IR+rho,
##   data=subset(results_global,logLik>max(logLik)-250))


## ----pairs_global,echo=FALSE,eval=TRUE,out.width="12cm"-----------------------
pairs(~logLik+Beta+mu_IR+rho,
  data=subset(results_global,logLik>max(logLik)-250))


## ----class_mifs_global--------------------------------------------------------
class(mifs_global)
class(mifs_global[[1]])


## ----mifs_global_plot,out.width="8cm",echo=FALSE------------------------------
plot(mifs_global)

