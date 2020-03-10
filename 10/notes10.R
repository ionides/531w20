## ----setup,echo=F,results=F,cache=F--------------------------------------
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



## ----prelims,cache=F-----------------------------------------------------
set.seed(594709947L)
require(ggplot2)
require(plyr)
require(reshape2)
require(pomp)


## ----load-ricker,cache=FALSE---------------------------------------------
ricker <- ricker()


## ----plot-ricker-code,eval=F,echo=T--------------------------------------
## plot(ricker)


## ----plot-ricker,out.width="12cm",echo=F---------------------------------
par(mai=c(0.8,0.8,0.1,0.1))
plot(ricker)


## ----sim-ricker1---------------------------------------------------------
simulated_ricker <- simulate(ricker)


## ----class_simulated_ricker----------------------------------------------
class(simulated_ricker)


## ----plot_simulated_ricker_code,eval=F,echo=T----------------------------
## plot(simulated_ricker)


## ----plot_simulated_ricker,eval=T,echo=F,out.width="12cm"----------------
par(mai=c(0.8,0.5,0.1,0.1))
plot(simulated_ricker)


## ----ricker_to_data_frame------------------------------------------------
y <- as.data.frame(ricker)
head(y,3)
head(simulate(ricker,format="data.frame"))


## ----vectorize_simulate--------------------------------------------------
x <- simulate(ricker,nsim=10)
class(x)
sapply(x,class)
x <- simulate(ricker,nsim=10,format="data.frame")
head(x,3)


## ----plot_sim_demo_code,echo=T,eval=F------------------------------------
## x <- simulate(ricker,nsim=9,format="data.frame",include.data=TRUE)
## ggplot(data=x,aes(x=time,y=y,group=.id,color=(.id=="data")))+
##   geom_line()+guides(color=FALSE)+
##   facet_wrap(~.id,ncol=2)


## ----plot_sim_demo,echo=F,eval=T,out.width="12cm"------------------------
par(mai=c(0.8,0.5,0.1,0.1))
x <- simulate(ricker,nsim=9,format="data.frame",include.data=TRUE)
ggplot(data=x,aes(x=time,y=y,group=.id,color=(.id=="data")))+
  geom_line()+guides(color=FALSE)+
  facet_wrap(~.id,ncol=2)


## ----traj-ricker---------------------------------------------------------
y <- trajectory(ricker)
dim(y)
dimnames(y)


## ----plot_ricker_traj_code,echo=T,eval=F---------------------------------
## plot(time(ricker),y["N",1,],type="l")


## ----plot_ricker_traj,echo=F,eval=T,out.width="12cm"---------------------
par(mai=c(0.8,0.5,0.1,0.1))
plot(time(ricker),y["N",1,],type="l")


## ----coef-ricker---------------------------------------------------------
coef(ricker)


## ----plot_at_different_parameters,out.width="10cm"-----------------------
theta <- coef(ricker)
theta[c("r","N.0")] <- c(5,3)
y <- trajectory(ricker,params=theta)
plot(time(ricker),y["N",1,],type="l")
x <- simulate(ricker,params=theta)
plot(x,var="y")


## ----change_ricker_coef_code,eval=F,echo=T-------------------------------
## coef(ricker,c("r","N.0","sigma")) <- c(39,0.5,1)
## coef(ricker)
## plot(simulate(ricker),var="y")


## ----change_ricker_coef,eval=T,echo=F,out.width="12cm"-------------------
par(mai=c(0.8,0.5,0.1,0.1))
coef(ricker,c("r","N.0","sigma")) <- c(39,0.5,1)
coef(ricker)
plot(simulate(ricker),var="y")


## ----bifdiag_code,eval=F,echo=T------------------------------------------
## p <- parmat(coef(ricker),500)
## p["r",] <- seq(from=2,to=40,length=500)
## y <- trajectory(ricker,params=p,times=200:1000)
## matplot(p["r",],y["N",,],pch=".",col='black',xlab='r',ylab='N',log='x')


## ----bifdiag,eval=T,echo=F,out.width="8cm"-------------------------------
par(mai=c(0.8,0.5,0.1,0.1))
p <- parmat(coef(ricker),500)
p["r",] <- seq(from=2,to=40,length=500)
y <- trajectory(ricker,params=p,times=200:1000)
matplot(p["r",],y["N",,],pch=".",col='black',xlab='r',ylab='N',log='x')


## ----pfilter1------------------------------------------------------------
pf <- pfilter(ricker,Np=1000)
class(pf)
plot(pf)
logLik(pf)


## ----pfilter2------------------------------------------------------------
pf <- pfilter(pf)
logLik(pf)


## ----pfilter3------------------------------------------------------------
pf <- pfilter(pf,Np=100)
logLik(pf)


## ----parus-data----------------------------------------------------------
dat <- read.csv("parus.csv")
head(dat)
plot(pop~year,data=dat,type='o')


## ----parus-pomp1---------------------------------------------------------
parus <- pomp(dat,times="year",t0=1959)


## ----parus-plot1-code,eval=F,echo=T--------------------------------------
## plot(parus)

## ----parus-plot1,eval=T,echo=F,out.width="12cm"--------------------------
par(mai=c(0.8,0.5,0.1,0.1))
plot(parus)


## ----parus-skel-defn-----------------------------------------------------
skel <- Csnippet("DN = r*N*exp(-N);")


## ----parus-add-skel------------------------------------------------------
parus <- pomp(parus,skeleton=map(skel),statenames="N",paramnames="r")


## ----parus-first-traj,out.width="8cm"------------------------------------
traj <- trajectory(parus,params=c(N.0=1,r=12), format="data.frame")
ggplot(data=traj,aes(x=year,y=N))+geom_line()


## ----parus-first-traj-vector,results='markup',out.width="8cm"------------
parus2 <- pomp(parus,skeleton=vectorfield(skel),statenames="N",paramnames="r")
traj2 <- trajectory(parus2,params=c(N.0=1,r=12),format="data.frame")
ggplot(data=traj2,aes(x=year,y=N))+geom_line()


## ----parus-sim-defn------------------------------------------------------
stochStep <- Csnippet("
  e = rnorm(0,sigma);
  N = r*N*exp(-N+e);
")
pomp(parus,rprocess=discrete.time.sim(step.fun=stochStep,delta.t=1),
     paramnames=c("r","sigma"),statenames=c("N","e")) -> parus


## ----ricker-first-sim_code,eval=F,echo=T---------------------------------
## sim <- simulate(parus,params=c(N.0=1,e.0=0,r=12,sigma=0.5),
##                 format="data.frame")
## plot(N~year,data=sim,type='o')

## ----ricker-first-sim,eval=T,echo=F,out.width="12cm"---------------------
par(mai=c(0.8,0.5,0.1,0.1))
sim <- simulate(parus,params=c(N.0=1,e.0=0,r=12,sigma=0.5),
                format="data.frame")
plot(N~year,data=sim,type='o')


## ----parus-rmeas-defn----------------------------------------------------
rmeas <- Csnippet("pop = rpois(phi*N);")


## ----parus-dmeas-defn----------------------------------------------------
dmeas <- Csnippet("lik = dpois(pop,phi*N,give_log);")


## ----parus-add-meas------------------------------------------------------
pomp(parus,rmeasure=rmeas,dmeasure=dmeas,statenames=c("N"),
  paramnames=c("phi")) -> parus


## ----ricker-add-params,out.width="10cm"----------------------------------
coef(parus) <- c(N.0=1,e.0=0,r=20,sigma=0.1,phi=200)
sims <- simulate(parus,nsim=3,format="data.frame",
  include.data=TRUE)
ggplot(data=sims,mapping=aes(x=year,y=pop))+geom_line()+
  facet_wrap(~.id)

