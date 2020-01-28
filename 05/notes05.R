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



## ----chi_squared--------------------------------------------------------------
qchisq(0.95,df=1)


## ----read_data----------------------------------------------------------------
dat <- read.table(file="huron_depth.csv",sep=",",header=TRUE)
dat$Date <- strptime(dat$Date,"%m/%d/%Y")
dat$year <- as.numeric(format(dat$Date, format="%Y"))
dat$month <- as.numeric(format(dat$Date, format="%m"))
head(dat,3)


## ----select_annual_code,echo=T,eval=F-----------------------------------------
## dat <- subset(dat,month==1)
## huron_depth <- dat$Average
## year <- dat$year
## plot(huron_depth~year,type="l")


## ----select_annual,echo=F,eval=T----------------------------------------------
dat <- subset(dat,month==1)
huron_depth <- dat$Average
year <- dat$year
plot(huron_depth~year,type="l")


## ----aic_table_code,echo=T,eval=F---------------------------------------------
## aic_table <- function(data,P,Q){
##   table <- matrix(NA,(P+1),(Q+1))
##   for(p in 0:P) {
##     for(q in 0:Q) {
##        table[p+1,q+1] <- arima(data,order=c(p,0,q))$aic
##     }
##   }
##   dimnames(table) <- list(paste("AR",0:P, sep=""),paste("MA",0:Q,sep=""))
##   table
## }
## huron_aic_table <- aic_table(huron_depth,4,5)
## require(knitr)
## kable(huron_aic_table,digits=2)


## ----aic_table,echo=F,eval=T--------------------------------------------------
aic_table <- function(data,P,Q){
  table <- matrix(NA,(P+1),(Q+1))
  for(p in 0:P) {
    for(q in 0:Q) {
       table[p+1,q+1] <- arima(data,order=c(p,0,q))$aic
    }
  }
  dimnames(table) <- list(paste("AR",0:P, sep=""),paste("MA",0:Q,sep=""))
  table
}
huron_aic_table <- aic_table(huron_depth,4,5)
require(knitr)
kable(huron_aic_table,digits=2)


## ----arma21fit----------------------------------------------------------------
huron_arma21 <- arima(huron_depth,order=c(2,0,1))
huron_arma21


## ----huron_roots--------------------------------------------------------------
AR_roots <- polyroot(c(1,-coef(huron_arma21)[c("ar1","ar2")]))
AR_roots


## ----huron_profile_code,echo=T,eval=F-----------------------------------------
## K <- 500
## ma1 <- seq(from=0.2,to=1.1,length=K)
## profile_loglik <- rep(NA,K)
## for(k in 1:K){
##    profile_loglik[k] <- logLik(arima(huron_depth,order=c(2,0,1),
##       fixed=c(NA,NA,ma1[k],NA)))
## }
## plot(profile_loglik~ma1,ty="l")


## ----huron_profile,echo=F,fig.width=5,fig.height=2.5,out.width="9cm"----------
par(mai=c(0.8,0.8,0.1,0.1))
K <- 500
ma1 <- seq(from=0.2,to=1.1,length=K)
profile_loglik <- rep(NA,K)
for(k in 1:K){
   profile_loglik[k] <- logLik(arima(huron_depth,order=c(2,0,1),
      fixed=c(NA,NA,ma1[k],NA)))
}
plot(profile_loglik~ma1,ty="l")


## ----simA_code,echo=T,eval=F--------------------------------------------------
## set.seed(57892330)
## J <- 1000
## params <- coef(huron_arma21)
## ar <- params[grep("^ar",names(params))]
## ma <- params[grep("^ma",names(params))]
## intercept <- params["intercept"]
## sigma <- sqrt(huron_arma21$sigma2)
## theta <- matrix(NA,nrow=J,ncol=length(params),
##    dimnames=list(NULL,names(params)))
## for(j in 1:J){
##    Y_j <- arima.sim(
##       list(ar=ar,ma=ma),
##       n=length(huron_depth),
##       sd=sigma
##    )+intercept
##    theta[j,] <- coef(arima(Y_j,order=c(2,0,1)))
## }
## hist(theta[,"ma1"],freq=FALSE)


## ----simA,echo=F,eval=T,fig.width=5,fig.height=2.5,out.width="10cm"-----------
par(mai=c(0.8,0.8,0.5,0.1))
set.seed(57892330)
J <- 1000
params <- coef(huron_arma21)
ar <- params[grep("^ar",names(params))]
ma <- params[grep("^ma",names(params))]
intercept <- params["intercept"]
sigma <- sqrt(huron_arma21$sigma2)
theta <- matrix(NA,nrow=J,ncol=length(params),
   dimnames=list(NULL,names(params)))
for(j in 1:J){
   Y_j <- arima.sim(
      list(ar=ar,ma=ma),
      n=length(huron_depth),
      sd=sigma
   )+intercept
   theta[j,] <- coef(arima(Y_j,order=c(2,0,1)))
}
hist(theta[,"ma1"],freq=FALSE) 


## ----density_code,eval=F,echo=T-----------------------------------------------
## plot(density(theta[,"ma1"],bw=0.05))


## ----density,echo=F,eval=T,fig.width=5,fig.height=2.5,out.width="10cm"--------
par(mai=c(0.8,0.8,0.5,0.1))
plot(density(theta[,"ma1"],bw=0.05))


## ----range--------------------------------------------------------------------
range(theta[,"ma1"])


## ----parallel-setup,cache=FALSE-----------------------------------------------
library(doParallel)
registerDoParallel()


## ----simB---------------------------------------------------------------------
J <- 1000
huron_ar1 <- arima(huron_depth,order=c(1,0,0))
params <- coef(huron_ar1)
ar <- params[grep("^ar",names(params))]
intercept <- params["intercept"]
sigma <- sqrt(huron_ar1$sigma2)
t1 <- system.time(
  huron_sim <- foreach(j=1:J) %dopar% {
     Y_j <- arima.sim(list(ar=ar),n=length(huron_depth),sd=sigma)+intercept
     try(coef(arima(Y_j,order=c(2,0,1))))
  }
) 


## ----out, cache=FALSE---------------------------------------------------------
sum(sapply(huron_sim, function(x) inherits(x,"try-error"))) 


## ----histB, cache=FALSE-------------------------------------------------------
ma1 <- unlist(lapply(huron_sim,function(x)
   if(!inherits(x,"try-error"))x["ma1"] else NULL ))
hist(ma1,breaks=50)  


## ----repeated_aic,echo=FALSE--------------------------------------------------
require(knitr)
kable(huron_aic_table,digits=1)

