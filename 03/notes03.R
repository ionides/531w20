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


## ----echo=F--------------------------------------------------------------
set.seed(2050320976)


## ----opts,include=FALSE,cache=FALSE--------------------------------------
options(
  keep.source=TRUE,
  encoding="UTF-8"
)


## ----stationarity_sim, echo=FALSE----------------------------------------
N <- 500
times <- 1:N
T1 <- 120
T2 <- 37
set.seed(73413)
y <- sin(2*pi*(times/T1 + runif(1))) +   sin(2*pi*(times/T2 + runif(1))) + rnorm(N)
x <- y[1:50]
oldpars <- par(mfrow=c(1,2))
plot(x,ty="l",xlab="")
plot(y,ty="l",xlab="")
par(oldpars)


## ----sinusoidal,echo=FALSE,out.width="10cm"------------------------------
np <- 500
U <- seq(from=0,to=1,length=np)
epsilon1 <- sin(2*pi*U)
epsilon2 <- sin(2*pi*2*U)
epsilon3 <- sin(2*pi*3*U)
matplot(U,cbind(epsilon1,epsilon2,epsilon3),col=c("black","red","blue"),lty=c(1,2,4),ylab="",ty="l",xlab="U")
abline(h=0,lty="dotted")
abline(v=c(1/4,1/2,3/4),lty="dotted")



## ----ar_arima_sim_code,echo=T,eval=F-------------------------------------
## set.seed(123456789)
## ar1 <- arima.sim(list(ar=0.6),n=100,sd=1)
## plot(ar1,type="l")


## ----ar_arima_sim,fig.width=4,fig.height=2.5,echo=F,eval=T,out.width="8cm"----
par(mai=c(0.8,0.8,0.1,0.1))
set.seed(123456789)
ar1 <- arima.sim(list(ar=0.6),n=100,sd=1)
plot(ar1,type="l")


## ----ar_sim_code,echo=T,eval=F-------------------------------------------
## set.seed(123456789)
## N <- 100
## X <- numeric(N)
## X[1] <- rnorm(1,sd=1.56)
## for(n in 2:N) X[n] <- 0.6 * X[n-1] + rnorm(1)
## plot(X,type="l")


## ----ar_sim,fig.width=4,fig.height=2.5,out.width="6cm",echo=F,eval=T-----
par(mai=c(0.5,0.8,0.1,0.1))
set.seed(123456789)
N <- 100
X <- numeric(N)
X[1] <- rnorm(1,sd=1.56)
for(n in 2:N) X[n] <- 0.6 * X[n-1] + rnorm(1)
plot(X,type="l")


## ----ma_sim_code,eval=F,echo=T-------------------------------------------
## N <- 100
## set.seed(123456789)
## X1 <- arima.sim(list(ma=c(1.5,1)),n=N,sd=1)
## set.seed(123456789)
## epsilon <- rnorm(N+2)
## X2 <- numeric(N)
## for(n in 1:N) X2[n] <- epsilon[n+2]+1.5*epsilon[n+1]+epsilon[n]
## plot(X1,type="l") ; plot(X2,type="l")


## ----ma_sim,eval=T,echo=F,fig.width=8,fig.height=2.5,out.width="12cm"----
oldpars <- par(mfrow=c(1,2))
par(mai=c(0.8,0.8,0.1,0.1))
N <- 100
set.seed(123456789)
X1 <- arima.sim(list(ma=c(1.5,1)),n=N,sd=1)
set.seed(123456789)
epsilon <- rnorm(N+2)
X2 <- numeric(N)
for(n in 1:N) X2[n] <- epsilon[n+2]+1.5*epsilon[n+1]+epsilon[n]
plot(X1,type="l") ; plot(X2,type="l")
par(oldpars)


## ----check---------------------------------------------------------------
all(X1==X2)


## ----noise_sim_code,echo=T,eval=F----------------------------------------
## N <- 100
## set.seed(123456789)
## epsilon <- rnorm(N)
## plot(epsilon,type="l")


## ----noise_sim,fig.width=4,fig.height=2.5,out.width="6cm",echo=F,eval=T----
par(mai=c(0.8,0.5,0.1,0.1))
N <- 100
set.seed(123456789)
epsilon <- rnorm(N)
plot(epsilon,type="l")


## ----sp500_code,echo=T,eval=F--------------------------------------------
## dat <- read.table("sp500.csv",sep=",",header=TRUE)
## N <- nrow(dat)
## sp500 <- dat$Close[N:1] # data are in reverse order in sp500.csv
## par(mfrow=c(1,2))
## plot(sp500,type="l") ; plot(log(sp500),type="l")


## ----sp500,echo=F,eval=T,fig.width=6,out.width="8cm"---------------------
dat <- read.table("sp500.csv",sep=",",header=TRUE)
N <- nrow(dat)
sp500 <- dat$Close[N:1] # data are in reverse order in sp500.csv
par(mfrow=c(1,2))
plot(sp500,type="l") ; plot(log(sp500),type="l")


## ----sp500params_code,echo=T,eval=F--------------------------------------
## mu <- mean(diff(log(sp500)))
## sigma <- sd(diff(log(sp500)))
## set.seed(95483123)
## X1 <- log(sp500[1])+cumsum(c(0,rnorm(N-1,mean=mu,sd=sigma)))
## set.seed(324324587)
## X2 <- log(sp500[1])+cumsum(c(0,rnorm(N-1,mean=mu,sd=sigma)))
## par(mfrow=c(1,2))
## plot(X1,type="l") ; plot(X2,type="l")


## ----sp500params,echo=T,eval=F,fig.width=4,fig.height=2.5,out.width="8cm"----
## par(mai=c(0.8,0.5,0.1,0.1))


## ----sp500_acf_code,eval=F,echo=T----------------------------------------
## z <- diff(log(sp500))
## acf(z)


## ----sp500_acf,echo=F,fig.width=4,fig.height=2.5,out.width="8cm"---------
par(mai=c(0.8,0.5,0.1,0.1))
z <- diff(log(sp500))
acf(z)


## ----sp500_abs_return_acf_code,echo=T,eval=F-----------------------------
## acf(abs(z-mean(z)),lag.max=200)


## ----sp500_abs_return_acf,echo=F,eval=T,fig.height=2.5,fig.width=4,out.width="6cm"----
par(mai=c(0.8,0.5,0.1,0.1))
acf(abs(z-mean(z)),lag.max=200)

