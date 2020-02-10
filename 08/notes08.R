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



## ----data_unadj----------------------------------------------------------
system("head unadjusted_unemployment.csv",intern=TRUE)
U1 <- read.table(file="unadjusted_unemployment.csv",
  sep=",",header=TRUE)
head(U1,3)


## ----reshape_code,echo=T,eval=F------------------------------------------
## u1 <- t(as.matrix(U1[2:13]))
## dim(u1) <- NULL
## date <- seq(from=1948,length=length(u1),by=1/12)
## plot(date,u1,type="l",ylab="Percent unemployment (unadjusted)")


## ----reshape,echo=F,out.width="10cm"-------------------------------------
u1 <- t(as.matrix(U1[2:13]))
dim(u1) <- NULL
date <- seq(from=1948,length=length(u1),by=1/12)
plot(date,u1,type="l",ylab="Percent unemployment (unadjusted)")


## ----data_adj_code,echo=F,eval=T,out.width="10cm"------------------------
U2 <- read.table(file="adjusted_unemployment.csv",sep=",",header=TRUE)
u2 <- t(as.matrix(U2[2:13]))
dim(u2) <- NULL
plot(date,u1,type="l",ylab="percent",col="black")
lines(date,u2,type="l",col="red")
title("Unemployment. Raw (black) and seasonally adjusted (red)")


## ----adjustment_spectrum_code,eval=F,echo=T------------------------------
## u1_ts <- ts(u1,start=1948,frequency=12)
## u2_ts <- ts(u2,start=1948,frequency=12)
## spectrum(ts.union(u1_ts,u2_ts),spans=c(3,5,3),
##   main="Unemployment. Raw (black) and seasonally adjusted (red)")


## ----adjustment_spectrum,eval=T,echo=F,out.width="10cm"------------------
u1_ts <- ts(u1,start=1948,frequency=12)
u2_ts <- ts(u2,start=1948,frequency=12)
spectrum(ts.union(u1_ts,u2_ts),spans=c(3,5,3),
  main="Unemployment. Raw (black) and seasonally adjusted (red)")


## ----plot.ts_code,echo=T,eval=F------------------------------------------
## plot(u1_ts)


## ----plot.ts,echo=F,out.width="10cm"-------------------------------------
par(mai=c(0.8,0.8,0.1,0.1))
plot(u1_ts)


## ----bls_filter----------------------------------------------------------
s <- spectrum(ts.union(u1_ts,u2_ts),plot=FALSE)


## ----s_names-------------------------------------------------------------
names(s)


## ----s_transfer_code,eval=F,echo=T---------------------------------------
## plot(s$freq,s$spec[,2]/s$spec[,1],type="l",log="y",
##   ylab="frequency ratio", xlab="frequency",
##   main="frequency response (dashed lines at 0.9 and 1.1)")
## abline(h=c(0.9,1.1),lty="dashed",col="red")


## ----s_transfer,eval=T,echo=F,out.width="10cm"---------------------------
plot(s$freq,s$spec[,2]/s$spec[,1],type="l",log="y",
  ylab="frequency ratio", xlab="frequency",  
  main="frequency response (dashed lines at 0.9 and 1.1)")
abline(h=c(0.9,1.1),lty="dashed",col="red")


## ----loess_code,echo=T,eval=F--------------------------------------------
## u1_loess <- loess(u1~date,span=0.5)
## plot(date,u1,type="l",col="red")
## lines(u1_loess$x,u1_loess$fitted,type="l")


## ----loess,echo=F,out.width="10cm"---------------------------------------

## ----loess_code----------------------------------------------------------
u1_loess <- loess(u1~date,span=0.5)
plot(date,u1,type="l",col="red")
lines(u1_loess$x,u1_loess$fitted,type="l")


## ----loess_transfer_code,echo=T,eval=F-----------------------------------
## s2 <- spectrum(ts.union(
##   u1_ts,ts(u1_loess$fitted,start=1948,frequency=12)),
##   plot=FALSE)
## plot(s2$freq,s2$spec[,2]/s$spec[,1],type="l",log="y",
##   ylab="frequency ratio", xlab="frequency", xlim=c(0,1.5),
##   main="frequency response (dashed line at 1.0)")
## abline(h=1,lty="dashed",col="red")


## ----loess_transfer,eval=T,echo=F,out.width="10cm"-----------------------
s2 <- spectrum(ts.union(
  u1_ts,ts(u1_loess$fitted,start=1948,frequency=12)),
  plot=FALSE)
plot(s2$freq,s2$spec[,2]/s$spec[,1],type="l",log="y",
  ylab="frequency ratio", xlab="frequency", xlim=c(0,1.5),
  main="frequency response (dashed line at 1.0)")
abline(h=1,lty="dashed",col="red")


## ----cycles_code,echo=T,eval=F-------------------------------------------
## u_low <- ts(loess(u1~date,span=0.5)$fitted,
##   start=1948,frequency=12)
## u_hi <- ts(u1 - loess(u1~date,span=0.1)$fitted,
##   start=1948,frequency=12)
## u_cycles <- u1 - u_hi - u_low
## plot(ts.union(u1, u_low,u_hi,u_cycles),
##   main="Decomposition of unemployment as trend + noise + cycles")


## ----cycles,echo=F,eval=T,fig.width=7,fig.height=5,out.width="10cm"------
par(mai=c(0.8,0.8,0.5,0.1))
u_low <- ts(loess(u1~date,span=0.5)$fitted,
  start=1948,frequency=12)
u_hi <- ts(u1 - loess(u1~date,span=0.1)$fitted,
  start=1948,frequency=12)
u_cycles <- u1 - u_hi - u_low
plot(ts.union(u1, u_low,u_hi,u_cycles),
  main="Decomposition of unemployment as trend + noise + cycles")


## ----freq_response,echo=F,out.width="10cm"-------------------------------
spec_cycle <- spectrum(ts.union(u1_ts,u_cycles),
  spans=c(3,3),
  plot=FALSE)
freq_response_cycle <- spec_cycle$spec[,2]/spec_cycle$spec[,1]
plot(spec_cycle$freq,freq_response_cycle,
  type="l",log="y",
  ylab="frequency ratio", xlab="frequency", xlim=c(0,1.2), ylim=c(5e-6,1.1),
  main="frequency response (dashed line at 1.0)")
abline(h=1,lty="dashed",col="red")  



## ----show_range,echo=F,out.width="10cm"----------------------------------
cut_fraction <- 0.5
plot(spec_cycle$freq,freq_response_cycle,
  type="l",log="y",
  ylab="frequency ratio", xlab="frequency", xlim=c(0,0.9), ylim=c(1e-4,1.1),
  main=paste("frequency response, showing region for ratio >", cut_fraction))
abline(h=1,lty="dashed",col="blue")  
freq_cycles <- range(spec_cycle$freq[freq_response_cycle>cut_fraction]) 
abline(v=freq_cycles,lty="dashed",col="blue") 
abline(h=cut_fraction,lty="dashed",col="blue")


## ----print_range---------------------------------------------------------
kable(matrix(freq_cycles,nrow=1,
  dimnames=list("frequency",c("low","hi"))),digits=3)


## ----zoomed_spectrum,echo=F,fig.width=6,fig.height=3,out.width="10cm"----
s1 <- spectrum(u1_ts,spans=c(3),plot=FALSE)
par(mai=c(1,0.8,0.1,0.1))
plot(s1,xlim=c(0,0.7),ylim=c(1e-2,max(s1$spec)),main="")

