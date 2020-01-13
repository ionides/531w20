## ----setup,echo=F,results=F,cache=F--------------------------------------
library(broman) # used for myround 


## ----echo=F--------------------------------------------------------------
set.seed(2050320976)


## ----opts,include=FALSE,cache=FALSE--------------------------------------
options(
  keep.source=TRUE,
  encoding="UTF-8"
)


## ----data----------------------------------------------------------------
global_temp <- read.table("Global_Temperature.txt",header=TRUE)
str(global_temp)


## ----data_plot_code,echo=T,eval=F----------------------------------------
## plot(Annual~Year,data=global_temp,ty="l")


## ----data_plot,echo=F,fig.width=4,fig.height=3.2,out.width="4in",cache=F----
par(mai=c(1.1,1,0.4,1))
plot(Annual~Year,data=global_temp,ty="l")


## ----glob_temp_lm--------------------------------------------------------
lm_fit <- lm(Annual~Year+I(Year^2),data=global_temp)


## ----glob_temp_summary---------------------------------------------------
summary(lm_fit)


## ----glob_temp_lm_plot_code,echo=T,eval=F--------------------------------
## yr <- 1880:2026
## Z <- cbind(1,yr,yr^2)
## beta <- coef(lm_fit)
## prediction <- Z%*%beta
## plot(Annual~Year,data=global_temp,ty="l",xlim=range(yr),
##   ylim=range(c(global_temp$Annual,prediction),na.rm=TRUE),
##   lty="dashed")
## lines(x=yr,y=prediction,col="red")


## ----glob_temp_lm_plot ,fig.width=4,fig.height=2,out.width="8cm",echo=F----
par(mai=c(0.8,1,0.1,1))
yr <- 1880:2026
Z <- cbind(1,yr,yr^2)
beta <- coef(lm_fit)
prediction <- Z%*%beta
plot(Annual~Year,data=global_temp,ty="l",xlim=range(yr),
  ylim=range(c(global_temp$Annual,prediction),na.rm=TRUE),
  lty="dashed")
lines(x=yr,y=prediction,col="red")


## ----acf_global_temp-----------------------------------------------------
acf(resid(lm_fit))


## ----eval=F,echo=T-------------------------------------------------------
## stats:::plot.acf


## ----eval=F,echo=T-------------------------------------------------------
## clim0 <- if (with.ci) qnorm((1 + ci)/2)/sqrt(x$n.used)

