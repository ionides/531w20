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



## ----eigen_code,echo=T,eval=F-------------------------------------------------
## N <- 100;  phi <- 0.8;  sigma <- 1
## V <- matrix(NA,N,N)
## for(m in 1:N) for(n in 1:N) V[m,n]<-sigma^2*phi^abs(m-n)/(1-phi^2)
## V_eigen <- eigen(V,symmetric=TRUE)
## matplot(V_eigen$vectors[,1:5],type="l")
## matplot(V_eigen$vectors[,6:9],type="l")


## ----eigen,echo=F,fig.width=6,fig.height=2.5,out.width="12cm"-----------------
oldpars <- par(mfrow=c(1,2))
par(mai=c(0.8,1.2,0.1,0.1))
N <- 100;  phi <- 0.8;  sigma <- 1
V <- matrix(NA,N,N)
for(m in 1:N) for(n in 1:N) V[m,n]<-sigma^2*phi^abs(m-n)/(1-phi^2)
V_eigen <- eigen(V,symmetric=TRUE)
matplot(V_eigen$vectors[,1:5],type="l")
matplot(V_eigen$vectors[,6:9],type="l")
par(oldpars)


## ----evals--------------------------------------------------------------------
round(V_eigen$values[1:9],2)


## ----weather_data_file,eval=F,echo=F------------------------------------------
## system("head ann_arbor_weather.csv",intern=TRUE)


## ----weather_data-------------------------------------------------------------
y <- read.table(file="ann_arbor_weather.csv",header=TRUE)
head(y[,1:9],3)


## ----replace_na---------------------------------------------------------------
low <- y$Low
low[is.na(low)] <- mean(low, na.rm=TRUE)


## ----periodogram--------------------------------------------------------------
spectrum(low, main="Unsmoothed periodogram")


## ----smoothed_periodogram_code,echo=T,eval=F----------------------------------
## spectrum(low, spans=c(3,5,3), main="Smoothed periodogram",
##   ylim=c(15,100))


## ----smoothed_periodogram,echo=F,fig.width=6,fig.height=3,out.width="12cm"----
par(mai=c(0.8,0.8,0.5,0.1))
spectrum(low, spans=c(3,5,3), main="Smoothed periodogram",
  ylim=c(15,100))


## ----taper_plot_code,echo=T,eval=F--------------------------------------------
## plot(spec.taper(rep(1,100)),type="l",
##   main="Default taper in R, for a time series of length 100")
## abline(v=c(10,90),lty="dotted",col="red")


## ----taper_plot,echo=F,fig.width=5,fig.height=2.5,out.width="8cm"-------------
par(mai=c(0.4,0.8,0.5,0.4)) 
plot(spec.taper(rep(1,100)),type="l",
  main="Default taper in R, for a time series of length 100")
abline(v=c(10,90),lty="dotted",col="red") 


## ----ar_periodogram_code,eval=F,echo=T----------------------------------------
## spectrum(low,method="ar",
##   main="Spectrum estimated via AR model picked by AIC")


## ----ar_periodogram,echo=F,fig.width=5,fig.height=2.5,out.width="11cm"--------
par(mai=c(0.4,0.8,0.5,0.4))
spectrum(low,method="ar",
  main="Spectrum estimated via AR model picked by AIC")


## ----poly_fit-----------------------------------------------------------------
lm0 <- lm(Low~1,data=y)
lm1 <- lm(Low~Year,data=y)
lm2 <- lm(Low~Year+I(Year^2),data=y)
lm3 <- lm(Low~Year+I(Year^2)+I(Year^3),data=y)
poly_aic <- matrix( c(AIC(lm0),AIC(lm1),AIC(lm2),AIC(lm3)), nrow=1,
   dimnames=list("<b>AIC</b>", paste("order",0:3)))
require(knitr)
kable(poly_aic,digits=1)


## ----plot_jan_temp,fig.width=5------------------------------------------------
plot(Low~Year,data=y,type="l")


## ----read_glob_temp-----------------------------------------------------------
Z <- read.table("Global_Temperature.txt",header=TRUE)
global_temp <- Z$Annual[Z$Year %in% y$Year]
lm_global <- lm(Low~global_temp,data=y)
AIC(lm_global)


## ----glob_temp_fit------------------------------------------------------------
summary(lm_global)

