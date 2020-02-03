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



## ----read_data,echo=F----------------------------------------------------
dat <- read.table(file="huron_depth.csv",sep=",",header=TRUE)
dat$Date <- strptime(dat$Date,"%m/%d/%Y")
dat$year <- as.numeric(format(dat$Date, format="%Y"))
dat$month <- as.numeric(format(dat$Date, format="%m"))


## ----plot_data-----------------------------------------------------------
head(dat)
huron_depth <- dat$Average
time <- dat$year + dat$month/12 # Note: we treat December 2011 as time 2012.0, etc
plot(huron_depth~time,type="l")


## ----sarima--------------------------------------------------------------
huron_sarma11x10 <- arima(huron_depth,
   order=c(1,0,1),
   seasonal=list(order=c(1,0,0),period=12)
)
huron_sarma11x10


## ----residuals-----------------------------------------------------------
acf(resid(huron_sarma11x10))


## ----data_subset---------------------------------------------------------
monthly_dat <- subset(dat, month==1)
huron <- monthly_dat$Average
year <- monthly_dat$year
plot(x=year,y=huron,type="l")


## ----h0_fit--------------------------------------------------------------
fit0 <- arima(huron,order=c(1,0,0))
fit0


## ----h1_fit--------------------------------------------------------------
fit1 <- arima(huron,order=c(1,0,0),xreg=year)
fit1

