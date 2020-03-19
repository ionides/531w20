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



## ----flu-data1-----------------------------------------------------------
bsflu <- read.table("bsflu_data.txt")
head(bsflu)








## ----rproc2--------------------------------------------------------------
sir_step <- Csnippet("
  double dN_SI = rbinom(S,1-exp(-Beta*I/N*dt));
  double dN_IR = rbinom(I,1-exp(-gamma*dt));
  S -= dN_SI;
  I += dN_SI - dN_IR;
  R += dN_IR;
  H += dN_IR;
")

sir_rinit <- Csnippet("
  S = N-1;
  I = 1;
  R = 0;
  H = 0;
")

pomp(sir,rprocess=euler.sim(sir_step,delta.t=1/6),rinit=sir_rinit,
     paramnames=c("Beta","gamma","N"),statenames=c("S","I","R","H")) -> sir


## ----zero1---------------------------------------------------------------
pomp(sir,zeronames="H") -> sir


## ----meas-model----------------------------------------------------------
dmeas <- Csnippet("lik = dbinom(B,H,rho,give_log);")
rmeas <- Csnippet("B = rbinom(H,rho);")


## ----sir_sim-------------------------------------------------------------
sims <- simulate(sir,params=c(Beta=1.5,gamma=1,rho=0.9,N=2600),
                 nsim=20,format="data.frame",include=TRUE)

ggplot(sims,mapping=aes(x=time,y=B,group=sim,color=sim=="data"))+
  geom_line()+guides(color=FALSE)


## ----seir-diagram,echo=FALSE,cache=FALSE---------------------------------
require(DiagrammeR)
DiagrammeR("graph LR; S(S) --> E; E(E) --> I; I(I) --> R(R);"
           ,height=200,width=600)


## ----bsflu-plot2---------------------------------------------------------
require(reshape2)
ggplot(data=melt(bsflu,id="day"),mapping=aes(x=day,y=value,color=variable))+
  geom_line()+geom_point()


## ----sirr-diagram,echo=FALSE,cache=FALSE---------------------------------
require(DiagrammeR)
DiagrammeR("graph LR; S(S) --> I; I(I) --> R1(R1); R1 --> R2(R2);"
           ,height=200,width=600)

