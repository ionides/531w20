library(pomp)

sir_data <- read.table("bsflu_data.txt")

sir_step <- Csnippet("
  double dN_SI = rbinom(S,1-exp(-Beta*I/N*dt));
  double dN_IR = rbinom(I,1-exp(-gamma*dt));
  S -= dN_SI;
  I += dN_SI - dN_IR;
  R += dN_IR;
  H += dN_IR;
")

sir_skel <- Csnippet("
  double dN_SIdt = Beta*S*I/N;
  double dN_IRdt = gamma*I;
  DS -= dN_SIdt;
  DI += dN_SIdt - dN_IRdt;
  DR += dN_IRdt;
  DH += dN_IRdt;
")

sir_rinit <- Csnippet("
  S = nearbyint(N)-1;
  I = 1;
  R = 0;
  H = 0;
")

sir_dmeas <- Csnippet("lik = dbinom(B,H,rho,give_log);")
sir_rmeas <- Csnippet("B = rbinom(H,rho);")

sir <- pomp(subset(sir_data,select=c(day,B)),
  time="day",t0=0,
  rprocess=euler(sir_step,delta.t=1/6),
  rmeasure=sir_rmeas,
  dmeasure=sir_dmeas,
  rinit=sir_rinit,
  skeleton = vectorfield(sir_skel),
  paramnames=c("Beta","gamma","N","rho"),
  statenames=c("S","I","R","H"),
  zeronames="H") 

if(sir_test <- 1){
  sims <- simulate(sir,params=c(Beta=1.5,gamma=1,rho=0.9,N=2600),
    nsim=20,format="data.frame",include=TRUE)
  ggplot(sims,mapping=aes(x=day,y=B,group=.id,color=.id=="data"))+
    geom_line()+guides(color=FALSE)
}



