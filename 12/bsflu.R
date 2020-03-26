library(pomp)
library(ggplot2)
library(dplyr)

sir_data <- read.table("bsflu_data.txt")

sir_step <- Csnippet("
  double dN_SI = rbinom(S,1-exp(-Beta*I/N*dt));
  double dN_IR = rbinom(I,1-exp(-mu_IR*dt));
  S -= dN_SI;
  I += dN_SI - dN_IR;
  R += dN_IR;
  H += dN_IR;
")

sir_skel <- Csnippet("
  double dN_SIdt = Beta*S*I/N;
  double dN_IRdt = mu_IR*I;
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
  paramnames=c("Beta","mu_IR","N","rho"),
  statenames=c("S","I","R","H"),
  accumvars="H",
  partrans=parameter_trans(log=c("Beta","mu_IR","N"),logit="rho")
) 

sir_test <- 2

if(sir_test==1){
  sims <- simulate(sir,params=c(Beta=1.5,mu_IR=1,rho=0.9,N=2600),
    nsim=20,format="data.frame",include=TRUE)
  ggplot(sims,mapping=aes(x=day,y=B,group=.id,color=.id=="data"))+
    geom_line()+guides(color=FALSE)
}

if(sir_test==2){
 coef(sir) <- c(Beta=1.5,mu_IR=1,rho=0.9,N=2600)
 x <- trajectory(sir) 
 y <- cbind(as.data.frame(sir),x=x["cases",1,])
 mutate(y,xlab=sprintf("x[%d]",time),
       ylab=sprintf("y[%d]",time)) -> y

 ggplot(data=y,
       mapping=aes(x=time,xend=time))+
  geom_point(aes(y=reports),color='black',alpha=0.5)+
  geom_point(aes(y=x),color='red',alpha=0.5)+
  geom_line(aes(y=reports),color='black',alpha=0.5)+
  geom_line(aes(y=x),color='red',alpha=0.5)+
  geom_text(aes(y=reports,label=ylab,vjust=ifelse(time>=10,2,-1)),
    parse=TRUE,color='black')+
  geom_text(aes(y=x,label=xlab,vjust=ifelse(time>=10,-1,2)),
    parse=TRUE,color='red')+
  geom_segment(aes(y=x,yend=reports),color='blue',linetype=2,alpha=0.3,
               arrow=grid::arrow(length=grid::unit(0.02,"npc")))+
  expand_limits(y=c(-20,320))+
  labs(y="")
}

