
library(pomp)

bsflu_statenames <- c("S","I","R1","R2")
bsflu_paramnames <- c("Beta","mu_IR","rho","mu_R1","mu_R2")

bsflu_dmeasure <- "
  lik = dpois(B,rho*R1+1e-6,give_log);
"

bsflu_rmeasure <- "
  B = rpois(rho*R1+1e-10);
  C = rpois(rho*R2);
"

bsflu_rprocess <- "
  double t1 = rbinom(S,1-exp(-Beta*I*dt));
  double t2 = rbinom(I,1-exp(-dt*mu_IR));
  double t3 = rbinom(R1,1-exp(-dt*mu_R1));
  double t4 = rbinom(R2,1-exp(-dt*mu_R2));
  S -= t1;
  I += t1 - t2;
  R1 += t2 - t3;
  R2 += t3 - t4;
"

bsflu_rinit <- "
 S=762;
 I=1;
 R1=0;
 R2=0;
"

bsflu2 <- pomp(
  data=bsflu_data,
  times="day",
  t0=0,
  rprocess=euler.sim(
    step.fun=Csnippet(bsflu_rprocess),
    delta.t=1/12
  ),
  rmeasure=Csnippet(bsflu_rmeasure),
  dmeasure=Csnippet(bsflu_dmeasure),
  partrans=parameter_trans(log=c("Beta","mu_IR","mu_R1","mu_R2"),logit="rho"),
  statenames=bsflu_statenames,
  paramnames=bsflu_paramnames,
  rinit=Csnippet(bsflu_initializer)
)


