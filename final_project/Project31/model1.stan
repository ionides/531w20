data {
int<lower=0> N; // number of data items
 // number of predictors
vector[N] x; // predictor matrix
real y[N]; // outcome vector
}
parameters {
 // intercept
real beta; 
real alpha;
real sigma;// coefficients for predictors
// error scale
}
model {
 
      beta ~ normal(0,10);
alpha ~normal(0,1);
sigma~ inv_gamma(1, 1);
y ~ normal(alpha+x * beta, sigma); // likelihood
   

}

