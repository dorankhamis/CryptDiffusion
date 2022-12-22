data {
   int num_slides;
   int N[num_slides]; // number of crypts per slide
   int n[num_slides]; // number of fufis per slide
}
parameters {
   real<lower=0, upper=1> mean_p; // of proportion population distribution
   real<lower=0> var_p; // of proportion population distribution
   real<lower=0, upper=1> p[num_slides]; // proportion of fufis
}
model {
   mean_p ~ normal(0.005, 0.05);
   var_p  ~ normal(0, 0.05);
   p      ~ normal(mean_p, var_p);
   n      ~ binomial(N, p);
}
