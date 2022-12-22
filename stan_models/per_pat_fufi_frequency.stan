data {
  int n_pat; // number of patients
  int NCrypts[n_pat];
  int NFufis[n_pat];
}

parameters {
  real<lower=0> mu_nf;
  real<lower=0> sig_nf;
  real<lower=0,upper=1> p_fufi_pat_i[n_pat];
}

model {
  mu_nf ~ std_normal();
  sig_nf ~ std_normal();
  p_fufi_pat_i ~ normal(mu_nf, sig_nf);

  for (i in 1:n_pat) {
    NFufis[i] ~ binomial(NCrypts[i], p_fufi_pat_i[i]);
  }
}
