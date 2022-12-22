
functions {
  
  real u_R (real r, real gamma_a, real Rdiff, real epsilon) {
    return(epsilon * gamma_a / (3.1415926535897932 * Rdiff*Rdiff) * exp(-r*r / (Rdiff*Rdiff)));
  }
  
  real Gamma_perturbation (real gamma_a, real Rdiff, real epsilon,
                           real Rin, real Rout, real theta) {
     // compute integral around circle radius R of the perturbation of whitespace density
     real xin = u_R(Rin, gamma_a, Rdiff, epsilon);
     real xout = u_R(Rout, gamma_a, Rdiff, epsilon);
     return(-0.5*theta*gamma_a*Rdiff*Rdiff * log((1 + xin)/(1  + xout)));
  }

  vector Gamma_total (vector gamma_a, vector Rdiff, real[] epsilon, int N,
                      vector Rin, vector Rout, real[] theta, int[] nbrhood_index) {
    // compute integral around circle radius R of the total whitespace density
    vector[N] Gamma_out = rep_vector(0., N);
    for (i in 1:N) {
      Gamma_out[i] = 0.5 * theta[i] * (Rout[i]*Rout[i] - Rin[i]*Rin[i]) * gamma_a[nbrhood_index[i]] + 
                       Gamma_perturbation(gamma_a[nbrhood_index[i]], Rdiff[nbrhood_index[i]], 
                                          epsilon[nbrhood_index[i]], Rin[i], Rout[i], theta[i]);
    }
    return(Gamma_out);
  }
  
}

data {
  
  int N; // number of WT+mut patches to fit
  int M; // number of mut patches/neighbourhoods
  int S; // number of unique slides
  
  vector[N] Gamma_dat; // total white space in patches
  vector[N] Rin; // inner radius of patch subsectors
  vector[N] Rout; // outer radius of patch subsectors
  real theta[N]; // angle subtended by patch subsectors
  // vector[S] gamma_a_wt30; // readings of white space per crypt in WT30 patches
  real epsilon [M]; // mass/epsilon parameter to density perturbation
  
  int nbrhood_index [N]; // picks out correct neighbourhood-specific values (maps N<->M)
  int slide_index [M]; // picks out correct slide-specific values (maps M<->S)
  int mark_index [M]; // picks out correct mark-specific time scale
  // real min_ages; // minimum ages for patches
  vector[M] pat_ages; // age of patient as hard bound on patch age
  
}

parameters {
  
  vector<lower=0,upper=1>[M] gamma_a; // the white-space fraction at r->inf in neighbourhood of mut patch
  vector<lower=5,upper=100>[M] t_p; // the age of the mutant patch
  real<lower=0> D; // the tissue-intrinsic diffusion coefficient
  real<lower=0> sd_Gacc;
  // vector<lower=0>[S] a_gpat;
  // vector<lower=0>[S] b_gpat;
  real<lower=0> a1_gpop;
  // real<lower=0> b1_gpop;
  real<lower=0> a2_gpop;
  // real<lower=0> b2_gpop;
  
}

transformed parameters {
  
  vector<lower=0>[M] Rdiff; // the diffusion distance of perturbation
  Rdiff = 4 * D * t_p;
  
}

model {
  
  // priors
  sd_Gacc ~ std_normal();
  a1_gpop ~ gamma(10, 1);
  // b1_gpop ~ gamma(1, 1);
  a2_gpop ~ gamma(10, 1);
  // b2_gpop ~ gamma(1, 1);
  // a_gpat ~ gamma(a1_gpop, b1_gpop);
  // b_gpat ~ gamma(a2_gpop, b2_gpop);
  D ~ std_normal();
  gamma_a ~ beta(a1_gpop, a2_gpop);
  for (i in 1:M) {
    // gamma_a[slide_index[i]] ~ beta(a_gpat[slide_index[i]], b_gpat[slide_index[i]]);
    t_p[i] ~ normal(pat_ages[i]/2, pat_ages[i]/4);
  }

  // likelihoods
  // gamma_a_wt30 ~ beta(a_gpat, b_gpat);
  Gamma_dat ~ normal(Gamma_total(gamma_a, Rdiff, epsilon, N, Rin, Rout, theta, nbrhood_index), sd_Gacc);
  
}

generated quantities {
  
  // vector[N+S] log_lik;
  vector[N] log_lik;
  vector[N] Gamma_gen;
  vector[N] mu_gen;
  // vector[S] gamma_a_pat;
  real gamma_a_pop;
  mu_gen = Gamma_total(gamma_a, Rdiff, epsilon, N,
                       Rin, Rout, theta, nbrhood_index);
  for (i in 1:N) {
    Gamma_gen[i] = normal_rng(mu_gen[i], sd_Gacc);
  }
  // for (i in 1:S) {
  //   gamma_a_pat[i] = beta_rng(a_gpat[i], b_gpat[i]);
  // }
  // gamma_a_pop = beta_rng(gamma_rng(a1_gpop, b1_gpop), gamma_rng(a2_gpop, b2_gpop));
  gamma_a_pop = beta_rng(a1_gpop, a2_gpop);

  for (n in 1:N) {
    log_lik[n] = normal_lpdf(Gamma_dat[n] | mu_gen[n], sd_Gacc);
  }
  // for (n in 1:S) {
  //   log_lik[N+n] = beta_lpdf(gamma_a_wt30[n] | a_gpat[n], b_gpat[n]);
  // }
  
}
