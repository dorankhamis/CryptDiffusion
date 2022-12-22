
functions {
  
  real Poly_interp4(vector y, real x, vector xx, int ii) {
    // Interpolate the value of a vector y (defined on a grid xx) 
    // at a point x, passing through the points xx[ii] -- x[ii+5].
    real a0;
    real a1;
    real a2;
    real a3;
    real a4;
    a0 = (x-xx[ii+1])*(x-xx[ii+2])*(x-xx[ii+3])*(x-xx[ii+4]) / ((xx[ii]-xx[ii+1])*(xx[ii]-xx[ii+2])*(xx[ii]-xx[ii+3])*(xx[ii]-xx[ii+4])) * y[ii];
    a1 = (x-xx[ii])*(x-xx[ii+2])*(x-xx[ii+3])*(x-xx[ii+4]) / ((xx[ii+1]-xx[ii])*(xx[ii+1]-xx[ii+2])*(xx[ii+1]-xx[ii+3])*(xx[ii+1]-xx[ii+4])) * y[ii+1];
    a2 = (x-xx[ii])*(x-xx[ii+1])*(x-xx[ii+3])*(x-xx[ii+4]) / ((xx[ii+2]-xx[ii])*(xx[ii+2]-xx[ii+1])*(xx[ii+2]-xx[ii+3])*(xx[ii+2]-xx[ii+4])) * y[ii+2];
    a3 = (x-xx[ii])*(x-xx[ii+1])*(x-xx[ii+2])*(x-xx[ii+4]) / ((xx[ii+3]-xx[ii])*(xx[ii+3]-xx[ii+1])*(xx[ii+3]-xx[ii+2])*(xx[ii+3]-xx[ii+4])) * y[ii+3];
    a4 = (x-xx[ii])*(x-xx[ii+1])*(x-xx[ii+2])*(x-xx[ii+3]) / ((xx[ii+4]-xx[ii])*(xx[ii+4]-xx[ii+1])*(xx[ii+4]-xx[ii+2])*(xx[ii+4]-xx[ii+3])) * y[ii+4];
    return(a0+a1+a2+a3+a4);
  }

  real integrate(vector f, vector x, int n) {
    // Integrate an array f on a grid x of size n
    int jj;
    real bma;
    real x1;
    real x2;
    real x3;
    real s = 0.0;
    for (i in 2:(n-5)) {
      jj = i-1;
      bma = x[i]-x[i-1];
      x1 = x[i-1] + bma/4.0; x2 = x[i-1] + bma/2.0; x3 = x[i-1] + 3.0*bma/4.0;
      s = s + bma/90.0 * (7.*f[i-1] + 32.0*Poly_interp4(f, x1, x, jj) + 12.0*Poly_interp4(f, x2, x, jj) + 32.0*Poly_interp4(f, x3, x, jj) + 7.0*f[i]);
    }
    jj = n-4;
    for (i in (n-4):n) {
      bma = x[i]-x[i-1];
      x1 = x[i-1] + bma/4.0; x2 = x[i-1] + bma/2.0; x3 = x[i-1] + 3.0*bma/4.0;
      s = s + bma/90.0 * (7.0*f[i-1] + 32.0*Poly_interp4(f, x1, x, jj) + 12.0*Poly_interp4(f, x2, x, jj) + 32.0*Poly_interp4(f, x3, x, jj) + 7.0*f[i]);
    }
    return(s);
  }
  
  vector gamma_on_grid(real gamma_a, int n) {
    vector[n] gamma_r = rep_vector(gamma_a, n);
    return(gamma_r);
  }
  
  vector create_grid(real lower, real upper, int npts) {
    vector[npts] grid;
    real dx = (upper - lower) / (npts - 1.0);
    grid[1] = lower;
    for (i in 2:npts) grid[i] = grid[i-1] + dx;
    return(grid);
  }
  
  vector Gamma_all (vector gamma_a, int N, vector Rin, vector Rout, 
                    real[] theta, int[] nbrhood_index, int n) {
    // compute integral around circle radius R of the total whitespace density
    vector[N] Gamma_out = rep_vector(0., N);
    vector[n] rgrid;
    vector[n] gamma_r;
    for (i in 1:N) {
      rgrid = create_grid(Rin[i], Rout[i], n);
      gamma_r = gamma_on_grid(gamma_a[nbrhood_index[i]], n);
      Gamma_out[i] = theta[i] * integrate(gamma_r .* rgrid, rgrid, n); // f(r,theta) r dr dtheta
    }
    return(Gamma_out);
  }  
  
}

data {
  
  int N; // number of WT+mut patches to fit
  int M; // number of mut patches/neighbourhoods
  
  vector[N] Gamma_dat; // total white space in patches
  vector[N] Rin; // inner radius of patch subsectors
  vector[N] Rout; // outer radius of patch subsectors
  real theta[N]; // angle subtended by patch subsectors
  
  int nbrhood_index [N]; // picks out correct neighbourhood-specific values (maps N<->M)
  int mark_index [M]; // picks out correct mark-specific values (maps M->[1,2], or [1,...,nmarks])
  int n_grid;
  
}

parameters {
  
  vector<lower=0,upper=1>[M] gamma_a; // the white-space fraction at r->inf in neighbourhood of mut patch
  real<lower=0> sd_Gacc;
  real<lower=0> a1_gpop;
  real<lower=0> a2_gpop;
  
}

model {
  
  // priors
  sd_Gacc ~ std_normal();
  a1_gpop ~ gamma(10, 1);
  a2_gpop ~ gamma(10, 1);
  gamma_a ~ beta(a1_gpop, a2_gpop);
  
  // likelihoods
  target += normal_lpdf(Gamma_dat | Gamma_all(gamma_a, N, Rin, Rout, theta, nbrhood_index, n_grid), sd_Gacc);
  
}

// generated quantities {
//   vector[N] log_lik;
//   vector[N] Gamma_gen;
//   vector[N] mu_gen;
//   real gamma_a_pop;
//   mu_gen = Gamma_all(gamma_a, D, am, awt, N, event_times, pat_ages, 
//                      Rin, Rout, theta, nbrhood_index, n_grid, num_events);
//   for (i in 1:N) {
//     Gamma_gen[i] = normal_rng(mu_gen[i], sd_Gacc);
//   }
//   gamma_a_pop = beta_rng(a1_gpop, a2_gpop);
// 
//   for (n in 1:N) {
//     log_lik[n] = normal_lpdf(Gamma_dat[n] | mu_gen[n], sd_Gacc);
//   }
// }
