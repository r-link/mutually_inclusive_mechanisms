// Initial model generated with brms 2.12.0 and modified manually.
//
// The model code was based on brms to take advantage of the
// automated generation of datasets for model fitting (compare script 05).
// A downside are the unspecific custom names of variables and parameters 
// generated by brms. We hope that the ample comments will make up for it. 
functions {
}
data {
  int<lower=1> N;  // number of observations
  int Y[N];  // response variable
  int<lower=1> K;  // number of population-level effects
  matrix[N, K] X;  // population-level design matrix
  // HSM  effects
  int<lower=1> N_spec; // number of species
  vector[N_spec] hsm_obs;         // observed hydraulic safety margins
  vector<lower=0>[N_spec] sd_hsm; // HSM measurement errors
  // NSC effects
  int<lower=1> N_nsc; // number of repeated measurements of change in relative sugar content
  int<lower=1> J_nsc[N_nsc];  // species indicator for the NSC measurements
  vector[N_nsc] nsc_obs;      // observed change in relative sugar content
  // data for plot-level effects
  int<lower=1> N_1;  // number of grouping levels
  int<lower=1> M_1;  // number of coefficients per level
  int<lower=1> J_1[N];  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_1_1;
  // data for species-level trait effects
  int<lower=1> N_2;  // number of grouping levels
  int<lower=1> M_2;  // number of coefficients per level
  int<lower=1> J_2[N];  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_2_1;
  vector[N] Z_2_2;
  vector[N] Z_2_3;
  vector[N] Z_2_4;
  vector[N] Z_2_5;
  vector[N] Z_2_6;
  int<lower=1> NC_2;  // number of group-level correlations
  // data for neighborhood effects
  int<lower=1> N_3;  // number of grouping levels
  int<lower=1> M_3;  // number of coefficients per level
  int<lower=1> J_3[N];  // grouping indicator per observation
  int<lower=1> N_pairs;  // number of neighbor pairs
  // group-level predictor values
  vector[N] Z_3_1;
  vector[N] Z_3_2;
  vector[N] Z_3_3;
  vector[N] Z_3_4;
  vector[N] Z_3_5;
  vector[N] Z_3_6;
  vector[N] Z_3_7;
  vector[N] Z_3_8;
  vector[N] Z_3_9;
  vector[N] Z_3_10;
  vector[N] Z_3_11;
  vector[N] Z_3_12;
  int pos_mat[N_3, N_3];  // matrix with positions for correlated neighborhood effects
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
  int Kc = K - 1;
  matrix[N, Kc] Xc;  // centered version of X without an intercept
  vector[Kc] means_X;  // column means of X before centering
  for (i in 2:K) {
    means_X[i - 1] = mean(X[, i]);
    Xc[, i - 1] = X[, i] - means_X[i - 1];
  }
}
parameters {
  vector[Kc] b;    // population-level effects
  real Intercept;  // temporary intercept for centered predictors
  real b_nsc;      // parameter for the effect of the change in relative sugar content
  real b_hsm;      // parameter for the effect of hydraulic safety margins
  vector[N_spec] nsc_true; // latent true species level averages of change in relative sugar content
  vector[N_spec] hsm_true; // latent true species level hydraulic safety margins
  real<lower=0> sd_nsc;    // NSC measurement variability
  vector[N_spec] gamma0;   // average neighbourhood effects
  real<lower=0> c; 
  // varying plot effects
  vector<lower=0>[M_1] sd_1;  // group-level standard deviations
  vector[N_1] z_1[M_1];  // standardized group-level effects
  // varying species-level effects
  vector<lower=0>[M_2] sd_2;  // group-level standard deviations
  matrix[M_2, N_2] z_2;  // standardized group-level effects
  cholesky_factor_corr[M_2] L_2;  // cholesky factor of correlation matrix
  // neighbourhood effects
  real<lower=0> sd_3;  // group-level standard deviations
  matrix[2, N_pairs] z_3_tri;  // standardized group-level effects for heterospecific interactions
  cholesky_factor_corr[2] L_3;  // cholesky factor of correlation matrix
  vector[N_3] z_3_dia;          // standardized group-level effects for intraspecific interactions
  real<lower=0, upper=0.5> phi[N_spec]; // species specific misclassification rates
}
transformed parameters {
  // plot effects
  vector[N_1] r_1_1;  // actual group-level effects
  // species-level predictor effects
  matrix[N_2, M_2] r_2;  // actual group-level effects
  // using vectors speeds up indexing in loops
  vector[N_2] r_2_1;
  vector[N_2] r_2_2;
  vector[N_2] r_2_3;
  vector[N_2] r_2_4;
  vector[N_2] r_2_5;
  vector[N_2] r_2_6;

  // neighborhood effects
  matrix[N_pairs, 2] r_3_0;  // correlated neighborhood effects
  matrix[N_3, N_3] r_3;  // matrix of neighborhood effects
  // using vectors speeds up indexing in loops
  vector[N_3] r_3_1;  
  vector[N_3] r_3_2;  
  vector[N_3] r_3_3;
  vector[N_3] r_3_4;
  vector[N_3] r_3_5;
  vector[N_3] r_3_6;
  vector[N_3] r_3_7;
  vector[N_3] r_3_8;
  vector[N_3] r_3_9;
  vector[N_3] r_3_10;
  vector[N_3] r_3_11;  
  vector[N_3] r_3_12;

  // initialize linear predictor term
  vector[N] mu = Intercept + Xc * b;

  // compute plot-level effects
  r_1_1 = (sd_1[1] * (z_1[1]));
  // compute species-level predictor effects
  r_2 = (diag_pre_multiply(sd_2, L_2) * z_2)';
  // get vector form
  r_2_1 = r_2[, 1];
  r_2_2 = r_2[, 2];
  r_2_3 = r_2[, 3];
  r_2_4 = r_2[, 4];
  r_2_5 = r_2[, 5];
  r_2_6 = r_2[, 6];

  // compute correlated neighborhood effects
  r_3_0 = ((sd_3 * L_3) * z_3_tri)';
  // compute reshaped neighborhood effects matrix
  for (i in 1:(N_3 - 1)){
    r_3[i, i] = sd_3 * z_3_dia[i];
    for (j in (i + 1):N_3){
      r_3[i, j] = r_3_0[pos_mat[i, j], 1];
      r_3[j, i] = r_3_0[pos_mat[j, i], 2];
    }
  }
  r_3[N_3, N_3] =  sd_3 * z_3_dia[N_3];
  // get vector form
  r_3_1 = r_3[, 1];
  r_3_2 = r_3[, 2];
  r_3_3 = r_3[, 3];
  r_3_4 = r_3[, 4];
  r_3_5 = r_3[, 5];
  r_3_6 = r_3[, 6];
  r_3_7 = r_3[, 7];
  r_3_8 = r_3[, 8];
  r_3_9 = r_3[, 9];
  r_3_10 = r_3[, 10];
  r_3_11 = r_3[, 11];
  r_3_12 = r_3[, 12];
  
  // add terms to the linear predictor
  for (n in 1:N) {
  mu[n] +=  + hsm_true[J_2[n]] * b_hsm + nsc_true[J_2[n]] * b_nsc +
  r_1_1[J_1[n]] * Z_1_1[n] + 
  r_2_1[J_2[n]] * Z_2_1[n] + r_2_2[J_2[n]] * Z_2_2[n] + r_2_3[J_2[n]] * Z_2_3[n] + r_2_4[J_2[n]] * Z_2_4[n] + r_2_5[J_2[n]] * Z_2_5[n] + r_2_6[J_2[n]] * Z_2_6[n] +  //r_2_7[J_2[n]] * Z_2_7[n] +
  (r_3_1[J_3[n]] + gamma0[1]) * Z_3_1[n] ^ c + (r_3_2[J_3[n]] + gamma0[2]) * Z_3_2[n] ^ c + (r_3_3[J_3[n]] + gamma0[3] )* Z_3_3[n] ^ c + 
  (r_3_4[J_3[n]] + gamma0[4]) * Z_3_4[n] ^ c + (r_3_5[J_3[n]] + gamma0[5]) * Z_3_5[n] ^ c + (r_3_6[J_3[n]] + gamma0[6]) * Z_3_6[n] ^ c + 
  (r_3_7[J_3[n]] + gamma0[7]) * Z_3_7[n] ^ c + (r_3_8[J_3[n]] + gamma0[8]) * Z_3_8[n] ^ c+ (r_3_9[J_3[n]] + gamma0[9]) * Z_3_9[n] ^ c +
  (r_3_10[J_3[n]] + gamma0[10]) * Z_3_10[n] ^ c + (r_3_11[J_3[n]] + gamma0[11]) * Z_3_11[n] ^ c + (r_3_12[J_3[n]] + gamma0[12])* Z_3_12[n] ^ c;
  }
}
model {
  // priors including all constants
  // main effects
  target += normal_lpdf(b | 0, 1.5);
  target += normal_lpdf(b_nsc | 0, 1.5);
  target += normal_lpdf(b_hsm | 0, 1.5);
  target += normal_lpdf(Intercept | 0, 1.5);
  target += normal_lpdf(gamma0 | 0, 1.5);
  target += normal_lpdf(log(c) | 0, 1);
  // NSC measurement model
  target += normal_lpdf(sd_nsc | 0, 0.2)
    - 1 * normal_lccdf(0 | 0, 0.2);
  target += normal_lpdf(nsc_true | 0, 1);  
  target += normal_lpdf(nsc_obs | nsc_true[J_nsc], sd_nsc);  
  // HSM measurement model
  target += normal_lpdf(hsm_obs | hsm_true, sd_hsm);  
  target += normal_lpdf(hsm_true | 0, 1);  
  // varying plot effects
  target += student_t_lpdf(sd_1 | 3, 0, 1.5)
    - 1 * student_t_lccdf(0 | 3, 0, 1.5);
  target += normal_lpdf(z_1[1] | 0, 1);
  // varying species-specific trait effects
  target += student_t_lpdf(sd_2 | 3, 0, 1.5)
    - M_2 * student_t_lccdf(0 | 3, 0, 1.5);
  target += normal_lpdf(to_vector(z_2) | 0, 1);
  target += lkj_corr_cholesky_lpdf(L_2 | 2);
  // varying neighborhood effects
  target += student_t_lpdf(sd_3 | 3, 0, 1.5)
    - 1 * student_t_lccdf(0 | 3, 0, 1.5);
  target += normal_lpdf(z_3_dia | 0, 1);
  target += normal_lpdf(to_vector(z_3_tri) | 0, 1);
  target += lkj_corr_cholesky_lpdf(L_3 | 2);
  // misclassification rates
  target += beta_lpdf(phi | 0.1, 1);
  // likelihood including all constants 
  // (custom likelihood using if-else statement to account for one inflation due
  // to misclassification)
  if (!prior_only) {
   for (n in 1:N){
     if (Y[n] == 1)
       target += log_sum_exp(bernoulli_lpmf(1 | phi[J_2[n]]),
                           bernoulli_lpmf(0 | phi[J_2[n]]) + bernoulli_logit_lpmf(Y[n] | mu[n]));        
     else
       target +=  bernoulli_lpmf(0 | phi[J_2[n]]) + bernoulli_logit_lpmf(Y[n] | mu[n]); 
    }
  }
}
generated quantities {
  // actual population-level intercept
  real  b_Intercept = Intercept - dot_product(means_X, b);
  // log likelihood for LOO CV
  vector[N] log_lik;
  // compute group-level correlations
  real<lower=-1,upper=1> cor_neighbor = multiply_lower_tri_self_transpose(L_3)[2,1];
  corr_matrix[M_2] Cor_2 = multiply_lower_tri_self_transpose(L_2);
  vector<lower=-1,upper=1>[NC_2] cor_2;
  // extract upper diagonal of correlation matrix
  for (k in 1:M_2) {
    for (j in 1:(k - 1)) {
      cor_2[choose(k - 1, 2) + j] = Cor_2[j, k];
    }
  }
  // extract log likelihood (needed for loo package)
  for (n in 1:N) {
    if (Y[n] == 1)
      log_lik[n] = log_sum_exp(bernoulli_lpmf(1 | phi[J_2[n]]),
                               bernoulli_lpmf(0 | phi[J_2[n]]) + bernoulli_logit_lpmf(Y[n] | mu[n]));        
     else
       log_lik[n] = bernoulli_lpmf(0 | phi[J_2[n]]) + bernoulli_logit_lpmf(Y[n] | mu[n]); 
  }
}
