
data {
  
  //Observation level data
  int<lower=1> N; //Number of observations
  int<lower=1> K; //Number of predictor variables at observation level
  vector[N] y; //Growth value of observation i
  matrix[N,K] X; //Matrix of predictors, including intercept
  int<lower=1> j[N]; //Indictor for individual tree that each observation belongs to
  
  int<lower=1> Nsp; //Number of species at observation level
  int sp_obs[N]; //Species indicator for each observation
  
  //Individual-level data
  int<lower=1> J; //Number of target individuals
  int<lower=1> J_1; //Number of target trees in species 1 (Jacarnada)
  int<lower=1> J_2; //Number of target trees in species 2 (Luehea)
  int<lower=1> J_3; //Number of target trees in species 3 (Cecropia)
  int<lower=1> J_1_end; //Last target tree for species 1
  int<lower=1> J_2_end; //Last target tree for species 2
  int<lower=1> J_3_end; //Last target tree for species 3
  int<lower=1> J_2_start; //First target tree for species 2
  int<lower=1> J_3_start; //First target tree for species 3
  
  int<lower=1> L; //Number of predictor variables at individual tree level
  int<lower=1> q[J]; //Indicator for quadrat that each individual tree belongs to
  int<lower=1> q_1[J_1]; //Indicator for quadrat that each individual tree belongs to in species 1 (Jacaranda)
  int<lower=1> q_2[J_2]; //Indicator for quatrat that each individual tree belongs to in species 2 (Luehea)
  int<lower=1> q_3[J_3]; //Indicator for quatrat that each individual tree belongs to in species 3 (Cecropia)
  vector[J] U; //Vector to hold intercept and inbreeding predictor values for individual-level regression
  
  int sp_tree[J]; //Species indicator for each tree
  
  //Quadrat-level data
  int<lower=1> Q; //Number of quadrats
  int<lower=1> M; //Number of predictor variables at quadrat level
  matrix[Q,M] V; //Matrix of predictors for regression at quadrat-level

}

parameters {
  //Observation-level parameters
  vector[J] A; //Vector to hold alpha parameters for observation-level regression
  matrix[K,Nsp] B; //Matrix to hold beta parameters for observation-level regression (both species in same vector being indexed separately by species)
  
  //Individual-level parameters
  vector[Q] theta; //Vector to hold theta parameters for individual-level regression (both species in same vector being indexed separately by species)
  row_vector[Nsp] G; //Vector to hold gamma values for individual-level regression (don't have to define number of rows in G because only 1 predictor)
  
  //Quadrat-level predictors
  vector[M] D; //Vector to hold raw delta parameters for quadrat-level regression (both species get the same parameters here)
  
  //Hyperparameters
  // vector[K] mu_beta; //Vector to hold mu_beta hyperparameters (1 for each species)
  vector<lower=0>[K] sigma_beta; //Vector to hold sigma_beta hyperparameters (1 for each species)
  real mu_gamma; //Hyperparameter for mean of gammas
  real<lower=0> sigma_gamma; //Hyperparameters for variance of gammas
  vector[M] mu_delta; //Vector to hold mu_delta hyperparameters
  vector<lower=0>[M] sigma_delta; //Vector to hold sigma_delta hyperparameters
  
  //Variance parameters
  real<lower=0> sigma_y;
  real<lower=0> sigma_alpha;
  real<lower=0> sigma_theta;
}

model {
  
  //Priors
  B[1,] ~ normal(0, sigma_beta[1]);
  B[2,] ~ normal(0, sigma_beta[2]);
  B[3,] ~ normal(0, sigma_beta[3]);
  B[4,] ~ normal(0, sigma_beta[4]);
  B[5,] ~ normal(0, sigma_beta[5]);
  // mu_beta ~ normal(0,3);
  sigma_beta ~ cauchy(0,3);
  
  G ~ normal(mu_gamma, sigma_gamma);
  mu_gamma ~ normal(0,3);
  sigma_gamma ~ cauchy(0,3);
  
  D ~ normal(mu_delta, sigma_delta);
  mu_delta ~ normal(0,3);
  sigma_delta ~ cauchy(0,3);
  
  sigma_y ~ cauchy(0,3);
  sigma_alpha ~ cauchy(0,3);
  sigma_theta ~ cauchy(0,3);
    
  //Likelihood statements
  
  //Observation level
  for (i in 1:N) {
    y[i] ~ normal(A[j[i]] + (X[i] * B[,sp_obs[i]]), sigma_y); //Regression at observation level
  }
  
  //Individual level
  for (r in 1:J_1_end) { //Loop over trees in species 1 (Jacaranda)
    A[r] ~ normal(theta[q[r]] + (G[1] * U[r]), sigma_alpha); //Regression at individual tree level
  }
  for (t in J_2_start:J_2_end) { //Loop over trees in species 2 (Cecropia)
    A[t] ~ normal(theta[q[t]] + (G[2] * U[t]), sigma_alpha);
  }
  for (u in J_3_start:J_3_end) { //Loop over trees in species 2 (Triplaris)
    A[u] ~ normal(theta[q[u]] + (G[3] * U[u]), sigma_alpha);
  }
  
  //Quadrat level
  for (s in 1:Q) {
    theta[s] ~ normal( (V[s] * D), sigma_theta); //Regression at quadrat level independent of species
  }
}

generated quantities{
  vector[N] y_pred;
  vector[N] log_lik;
  for (i in 1:N) {
    y_pred[i] = normal_rng(A[j[i]] + (X[i] * B[,sp_obs[i]]), sigma_y);
    log_lik[i] = normal_lpdf(y[i] | A[j[i]] + (X[i] * B[,sp_obs[i]]), sigma_y);
  }
}
