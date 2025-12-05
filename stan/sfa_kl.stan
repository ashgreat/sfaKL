// Stan Model for Kumbhakar-Lai (2021) Multi-Output Multi-Input SFA
// with Input- and Output-Specific Inefficiency
//
// Reference: Kumbhakar, S.C. & Lai, H.P. (2021). A multi-output multi-input
// stochastic frontier system with input- and output-specific inefficiency.
// Economics Letters.
//
// The composite error follows a Closed Skew Normal (CSN) distribution.

functions {
  // Log probability of observing epsilon given CSN parameters
  // epsilon ~ CSN(0, Theta, Psi, 0, Delta)
  // f(epsilon) = phi(epsilon; 0, Theta) * Phi(Psi * epsilon; 0, Delta) / Phi(0; 0, Sigma)
  real csn_lpdf(vector epsilon, matrix Theta, matrix Psi, matrix Delta, matrix Sigma) {
    int n_eq = rows(epsilon);
    int n_ineff = cols(Psi);
    
    // Term 1: log phi(epsilon; 0, Theta)
    real term1 = multi_normal_lpdf(epsilon | rep_vector(0.0, n_eq), Theta);
    
    // Term 2: log Phi(Psi' * epsilon; 0, Delta)
    // Psi' * epsilon gives the argument to the CDF
    vector[n_ineff] cdf_arg = Psi' * epsilon;
    
    // For multivariate normal CDF, we use the approximation via multi_normal_cdf
    // Stan doesn't have multi_normal_cdf directly, so we use a latent variable approach
    // This is handled in the model block via data augmentation
    
    // We return only term1 here; the CDF terms are handled via latent variables
    return term1;
  }
}

data {
  int<lower=1> N;                    // Number of observations
  int<lower=2> J;                    // Number of inputs
  int<lower=1> M;                    // Number of outputs
  
  // Shares (N x n_eq where n_eq = J-1 + M)
  matrix[N, J-1] S;                  // Input cost shares S_2, ..., S_J (negative of paper's -S)
  matrix[N, M] R;                    // Output revenue shares R_1, ..., R_M
  
  // Log price ratios
  matrix[N, J-1] ln_w;               // ln(w_j / w_1) for j=2,...,J
  matrix[N, M] ln_p;                 // ln(p_m / w_1) for m=1,...,M
}

transformed data {
  int n_eq = J - 1 + M;              // Number of equations
  int n_ineff = J + M;               // Number of inefficiency terms (mu_1,...,mu_J, delta_1,...,delta_M)
  
  // Combine shares into single matrix for convenience
  matrix[N, n_eq] Y;
  for (i in 1:N) {
    for (j in 1:(J-1)) {
      Y[i, j] = -S[i, j];            // -S_j as in paper eq (12)
    }
    for (m in 1:M) {
      Y[i, J-1+m] = R[i, m];         // R_m as in paper eq (13)
    }
  }
}

parameters {
  // Profit function parameters
  vector[J-1] alpha;                 // alpha_2, ..., alpha_J (alpha_1 derived)
  vector[M] beta;                    // beta_1, ..., beta_M
  
  // Second-order terms (using Cholesky for symmetry)
  // A: (J-1) x (J-1) symmetric matrix
  // B: M x M symmetric matrix
  // Gamma: (J-1) x M matrix (no symmetry constraint)
  vector[(J-1)*(J-2)/2 + (J-1)] A_lower_tri;  // Lower triangular of A
  vector[M*(M-1)/2 + M] B_lower_tri;          // Lower triangular of B
  matrix[J-1, M] Gamma;
  
  // Covariance parameters (Cholesky factors for positive definiteness)
  cholesky_factor_cov[J] L_Sigma_mu;          // Cholesky of Sigma_mu (input inefficiency)
  cholesky_factor_cov[M] L_Sigma_delta;       // Cholesky of Sigma_delta (output inefficiency)
  cholesky_factor_cov[n_eq] L_Omega;          // Cholesky of Omega (noise)
  
  // Latent inefficiencies (for each observation)
  // mu_i <= 0, delta_i >= 0
  // We parameterize as: mu_raw ~ N(0, Sigma_mu), mu = -|mu_raw|
  //                     delta_raw ~ N(0, Sigma_delta), delta = |delta_raw|
  matrix[N, J] mu_raw;               // Raw input inefficiencies (before truncation)
  matrix[N, M] delta_raw;            // Raw output inefficiencies (before truncation)
}

transformed parameters {
  // Reconstruct symmetric matrices A and B from lower triangular
  matrix[J-1, J-1] A;
  matrix[M, M] B;
  
  {
    int idx = 1;
    for (i in 1:(J-1)) {
      for (j in 1:i) {
        A[i, j] = A_lower_tri[idx];
        A[j, i] = A_lower_tri[idx];
        idx += 1;
      }
    }
  }
  
  {
    int idx = 1;
    for (i in 1:M) {
      for (j in 1:i) {
        B[i, j] = B_lower_tri[idx];
        B[j, i] = B_lower_tri[idx];
        idx += 1;
      }
    }
  }
  
  // Covariance matrices from Cholesky factors
  matrix[J, J] Sigma_mu = multiply_lower_tri_self_transpose(L_Sigma_mu);
  matrix[M, M] Sigma_delta = multiply_lower_tri_self_transpose(L_Sigma_delta);
  matrix[n_eq, n_eq] Omega = multiply_lower_tri_self_transpose(L_Omega);
  
  // Block diagonal Sigma = diag(Sigma_mu, Sigma_delta)
  matrix[n_ineff, n_ineff] Sigma = rep_matrix(0.0, n_ineff, n_ineff);
  Sigma[1:J, 1:J] = Sigma_mu;
  Sigma[(J+1):n_ineff, (J+1):n_ineff] = Sigma_delta;
  
  // Truncated inefficiencies
  matrix[N, J] mu;                   // mu <= 0
  matrix[N, M] delta;                // delta >= 0
  
  for (i in 1:N) {
    for (j in 1:J) mu[i, j] = -abs(mu_raw[i, j]);
    for (m in 1:M) delta[i, m] = abs(delta_raw[i, m]);
  }
  
  // Construct H matrix (n_eq x n_ineff) as in paper eq (16)
  // H = [-H1, -Gamma; -H2, -B]
  // H1 = -A*(-l_{J-1}, I_{J-1}) + Gamma*(l_M, 0)
  // H2 = -Gamma'*(-l_{J-1}, I_{J-1}) + B*(l_M, 0)
  matrix[n_eq, n_ineff] H;
  
  {
    // Build H1: (J-1) x J
    matrix[J-1, J] H1;
    for (j in 1:(J-1)) {
      H1[j, 1] = 0;  // Initialize
      for (k in 1:(J-1)) {
        if (k == 1) {
          H1[j, 1] = A[j, k];  // -A * (-1) = A
        }
        H1[j, k+1] = -A[j, k];  // -A * I
      }
      // Add Gamma * (l_M, 0) contribution
      for (m in 1:M) {
        H1[j, 1] += Gamma[j, m];
      }
    }
    
    // Build H2: M x J
    matrix[M, J] H2;
    for (m in 1:M) {
      H2[m, 1] = 0;
      for (k in 1:(J-1)) {
        if (k == 1) {
          H2[m, 1] = Gamma[k, m];  // -Gamma' * (-1)
        }
        H2[m, k+1] = -Gamma[k, m];  // -Gamma' * I
      }
      // Add B * (l_M, 0) contribution
      for (n in 1:M) {
        H2[m, 1] += B[m, n];
      }
    }
    
    // Assemble H
    for (j in 1:(J-1)) {
      for (k in 1:J) H[j, k] = -H1[j, k];
      for (m in 1:M) H[j, J+m] = -Gamma[j, m];
    }
    for (m in 1:M) {
      for (k in 1:J) H[J-1+m, k] = -H2[m, k];
      for (n in 1:M) H[J-1+m, J+n] = -B[m, n];
    }
  }
  
  // Theta = Omega + H * Sigma * H' (total covariance of epsilon)
  matrix[n_eq, n_eq] Theta = Omega + H * Sigma * H';
}

model {
  // === PRIORS ===
  
  // Profit function parameters (weakly informative)
  alpha ~ normal(0, 1);
  beta ~ normal(0, 1);
  A_lower_tri ~ normal(0, 0.5);
  B_lower_tri ~ normal(0, 0.5);
  to_vector(Gamma) ~ normal(0, 0.5);
  
  // Covariance Cholesky factors (LKJ prior for correlation, half-normal for scales)
  // This is the KEY for identification - informative priors on variances
  L_Sigma_mu ~ lkj_corr_cholesky(2.0);
  L_Sigma_delta ~ lkj_corr_cholesky(2.0);
  L_Omega ~ lkj_corr_cholesky(2.0);
  
  // Diagonal elements (standard deviations) - informative to help identification
  for (j in 1:J) L_Sigma_mu[j, j] ~ normal(0.15, 0.1);
  for (m in 1:M) L_Sigma_delta[m, m] ~ normal(0.15, 0.1);
  for (k in 1:n_eq) L_Omega[k, k] ~ normal(0.2, 0.1);
  
  // Latent inefficiencies (half-normal via folded normal)
  for (i in 1:N) {
    mu_raw[i] ~ multi_normal_cholesky(rep_vector(0.0, J), L_Sigma_mu);
    delta_raw[i] ~ multi_normal_cholesky(rep_vector(0.0, M), L_Sigma_delta);
  }
  
  // === LIKELIHOOD ===
  
  // For each observation, compute the composite error and its likelihood
  for (i in 1:N) {
    // Deterministic parts of share equations
    vector[J-1] det_S;
    vector[M] det_R;
    
    for (j in 1:(J-1)) {
      det_S[j] = alpha[j];
      for (k in 1:(J-1)) det_S[j] += A[j, k] * ln_w[i, k];
      for (m in 1:M) det_S[j] += Gamma[j, m] * ln_p[i, m];
    }
    
    for (m in 1:M) {
      det_R[m] = beta[m];
      for (n in 1:M) det_R[m] += B[m, n] * ln_p[i, n];
      for (j in 1:(J-1)) det_R[m] += Gamma[j, m] * ln_w[i, j];
    }
    
    // Inefficiency contributions (u and omega from paper)
    vector[J-1] u_i;
    vector[M] omega_i;
    
    {
      // mu_diff = mu_{-1} - l_{J-1} * mu_1
      vector[J-1] mu_diff;
      for (j in 1:(J-1)) mu_diff[j] = mu[i, j+1] - mu[i, 1];
      
      // delta_diff = delta - l_M * mu_1
      vector[M] delta_diff;
      for (m in 1:M) delta_diff[m] = delta[i, m] - mu[i, 1];
      
      // u_j = -(A * mu_diff + Gamma * delta_diff)
      u_i = -(A * mu_diff + Gamma * delta_diff);
      
      // omega_m = -(B * delta_diff + Gamma' * mu_diff)
      omega_i = -(B * delta_diff + Gamma' * mu_diff);
    }
    
    // Composite errors (should be close to noise v)
    vector[n_eq] epsilon_i;
    for (j in 1:(J-1)) epsilon_i[j] = Y[i, j] - det_S[j] - u_i[j];
    for (m in 1:M) epsilon_i[J-1+m] = Y[i, J-1+m] - det_R[m] - omega_i[m];
    
    // epsilon_i | mu, delta ~ N(0, Omega)
    epsilon_i ~ multi_normal_cholesky(rep_vector(0.0, n_eq), L_Omega);
  }
}

generated quantities {
  // Technical efficiencies for each input and output
  matrix[N, J] TE_input;             // exp(-|mu|) in [0, 1]
  matrix[N, M] TE_output;            // 1/exp(|delta|) in [0, 1]
  
  for (i in 1:N) {
    for (j in 1:J) TE_input[i, j] = exp(mu[i, j]);  // mu <= 0, so exp(mu) <= 1
    for (m in 1:M) TE_output[i, m] = exp(-delta[i, m]);  // delta >= 0
  }
  
  // Log-likelihood for model comparison
  vector[N] log_lik;
  for (i in 1:N) {
    vector[J-1] det_S;
    vector[M] det_R;
    
    for (j in 1:(J-1)) {
      det_S[j] = alpha[j];
      for (k in 1:(J-1)) det_S[j] += A[j, k] * ln_w[i, k];
      for (m in 1:M) det_S[j] += Gamma[j, m] * ln_p[i, m];
    }
    
    for (m in 1:M) {
      det_R[m] = beta[m];
      for (n in 1:M) det_R[m] += B[m, n] * ln_p[i, n];
      for (j in 1:(J-1)) det_R[m] += Gamma[j, m] * ln_w[i, j];
    }
    
    vector[J-1] u_i;
    vector[M] omega_i;
    {
      vector[J-1] mu_diff;
      for (j in 1:(J-1)) mu_diff[j] = mu[i, j+1] - mu[i, 1];
      vector[M] delta_diff;
      for (m in 1:M) delta_diff[m] = delta[i, m] - mu[i, 1];
      u_i = -(A * mu_diff + Gamma * delta_diff);
      omega_i = -(B * delta_diff + Gamma' * mu_diff);
    }
    
    vector[n_eq] epsilon_i;
    for (j in 1:(J-1)) epsilon_i[j] = Y[i, j] - det_S[j] - u_i[j];
    for (m in 1:M) epsilon_i[J-1+m] = Y[i, J-1+m] - det_R[m] - omega_i[m];
    
    log_lik[i] = multi_normal_cholesky_lpdf(epsilon_i | rep_vector(0.0, n_eq), L_Omega);
  }
}
