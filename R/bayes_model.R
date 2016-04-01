#' The stan model used in the Bayesian fit
#'
#' This is the literal model used for fitting in Stan
#' @export

bayesStanModel <- "data {
    int<lower=0> P; // number of elements of the composition
    int<lower=0> N; // the number of subjects sampled
    matrix<lower=0,upper=1>[N,P] C; // the relative abundance of a feature in a sample
    vector[P] nu;  // the mean prior parameter for mu
    cov_matrix[P] Lambda; // the scale prior parameter for mu
    real<lower=1> eta; // the lkj parameter
    real<lower=0> alpha[P]; // the shape parameters of sigma
    real<lower=0> beta[P]; // the scale parameters of sigma
}

parameters {
    vector[P] mu; // the lognormal centrality parameter
    cholesky_factor_corr[P] L; // the cholesky decomposition of the correlation parameter
    vector<lower=0>[P] sigma;  // the standard deviations for the lognormal scale
}

transformed parameters {
  cov_matrix[P] Sigma; // the lognormal scale parameter
  corr_matrix[P] Rho; // the lognormal correlation parameter
  vector[P] ln_mu;    // the lognormal mean
  cov_matrix[P] ln_Sigma; // the lognormal covariance
  vector[P]  R_denom;  //
  corr_matrix[P] ln_Rho; // the lognormal correlation

  Rho <- L * L';
  Sigma <- quad_form_diag(Rho, sigma);

  ln_mu <- exp(mu + 0.5 * diagonal(Sigma));
  ln_Sigma <- tcrossprod(rep_matrix(ln_mu, 1)) .* (exp(Sigma) - 1);
  for (i in 1:P){
    R_denom[i] <- inv_sqrt(exp(Sigma[i, i]) - 1);
  }
  // This method of calculating ln_Rho eliminates ln_mu, so doesn't involve
  // such large numbers
  ln_Rho   <- quad_form_diag(exp(Sigma) - 1, R_denom);
}

model {
    matrix[N,P] alpha_star;
    matrix[P,P] invSigma;
    real sigma_sq_star;
    vector[N] mu_star;
    real inc_1;
    real inc_2;
    vector[N] inc_3i;
    vector[N] inc_4i;

    mu ~ multi_normal(nu, Lambda);
    L ~ lkj_corr_cholesky(eta);
    for(k in 1:P){
      sigma[k] ~ gamma(alpha[k], beta[k]);
    }

    invSigma <- inverse(Sigma);
    sigma_sq_star <- inv(sum(invSigma));

    alpha_star   <- rep_matrix(to_row_vector(mu), N) - log(C);
    mu_star <- rows_dot_product(alpha_star * invSigma, rep_matrix(1, N, P)) * sigma_sq_star;

    inc_1 <- 0.5 * N * log(sigma_sq_star);
    inc_2 <- - 0.5 * N * log_determinant(Sigma);

    increment_log_prob(inc_1 + inc_2);

    for (i in 1:N){
        inc_3i[i] <-  - 0.5 * quad_form(invSigma, alpha_star[i]' );
        inc_4i[i] <-  0.5 * mu_star[i] * mu_star[i] / sigma_sq_star;
                  }
    increment_log_prob( sum(inc_3i + inc_4i) );
}"
