#' The stan model used in the Bayesian fit
#'
#' This is the literal model used for fitting in Stan
#'
#' @return The BAnOCC model as a string to be compiled with
#'     \code{rstan::stan_model}
#' @examples
#' data(compositions_null)
#' \dontrun{
#'   compiled_banocc_model <- rstan::stan_model(model_code = banocc_model)
#' }
#' 
#' @export

banocc_model <- "data {
  int<lower=0> P; // number of elements of the composition
  int<lower=0> N; // the number of subjects sampled
  matrix<lower=0,upper=1>[N,P] C; // matrix of N comp. samples, P features
  vector[P] n;  // the mean prior parameter for m
  cov_matrix[P] L; // the scale prior parameter for m
  real<lower=0> a; // the shape prior parameter for lambda
  real<lower=0> b; // the rate prior parameter for lambda
}
transformed data{
  matrix[P,P] L_chol;
  matrix[P, N] log_C_t;
  matrix[N, P] one_mat;
  int<lower=0> P_tri;

  L_chol = cholesky_decompose(L);
  log_C_t = log(C)';
  P_tri = P * (P - 1) / 2;
  one_mat = rep_matrix(1, N, P);
}
parameters {
  vector[P] m_raw;
  cov_matrix[P] O; // the log-basis precision matrix
  real<lower=0> lambda;
}
transformed parameters{
  vector[P] m; // the lognormal centrality parameter

  m = n + L_chol * m_raw;
}
model {
  matrix[P,N] alpha_star;
  real s_sq_star;
  vector[N] m_star;
  vector[P] O_diag; // the diagonal values of O
  vector[P_tri] O_tri; // the upper-triangle values of O
  vector[4] lik; // the likelihood value has four parts

  m_raw ~ normal(0, 1); // implies: m ~ multi_normal(n, L)
  lambda ~ gamma(a, b);

  O_diag = diagonal(O);
  O_diag ~ exponential(lambda/2);

  for (k in 1:(P - 1)){
    vector[P - k] O_tri_row;
    O_tri_row = sub_col(O, k + 1, k, P - k);
    O_tri_row ~ double_exponential(0, lambda);
  }

  lik[1] = 0.5 * N * log_determinant(O);

  s_sq_star = inv(sum(O));
  lik[2] = 0.5 * N * log(s_sq_star);

  alpha_star   = rep_matrix(m, N) - log_C_t;
  lik[3] = -0.5 * trace_quad_form(O, alpha_star);

  m_star = rows_dot_product(alpha_star' * O, one_mat);
  // This is a more efficient way of calculating
  //  (1^T * O * alpha_i) * s_sq_star)^2 / s_sq_star
  lik[4] = 0.5 * dot_self(m_star) * s_sq_star;

  target += sum(lik);
}"
