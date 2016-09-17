#' The stan model used in the Bayesian fit
#'
#' This is the literal model used for fitting in Stan
#'
#' @return The BAnOCC model as a string to be compiled with
#'     \code{rstan::stan_model}
#' 
#' @export

banocc_model <- "data {
    int<lower=0> P; // number of elements of the composition
    int<lower=0> N; // the number of subjects sampled
    matrix<lower=0,upper=1>[N,P] C; // matrix of N comp. samples, P features
    vector[P] n;  // the mean prior parameter for mu
    cov_matrix[P] L; // the scale prior parameter for mu
    real<lower=1> eta; // the lkj parameter
    real<lower=0> a[P]; // the shape parameters of sigma
    real<lower=0> b[P]; // the scale parameters of sigma
}

parameters {
    vector[P] m; // the lognormal centrality parameter
    cholesky_factor_corr[P] WChol; // the cholesky decomposition of the correlation parameter
    vector<lower=0>[P] s;  // the standard deviations for the lognormal scale
}

transformed parameters {
  cov_matrix[P] S; // the lognormal scale parameter
  corr_matrix[P] W; // the lognormal correlation parameter

  W = WChol * WChol';
  S = quad_form_diag(W, s);
}

model {
    matrix[N,P] alpha_star;
    matrix[P,P] invS;
    real s_sq_star;
    vector[N] m_star;
    real inc_1;
    real inc_2;
    vector[N] inc_3i;
    vector[N] inc_4i;

    m ~ multi_normal(n, L);
    WChol ~ lkj_corr_cholesky(eta);
    for(k in 1:P){
      s[k] ~ gamma(a[k], b[k]);
    }

    invS = inverse(S);
    s_sq_star = inv(sum(invS));

    alpha_star   = rep_matrix(to_row_vector(m), N) - log(C);
    m_star = rows_dot_product(alpha_star * invS, rep_matrix(1, N, P)) * s_sq_star;

    inc_1 = 0.5 * N * log(s_sq_star);
    inc_2 = - 0.5 * N * log_determinant(S);

    target += inc_1 + inc_2;

    for (i in 1:N){
        inc_3i[i] =  - 0.5 * quad_form(invS, alpha_star[i]' );
        inc_4i[i] =  0.5 * m_star[i] * m_star[i] / s_sq_star;
                  }
    target += sum(inc_3i + inc_4i);
}"
