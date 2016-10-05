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
    vector[P] n;  // the mean prior parameter for mu
    cov_matrix[P] L; // the scale prior parameter for mu
    real<lower=0> lambda; // the glasso prior parameter
}

parameters {
    vector[P] m; // the lognormal centrality parameter
    corr_matrix[P] W; // The log-basis correlation matrix
    vector<lower=0>[P] s; // the log-basis standard deviations
}

transformed parameters {
  cov_matrix[P] O; // the log-basis precision matrix
  cov_matrix[P] S; // the log-basis covariance matrix

  S = quad_form_diag(W, s);
  O = inverse(S);
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
    for(k in 1:P){
      O[k,k] ~ exponential(lambda/2);
      for (i in (k+1):P){
        O[i, k] ~ double_exponential(0, lambda);
      }
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
