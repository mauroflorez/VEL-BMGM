#' Update Beta_j for Bayesian Logistic Model
#'
#' Updates the beta_j coefficient vector for a single predictor using Bayesian updates.
#'
#' @param X Design matrix of predictors.
#' @param j Index of predictor to update.
#' @param z_y Adjusted response vector.
#' @param Cj_var List or array of GP covariance matrices.
#' @param gamma Inclusion indicators for predictors.
#' @param omega_pg Polya-Gamma weights.
#' @param beta_0 Intercept.
#' @param beta_j_current Current beta_j matrix.
#'
#' @return Updated beta_j vector for predictor j.
#' @export
update_beta_j <- function(X, j, z_y, Cj_var, gamma, omega_pg, beta_0, beta_j_current){
  n <- nrow(X)
  p <- ncol(X)
  Xj <- diag(X[,j])
  if(gamma[j] == 0){
    return(rep(0, n))
  }
  included_vars <- which(gamma == 1 & (1:p != j))
  eta_y <- beta_0 + rowSums(beta_j_current[, included_vars, drop = FALSE] * X[,included_vars, drop = FALSE])
  epsilon_j <- z_y - eta_y
  V_omega_inv <- chol(diag(X[,j]^2 * omega_pg) + solve(Cj_var[,,j]))
  V_omega <- chol2inv(V_omega_inv)
  m_omega <- V_omega %*% (X[,j] * omega_pg * epsilon_j)
  beta_j_new <- mvnfast::rmvn(1, mu = c(m_omega), sigma = V_omega)
  return(beta_j_new)
}

#' Update Intercept Beta_0 for Bayesian Logistic Model
#'
#' @param X Design matrix of predictors.
#' @param z_y Adjusted response vector.
#' @param beta_j_current Current beta_j matrix.
#' @param gamma Inclusion indicators for predictors.
#' @param omega Observation weights.
#' @param tau Precision for prior on beta_0.
#'
#' @return Updated beta_0 scalar.
#' @export
update_beta_0 <- function(X, z_y, beta_j_current, gamma, omega, tau){
  active_vars <- which(gamma == 1)
  eta_y <- rowSums(beta_j_current[, active_vars, drop = FALSE] * X[, active_vars, drop = FALSE])
  residual <- z_y - eta_y
  var_0 <- 1/(sum(omega) + tau)
  mean_0 <- var_0 * sum(omega * residual)
  beta_0_new <- rnorm(1, mean = mean_0, sd = sqrt(var_0))
  return(beta_0_new)
}

#' Update Tau_w Hyperparameter for GP Priors
#'
#' @param lastgammak_all Inclusion indicators (matrix).
#' @param lastw_all GP weights (matrix).
#' @param a_tau_w Shape parameter for tau_w prior.
#' @param b_tau_w Rate parameter for tau_w prior.
#'
#' @return Updated tau_w value.
#' @export
update_tau_w <- function(lastgammak_all, lastw_all, a_tau_w, b_tau_w) {
  num_active <- sum(lastgammak_all)
  sum_w2 <- sum(lastgammak_all * (lastw_all^2))
  tau_w_new <- rgamma(1, shape = a_tau_w + 0.5 * num_active,
                      rate = b_tau_w + 0.5 * sum_w2)
  return(tau_w_new)
}
