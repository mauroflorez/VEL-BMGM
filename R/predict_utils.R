#' Predict with Bayesian Mixed Graphical Model
#'
#' Generates predictions for new data from a fitted VEL.BMGM model.
#'
#' @param fit Model fit object from bmgm_gp().
#' @param X_new New predictor matrix (for test data).
#' @param Z Covariate matrix (includes both train and test).
#' @param beta_j_est Posterior mean or estimate of beta_j.
#' @param PPI_Z Posterior probabilities for kernel components.
#' @param mcmc_samples Number of posterior samples to use (default 1000).
#' @param test Indices of testing data
#' @param train Indices of training data
#' @return Vector of predicted probabilities (for binary Y).
#' @export
predict_jbmgm <- function(fit, X_new, test, train, Z, beta_j_est, PPI_Z, mcmc_samples = 1000){
  n <- nrow(Z)
  K <- ncol(Z)
  p <- ncol(X_new)
  n_test <- nrow(X_new)
  n_train <- n - n_test
  nburn <- fit$nburn
  nsample <- fit$nsample
  X_centered <- scale(X_new, center = TRUE, scale = FALSE)

  Z_norm <- apply(Z, 2, function(x) scales::rescale(x, to = c(0, 1)))
  Kernels_Z_train <- array(0, dim = c(n_train, n_train, K))
  Kernels_Z_test_train <- array(0, dim = c(n_test, n_train, K))

  for(k in 1:K){
    dists_test_train <- as.matrix(proxy::dist(Z_norm[test, k, drop = FALSE],
                                              Z_norm[train, k, drop = FALSE],
                                              method = "euclidean"))^2
    Kernels_Z_train[,,k] <- kernel_se(Z_norm[train, k], lengthscale = 0.5, sigma_f = 1)
    Kernels_Z_test_train[,,k] <- exp(-0.5 * dists_test_train / 0.5^2)
  }

  calc_Cj <- function(w, r){
    CZZstar <- Reduce('+', lapply(1:K, function(x) w[x]^2 * Kernels_Z_test_train[,,x]))
    CZZ <- Reduce('+', lapply(1:K, function(x) w[x]^2 * Kernels_Z_train[,,x])) + 1/r*diag(n_train)
    CZZ_inv <- solve(CZZ)
    return(CZZstar %*% CZZ_inv)
  }

  sample_est <- sample(1:nsample, mcmc_samples)
  beta_j_star <- array(0, dim = c(n_test, p, length(sample_est)))
  eta_i <- matrix(0, nrow = length(sample_est), ncol = n_test)

  for(m in 1:length(sample_est)){
    for(j in 1:p){
      if(fit$post_gamma[sample_est[m] + nburn,j] == 1){
        ind <- (PPI_Z[j,] > 0.5)*1
        w_t <- fit$post_w[j,,sample_est[m]+nburn]
        r <- fit$post_r[sample_est[m]+nburn,j]
        beta_j_star[,j,m] <- calc_Cj(w_t*ind, r)%*%beta_j_est[,j]
      } else{
        beta_j_star[,j,m] <- rep(0, n_test)
      }
    }
  }

  eta_i <- matrix(0, nrow = length(sample_est), ncol = n_test)
  for (m in 1:length(sample_est)) {
    beta_0_m <- fit$post_beta0[sample_est[m] + nburn]
    eta_i[m, ] <- beta_0_m + rowSums(X_centered * beta_j_star[,, m])
  }

  Y_pred_prob <- LaplacesDemon::invlogit(eta_i)
  Y_pred_prob_mean <- apply(Y_pred_prob, 2, mean)
  y_pred_our <- ifelse(Y_pred_prob_mean > 0.5, 1, 0)

  return(Y_pred_prob_mean)
}
