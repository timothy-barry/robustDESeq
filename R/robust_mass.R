#' Run robust NB regression
#'
#' @param Y_list a list of response vectors
#' @param x the treatment vector
#' @param Z the covariate matrix
#' @param side the side of the test, one of "left", "right", or "two_tailed"
#' @param h tuning parameter for the anytime-valid permutation test
#' @param alpha nominal FDR of the BH procedure
#' @param method the method to use to carry out the NB regression under the null model, either "MASS" or "VGAM"
#'
#' @return a data frame containing the rejection set
#' @export
#'
#' @examples
#' ############################
#' # GENERATE MULTIPLE DATASETS
#' ############################
#' set.seed(10)
#' n <- 1000L
#' m <- 500L
#' theta_gt <- 5
#' Z <- MASS::mvrnorm(n = n, mu = c(-0.5, 0.5), Sigma = toeplitz(c(1, 0.5)))
#' # x <- rbinom(n = n, size = 1, prob = binomial()$linkinv(-1 + as.numeric(Z %*% c(0.8, 0.7))))
#' x <- rbinom(n = n, size = 1, prob = 0.3)
#' colnames(Z) <- c("1", "2")
#' family_object <- MASS::negative.binomial(theta_gt)
#' design_matrix <- cbind(x, Z)
#' colnames(design_matrix) <- c("x", "z1", "z2")
#' null_coefs <- log(c(7, 1.0, 0.8, 1.1))
#' alt_coefs <- log(c(7, 1.1, 0.8, 1.1))
#' under_null <- sample(c(rep(TRUE, 0.9 * m), rep(FALSE, 0.1 * m)))
#' Y_list <- sapply(X = seq_len(m), FUN = function(i) {
#'   generate_glm_data(design_matrix = design_matrix,
#'   coefficients = if (under_null[i]) null_coefs else alt_coefs,
#'   family_object = family_object,
#'   add_intercept = TRUE)
#' }, simplify = FALSE)
#' h <- 15L; alpha <- 0.1; theta <- 30
#'
#' ########################
#' # STANDARD NB REGRESSION
#' ########################
#' standard_res <- run_standard_nb_regression(Y_list = Y_list, x = x, Z = Z, side = "right")
#' get_result_metrics(standard_res, under_null)
#'
#' ######################
#' # ROBUST NB REGRESSION
#' ######################
#' robust_res <- run_robust_nb_regression(Y_list = Y_list, x = x, Z = Z, side = "right")
#' get_result_metrics(robust_res, under_null)
#'
#' ###########################
#' # REGRESSING OUT COVARIATES
#' ###########################
#' residual_res <- run_regress_out_covariates_test(Y_list = Y_list, x = x, Z = Z, side = "right")
#' get_result_metrics(residual_res, under_null)
run_robust_nb_regression <- function(Y_list, x, Z, side = "two_tailed", h = 15L, alpha = 0.1, method = "MASS", size_factors = NULL, theta = NULL, max_iterations = 200000L) {
  if (2/choose(length(x), sum(x)) > 5e-4) {
    warning("Your sample size may be too small for the permutation test to make any significant hits. Consider using a different method (e.g., DESeq2) or increasing your sample size.")
  }
  Z_model <- cbind(1, Z)
  colnames(Z_model) <- rownames(Z_model) <- NULL
  side_code <- get_side_code(side)
  my_offsets <- if (!is.null(size_factors)) log(size_factors) else NULL

  # fit null GLMs and perform precomputation
  precomp_list <- lapply(X = Y_list, FUN = function(y) {
    # try to fit the NB glm and compute the p-value; otherwise, fit Pois GLM
    if (is.null(theta)) {
      fit <- tryCatch({
        if (is.null(size_factors)) {
          suppressWarnings(MASS::glm.nb(y ~ Z))
        } else {
          suppressWarnings(MASS::glm.nb(y ~ Z + offset(my_offsets)))
        }
      }, error = function(e) {
        fit <- stats::glm(y ~ Z, family = poisson(), offset = my_offsets)
        fit$theta <- 10; fit
      })
      theta <- fit$theta
    } else {
      fit <- stats::glm(y ~ Z, family = MASS::negative.binomial(theta), offset = my_offsets)
    }
    compute_precomputation_pieces_v2(y, Z_model, theta, fit)
  })

  result <- run_adaptive_permutation_test(precomp_list, x, side_code, h, alpha, max_iterations, "compute_score_stat")
  as.data.frame(result) |> setNames(c("p_value", "rejected", "n_losses", "stop_time"))
}


compute_precomputation_pieces <- function(y, Z_model, coefs, theta) {
  mu <- as.numeric(exp(Z_model %*% coefs))
  denom <- 1 + mu / theta
  w <- mu / denom
  a <- (y - mu) / denom
  wZ <- w * Z_model
  Zt_wZ <- t(Z_model) %*% wZ
  P_decomp <- eigen(Zt_wZ, symmetric = TRUE)
  U <- P_decomp$vectors
  Lambda_minus_half <- 1 / sqrt(P_decomp$values)
  D <- (Lambda_minus_half * t(U)) %*% t(wZ)
  D_list <- apply(D, 1L, function(row) row, simplify = FALSE)
  out <- list(a = a, w = w, D_list = D_list)
  return(out)
}


compute_precomputation_pieces_v2 <- function(y, Z_model, theta, fit) {
  mu <- stats::setNames(fit$fitted.values, NULL)
  w <- stats::setNames(fit$weights, NULL)
  qr <- fit$qr
  a <- (y - mu) / (1 + mu / theta)
  wZ <- w * Z_model
  R <- qr.R(qr)
  R_t_inv <- backsolve(r = R, upper.tri = TRUE, x = diag(nrow(R)), transpose = TRUE)
  D <- R_t_inv %*% t(wZ)
  D_list <- apply(D, 1L, function(row) row, simplify = FALSE)
  out <- list(a = a, w = w, D_list = D_list)
  return(out)
}
