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
run_robust_nb_regression <- function(Y_list, x, Z, side = "two_tailed", h = 15L, alpha = 0.1, method = "MASS", theta = NULL, adaptive_permutation_test = TRUE) {
  Z_model <- cbind(1, Z)
  colnames(Z_model) <- rownames(Z_model) <- NULL
  side_code <- get_side_code(side)

  # fit null GLMs and perform precomputation
  if (method == "MASS") {
    precomp_list <- lapply(X = Y_list, FUN = function(y) {
      if (is.null(theta)) {
        fit <- MASS::glm.nb(y ~ Z)
        theta <- fit$theta
      } else {
        fit <- stats::glm(y ~ Z, family = MASS::negative.binomial(theta))
      }
      coefs <- fit$coefficients
      compute_precomputation_pieces(y, Z_model, coefs, theta)
    })
  } else if (method == "VGAM") {
    precomp_list <- lapply(X = Y_list, FUN = function(y) {
      fit <- VGAM::vglm(y ~ Z, family = VGAM::negbinomial())
      coefs <- stats::setNames(fit@coefficients[-2], NULL)
      theta <- stats::setNames(exp(fit@coefficients[2]), NULL)
      compute_precomputation_pieces(y, Z_model, coefs, theta)
    })
  } else {
    stop("Method not recognized.")
  }

  # run the permutation test
  if (adaptive_permutation_test) {
    result <- run_adaptive_permutation_test(precomp_list, x, side_code, h, alpha, "compute_score_stat")
    p_values <- result$p_values; rejected <- result$rejected
  } else {
    p_values <- run_permutation_test(precomp_list, x, side_code, round(10 * m/alpha), "compute_score_stat")
    rejected <- stats::p.adjust(p_values, method = "BH") < alpha
  }
  get_result_df(p_values = p_values, rejected = rejected)
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
