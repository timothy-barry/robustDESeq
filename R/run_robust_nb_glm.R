#' Run robust NB regression
#'
#' @param Y_list a list of response vectors
#' @param x the treatment vector
#' @param Z the covariate matrix
#'
#' @return a data frame containing the rejection set
#' @export
#'
#' @examples
#' ############################
#' # GENERATE MULTIPLE DATASETS
#' ############################
#' n <- 100
#' m <- 100L
#' x <- rbinom(n = n, size = 1, prob = 0.25)
#' Z <- MASS::mvrnorm(n = n, mu = c(-0.5, 0.5), Sigma = toeplitz(c(1, 0.5)))
#' colnames(Z) <- c("1", "2")
#' family_object <- MASS::negative.binomial(5)
#' design_matrix <- cbind(x, Z)
#' colnames(design_matrix) <- c("x", "z1", "z2")
#' null_coefs <- log(c(7, 1.0, 0.8, 1.1))
#' alt_coefs <- log(c(7, 1.3, 0.8, 1.1))
#' under_null <- sample(size = m, c(TRUE, FALSE), replace = TRUE, prob = c(0.9, 0.1))
#' Y_list <- sapply(X = seq_len(m), FUN = function(i) {
#'   generate_glm_data(design_matrix = design_matrix,
#'   coefficients = if (under_null[i]) null_coefs else alt_coefs,
#'   family_object = family_object,
#'   add_intercept = TRUE)
#' }, simplify = FALSE)
#' h <- 15L; alpha <- 0.1; theta <- NULL
#'
#' ########################
#' # STANDARD NB REGRESSION
#' ########################
#' standard_res <- run_standard_nb_regression(Y_list = Y_list, x = x, Z = Z)
#' mean(under_null[standard_res |> dplyr::filter(rejected) |> dplyr::pull(hyp_idx)])
#'
#' ######################
#' # ROBUST NB REGRESSION
#' ######################
#' robust_res <- run_robust_nb_regression(Y_list = Y_list, x = x, Z = Z)
#' mean(under_null[robust_res |> dplyr::filter(rejected) |> dplyr::pull(hyp_idx)])
run_robust_nb_regression <- function(Y_list, x, Z, h = 15L, alpha = 0.1, theta = NULL) {
  # perform the precomputation on each y vector
  precomp_list <- lapply(X = Y_list, FUN = function(y) {
    # fit the GLM
    fit <- MASS::glm.nb(y ~ Z)
    # compute the precomputation pieces
    compute_precomputation_pieces(fit)
  })
  # run the permutation test
  result <- run_adaptive_permutation_test(precomp_list, x, h, alpha)
  # result <- run_adaptive_permutation_test(precomp_list, x, h, alpha)
  #df <- data.frame(p_value = result$p_values,
  #                 rejected = result$rejected,
  #                 hyp_idx = seq_len(m)) |> dplyr::arrange(p_value)
}


compute_precomputation_pieces <- function(fit) {
  # extract components of fit
  y <- stats::setNames(fit$y, NULL)
  Z <- stats::model.matrix(fit)
  rownames(Z) <- NULL
  coefs <- fit$coefficients
  mu <- as.numeric(exp(Z %*% coefs))
  theta <- fit$theta
  denom <- 1 + mu / theta
  w <- mu / denom
  a <- (y - mu) / denom
  wZ <- w * Z
  Zt_wZ <- t(Z) %*% wZ
  P_decomp <- eigen(Zt_wZ, symmetric = TRUE)
  U <- P_decomp$vectors
  Lambda_minus_half <- 1 / sqrt(P_decomp$values)
  D <- (Lambda_minus_half * t(U)) %*% t(wZ)
  out <- list(a = a, w = w, D = D)
  return(out)
}


run_standard_nb_regression <- function(Y_list, x, Z, alpha = 0.1) {
  p_values <- sapply(X = Y_list, FUN = function(y) {
    fit <- MASS::glm.nb(y ~ x + Z)
    s <- summary(fit)
    coef(s)["x", "Pr(>|z|)"]
  })
  rejected <- stats::p.adjust(p_values, method = "BH") < 0.1
  data.frame(p_value = p_values,
             rejected = rejected,
             hyp_idx = seq_len(m)) |> dplyr::arrange(p_value)
}
