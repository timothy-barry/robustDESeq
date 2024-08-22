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
#' n <- 1000L
#' m <- 500L
#' theta <- 10
#' Z <- MASS::mvrnorm(n = n, mu = c(-0.5, 0.5), Sigma = toeplitz(c(1, 0.5)))
#' # x <- rbinom(n = n, size = 1, prob = binomial()$linkinv(-1 + as.numeric(Z %*% c(0.8, 0.7))))
#' x <- rbinom(n = n, size = 1, prob = 0.3)
#' colnames(Z) <- c("1", "2")
#' family_object <- MASS::negative.binomial(theta)
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
#' h <- 15L; alpha <- 0.1; theta <- NULL
#'
#' ########################
#' # STANDARD NB REGRESSION
#' ########################
#' standard_res <- run_standard_nb_regression(Y_list = Y_list, x = x, Z = Z)
#' n_true_discoveries_standard <- sum(!under_null[standard_res |> dplyr::filter(rejected) |> dplyr::pull(hyp_idx)])
#' fdp_standard <- mean(under_null[standard_res |> dplyr::filter(rejected) |> dplyr::pull(hyp_idx)])
#'
#' ######################
#' # ROBUST NB REGRESSION
#' ######################
#' robust_res <- run_robust_nb_regression(Y_list = Y_list, x = x, Z = Z)
#' n_true_discoveries_robust <- sum(!under_null[robust_res |> dplyr::filter(rejected) |> dplyr::pull(hyp_idx)])
#' fdp_robust <- mean(under_null[robust_res |> dplyr::filter(rejected) |> dplyr::pull(hyp_idx)])
run_robust_nb_regression <- function(Y_list, x, Z, h = 15L, alpha = 0.1, method = "MASS") {
  Z_model <- cbind(1, Z)
  colnames(Z_model) <- rownames(Z_model) <- NULL

  # fit null GLMs and perform precomputation
  if (method == "MASS") {
    precomp_list <- lapply(X = Y_list, FUN = function(y) {
      fit <- MASS::glm.nb(y ~ Z)
      coefs <- fit$coefficients
      theta <- fit$theta
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
  result <- run_adaptive_permutation_test(precomp_list, x, h, alpha, "compute_score_stat")
  df <- data.frame(p_value = result$p_values,
                   rejected = result$rejected,
                   hyp_idx = seq_len(m)) |> dplyr::arrange(p_value)
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


run_standard_nb_regression <- function(Y_list, x, Z, alpha = 0.1, method = "MASS") {
  if (method == "MASS") {
    p_values <- sapply(X = Y_list, FUN = function(y) {
      fit <- MASS::glm.nb(y ~ x + Z)
      s <- summary(fit)
      coef(s)["x", "Pr(>|z|)"]
    })
  } else if (method == "VGAM") {
    p_values <- sapply(X = Y_list, FUN = function(y) {
      fit <- VGAM::vglm(y ~ x + Z, family = VGAM::negbinomial())
      s <- VGAM::summaryvglm(fit)
      s@coef3["x", "Pr(>|z|)"]
    })
  } else {
    stop("Method not recognized.")
  }
  rejected <- stats::p.adjust(p_values, method = "BH") < 0.1
  data.frame(p_value = p_values,
             rejected = rejected,
             hyp_idx = seq_len(m)) |> dplyr::arrange(p_value)
}
