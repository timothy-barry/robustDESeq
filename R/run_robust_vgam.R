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
#' ######################
#' # GENERATE ONE DATASET
#' ######################
#' library(VGAM)
#' set.seed(10)
#' n <- 5000L
#' m <- 100L
#' theta_gt <- 5
#' x <- rbinom(n = n, size = 1, prob = 0.3)
#' null_coefs <- log(c(7, 1.0, 0.8, 1.1))
#' alt_coefs <- log(c(7, 1.1, 0.8, 1.1))
#' under_null <- sample(c(rep(TRUE, 0.9 * m), rep(FALSE, 0.1 * m)))
#' dat_list <- lapply(seq_len(100), function(i) {
#'   Z <- MASS::mvrnorm(n = n, mu = c(-0.5, 0.5), Sigma = toeplitz(c(1, 0.5)))
#'   colnames(Z) <- c("1", "2")
#'   family_object <- MASS::negative.binomial(theta_gt)
#'   design_matrix <- cbind(x, Z)
#'   colnames(design_matrix) <- c("x", "z1", "z2")
#'   y <- generate_glm_data(design_matrix = design_matrix,
#'       coefficients = if (under_null[i]) null_coefs else alt_coefs,
#'       family_object = family_object,
#'       add_intercept = TRUE
#'     )
#'     list(Z = Z, y = y)
#' })
#' y_list <- lapply(X = dat_list, FUN = function(l) l[["y"]])
#' Z_list <- lapply(X = dat_list, FUN = function(l) l[["Z"]])
#'
#' # standard VGAM
#' robust_res <- run_robust_vgam(y_list, x, Z_list)
#' get_result_metrics(robust_res, under_null)
run_robust_vgam <- function(y_list, x, Z_list, side = "two_tailed", h = 15L, alpha = 0.1, max_iterations = 200000L, custom_permutation_list = list()) { # Z does not contain an intercept
  side_code <- get_side_code(side)
  precomp_list <- lapply(X = seq_along(y_list), FUN = function(i) {
    y <- y_list[[i]]
    Z <- Z_list[[i]]
    Z_model <- cbind(1, Z)
    colnames(Z_model) <- rownames(Z_model) <- NULL
    # fit the GLM to esetimate theta
    suppressWarnings(fit_init <- VGAM::vglm(y ~ Z, family = VGAM::negbinomial()))
    theta <- stats::setNames(exp(fit_init@coefficients[2]), NULL)
    # rerun the GLM with fixed theta
    suppressWarnings(fit <- VGAM::vglm(y ~ Z, family = VGAM::negbinomial.size(size = theta)))
    # get the precomputation pieces
    R <- fit@R
    mu <- as.numeric(fit@fitted.values)
    w <- mu/(1 + mu/theta)
    compute_precomputation_pieces_vgam(y, Z_model, theta, mu, w, R)
  })
  result <- run_adaptive_permutation_test_v2(precomp_list, x, side_code, h, alpha, max_iterations, "compute_score_stat", custom_permutation_list)
  as.data.frame(result) |> setNames(c("p_value", "rejected", "n_losses", "stop_time"))
}


compute_precomputation_pieces_vgam <- function(y, Z_model, theta, mu, w, R) {
  a <- (y - mu) / (1 + mu / theta)
  wZ <- w * Z_model
  R_t_inv <- backsolve(r = R, upper.tri = TRUE, x = diag(nrow(R)), transpose = TRUE)
  D <- R_t_inv %*% t(wZ)
  D_list <- apply(D, 1L, function(row) row, simplify = FALSE)
  out <- list(a = a, w = w, D_list = D_list)
  return(out)
}

#' @export
run_robust_poisson <- function(y_list, x, Z_list, theta = 1000, side = "two_tailed", h = 15L, alpha = 0.1, max_iterations = 200000L, custom_permutation_list = list()) {
  side_code <- get_side_code(side)
  precomp_list <- lapply(X = seq_along(y_list), FUN = function(i) {
    y <- y_list[[i]]
    Z <- Z_list[[i]]
    Z_model <- cbind(1, Z)
    colnames(Z_model) <- rownames(Z_model) <- NULL
    # fit the GLM to esetimate thet
    suppressWarnings(fit <- glm(y ~ Z, family = poisson()))
    # compute the precomputation pieces
    compute_precomputation_pieces_mass(y, Z_model, theta, fit)
  })
  result <- run_adaptive_permutation_test_v2(precomp_list, x, side_code, h, alpha, max_iterations, "compute_score_stat", custom_permutation_list)
  as.data.frame(result) |> setNames(c("p_value", "rejected", "n_losses", "stop_time"))
}
