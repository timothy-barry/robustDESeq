#' Run standard NB regression
#'
#' @inheritParams run_robust_nb_regression
#'
#' @return a data frame containing the p-values and rejections
#' @export
run_standard_nb_regression <- function(Y_list, x, Z, side = "two_tailed", alpha = 0.1, method = "MASS") {
  if (!(side %in% c("left", "right", "two_tailed"))) stop("`side` not recognized.")
  if (method == "MASS") {
    p_values <- sapply(X = Y_list, FUN = function(y) {
      fit <- MASS::glm.nb(y ~ x + Z)
      s <- summary(fit)
      z <- coef(s)["x", "z value"]
      compute_gaussian_p_value(z, side)
    })
  } else if (method == "VGAM") {
    p_values <- sapply(X = Y_list, FUN = function(y) {
      fit <- VGAM::vglm(y ~ x + Z, family = VGAM::negbinomial())
      s <- VGAM::summaryvglm(fit)
      z <- s@coef3["x", "z value"]
      compute_gaussian_p_value(z, side)
    })
  } else {
    stop("Method not recognized.")
  }
  rejected <- stats::p.adjust(p_values, method = "BH") < alpha
  get_result_df(p_values, rejected)
}


#' Run regress out residuals test
#'
#' @inheritParams run_robust_nb_regression
#'
#' @return a data frame containing the p-values and rejections
#' @export
run_regress_out_covariates_test <- function(Y_list, x, Z, side = "two_tailed", h = 15L, alpha = 0.1, resid_type = c("response", "deviance", "pearson")[1]) {
  if (!(resid_type %in% c("response", "deviance", "pearson"))) stop("Residual type not recognized.")
  side_code <- get_side_code(side)
  # fit null GLMs and perform precomputation
    precomp_list <- lapply(X = Y_list, FUN = function(y) {
      fit <- MASS::glm.nb(y ~ Z)
      r <- stats::setNames(stats::resid(fit, type = resid_type), NULL)
      list(r)
    })
    result <- run_adaptive_permutation_test(precomp_list, x, side_code, h, alpha, "compute_mean_over_treated_units")
    df <- data.frame(p_value = result$p_values,
                     rejected = result$rejected,
                     hyp_idx = seq_len(m)) |> dplyr::arrange(p_value)
}


#' Run linear model
#'
#' @inheritParams run_robust_nb_regression
#'
#' @return a data frame containing the p-values and rejections
#' @export
run_linear_model <- function(Y_list, x, Z, side = "two_tailed", alpha = 0.1) {
  if (!(side %in% c("left", "right", "two_tailed"))) stop("`side` not recognized.")
  p_values <- sapply(X = Y_list, FUN = function(y) {
    # normalize y
    log_y <- log(y + 1)
    y_norm <- (log_y - mean(log_y))/sd(log_y)
    fit <- stats::lm(y_norm ~ x + Z)
    s <- summary(fit)
    t <- coef(s)["x", "t value"]
    deg_of_freedom <- fit$df.residual
    compute_t_dist_p_value(t, deg_of_freedom, side)
  })
  rejected <- stats::p.adjust(p_values, method = "BH") < alpha
  get_result_df(p_values, rejected)
}
