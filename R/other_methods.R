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


#' Title
#'
#' @param Y_list
#' @param x
#' @param Z
#' @param h
#' @param alpha
#' @param resid_type
#'
#' @return
#' @export
#'
#' @examples
run_regress_out_covariates_test <- function(Y_list, x, Z, h = 15L, alpha = 0.1, resid_type = c("response", "deviance", "pearson")[1]) {
  if (!(resid_type %in% c("response", "deviance", "pearson"))) stop("Residual type not recognized.")
  # fit null GLMs and perform precomputation
    precomp_list <- lapply(X = Y_list, FUN = function(y) {
      fit <- MASS::glm.nb(y ~ Z)
      r <- stats::setNames(stats::resid(fit, type = resid_type), NULL)
      list(r)
    })
    result <- run_adaptive_permutation_test(precomp_list, x, h, alpha, "compute_mean_over_treated_units")
    df <- data.frame(p_value = result$p_values,
                     rejected = result$rejected,
                     hyp_idx = seq_len(m)) |> dplyr::arrange(p_value)
}


#' Run Mann-Whitney test
#'
#' @param x binary vector of treatment indicators
#' @param y observations for the second sample
#'
#' @return the z-score
#' @export
#'
#' @examples
#' # x <- rbinom(n = 500, size = 1, prob = 0.5)
#' # y <- MASS::rnegbin(n = n, mu = 50, theta = 10)
#' n <- 1000
#' # s <- rnorm(n = n, mean = 50)
#' s <- MASS::rnegbin(n = n, mu = 40, theta = 10)
#' x <- s[1:500]
#' y <- s[501:1000]
#' run_mann_whitney_test_simple(x, y)
#' wilcox.test(x, y)
run_mann_whitney_test_simple <- function(s1, s2) {
  combined <- c(x, y)
  r <- rank(combined)
  n.x <- length(x)
  n.y <- length(y)
  STATISTIC <- sum(r[seq_along(x)]) - n.x * (n.x + 1)/2
  TIES <- (length(r) != length(unique(r)))
  NTIES <- table(r)
  z <- STATISTIC - n.x * n.y/2
  SIGMA <- sqrt((n.x * n.y/12) * ((n.x + n.y + 1) - sum(NTIES^3 - NTIES)/((n.x + n.y) * (n.x + n.y - 1))))
  CORRECTION <- switch(alternative, two.sided = sign(z) * 0.5, greater = 0.5, less = -0.5)
  z <- (z - CORRECTION)/SIGMA
}
