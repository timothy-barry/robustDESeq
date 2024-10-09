#' Run Mann-Whitney test (asymptotic)
#'
#' @inheritParams run_robust_nb_regression
#' @param implementation which implementation to use for the MW test, either "custom" (for the in-house implementation) or "r" (for the standard R implementation)
#'
#' @return a data frame containing the results
#' @export
run_mann_whitney_test_asymptotic <- function(Y_list, x, Z, side = "two_tailed", implementation = "custom", alpha = 0.1) {
  if (!(implementation %in% c("custom", "r"))) stop("`implementation` not recognized.")
  if (implementation == "r") {
    alternative <- switch(side, two_tailed = "two.sided", right = "greater", left = "less")
  }
  p_values <- if (implementation == "custom") {
    sapply(X = Y_list, FUN = function(y) {
      s_trt <- y[x == 1L]
      s_cntrl <- y[x == 0L]
      combined <- c(s_trt, s_cntrl)
      r <- rank(combined)
      n_s_trt <- length(s_trt)
      n_s_cntrl <- length(s_cntrl)
      statistic <- sum(r[seq_along(s_trt)]) - n_s_trt * (n_s_trt + 1)/2
      n_ties <- table(r)
      statistic_zero <- statistic - n_s_trt * n_s_cntrl/2
      sigma <- sqrt((n_s_trt * n_s_cntrl/12) * ((n_s_trt + n_s_cntrl + 1) - sum(n_ties^3 - n_ties)/((n_s_trt + n_s_cntrl) * (n_s_trt + n_s_cntrl - 1))))
      correction <- switch(side, two_tailed = sign(statistic_zero) * 0.5, right = 0.5, left = -0.5)
      z <- (statistic_zero - correction)/sigma
      compute_gaussian_p_value(z, side = side)
    })
  } else {
    sapply(X = Y_list, FUN = function(y) {
      s_trt <- y[x == 1L]
      s_cntrl <- y[x == 0L]
      fit <- stats::wilcox.test(x = s_trt, y = s_cntrl, alternative = alternative, correct = TRUE, exact = FALSE)
      fit$p.value
    })
  }
  rejected <- stats::p.adjust(p_values, method = "BH") < alpha
  data.frame(p_value = p_values, rejected = rejected)
}


#' Run Mann-Whitney test (finite-sample)
#'
#' @inheritParams run_robust_nb_regression
#'
#' @return  a data frame containing the results
#' @export
run_mann_whitney_test_permutations <- function(Y_list, x, Z, side = "two_tailed", h = 15L, B = NULL, alpha = 0.1, adaptive_permutation_test = TRUE) {
  if (2/choose(length(x), sum(x)) > 5e-4) {
    warning("Your sample size may be too small for the permutation test to make any significant hits. Consider using a different method (e.g., DESeq2) or increasing your sample size.")
  }
  side_code <- get_side_code(side)
  # iterate over hypotheses, performing precomputation
  precomp_list <- lapply(X = Y_list, FUN = function(y) {
    r <- rank(y)
    n_s1 <- as.double(sum(x == 1L))
    n_s2 <- as.double(length(x) - n_s1)
    n_ties <- as.numeric(table(r))
    sigma <- sqrt((n_s1 * n_s2/12) * ((n_s1 + n_s2 + 1) - sum(n_ties^3 - n_ties)/((n_s1 + n_s2) * (n_s1 + n_s2 - 1))))
    list(r = r, sigma = sigma, side_code = side_code)
  })
  # run the permutation test
  if (adaptive_permutation_test) {
    result <- run_adaptive_permutation_test(precomp_list, x, side_code, h, alpha, "compute_mw_test_statistic")
    out <- as.data.frame(result) |> setNames(c("p_value", "rejected", "n_losses", "stop_time"))
  } else {
    if (is.null(B)) B <- round(10 * length(Y_list)/alpha)
    p_values <- run_permutation_test(precomp_list, x, side_code, B, "compute_mw_test_statistic")
    rejected <- stats::p.adjust(p_values, method = "BH") < alpha
    out <- data.frame(p_value = p_values, rejected = rejected)
  }
  return(out)
}
