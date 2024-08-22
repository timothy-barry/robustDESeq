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
