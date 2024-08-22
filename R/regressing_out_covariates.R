run_regress_out_covariates_test <- function(Y_list, x, Z, h = 15L, alpha = 0.1, resid_type = c("response", "deviance", "pearson")[1]) {
  if (!(resid_type %in% c("response", "deviance", "pearson"))) stop("Residual type not recognized.")
  # fit null GLMs and perform precomputation
    precomp_list <- lapply(X = Y_list, FUN = function(y) {
      fit <- MASS::glm.nb(y ~ Z)
      list(stats::resid(fit, type = resid_type))
    })
    result <- run_adaptive_permutation_test(precomp_list, x, h, alpha)

}
