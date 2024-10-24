#' Run standard NB regression
#'
#' @inheritParams run_robust_nb_regression
#'
#' @return a data frame containing the p-values and rejections
#' @export
run_standard_nb_regression <- function(Y_list, x, Z, side = "two_tailed", alpha = 0.1, method = "MASS", size_factors = NULL, theta = NULL) {
  if (!(side %in% c("left", "right", "two_tailed"))) stop("`side` not recognized.")
  my_offsets <- if (!is.null(size_factors)) log(size_factors) else NULL
  p_values <- sapply(X = Y_list, FUN = function(y) {
    # try to fit the glm and compute the p-value; otherwise, return 1
    p <- tryCatch({
      if (is.null(theta)) {
        if (is.null(size_factors)) {
          suppressWarnings(fit <- MASS::glm.nb(y ~ x + Z))
        } else {
          suppressWarnings(fit <- MASS::glm.nb(y ~ x + Z + offset(my_offsets)))
        }
      } else {
        fit <- stats::glm(y ~ x + Z, family = MASS::negative.binomial(theta), offset = my_offsets)
      }
      fit$family <- stats::poisson()
      s <- summary(fit)
      z <- coef(s)["x", "z value"]
      compute_gaussian_p_value(z, side)
    }, error = function(e) 1)
    return(p)
  })
  rejected <- stats::p.adjust(p_values, method = "BH") < alpha
  data.frame(p_value = p_values, rejected = rejected)
}


#' Run regress out residuals test
#'
#' @inheritParams run_robust_nb_regression
#'
#' @return a data frame containing the p-values and rejections
#' @export
run_regress_out_covariates_test <- function(Y_list, x, Z, side = "two_tailed", h = 15L, alpha = 0.1, resid_type = c("response", "deviance", "pearson", "all")[1], theta = NULL, max_iterations = 50000L) {
  resid_types <- c("response", "deviance", "pearson")
  if (!(resid_type %in% c(resid_types, "all"))) stop("Residual type not recognized.")
  side_code <- get_side_code(side)
  # fit null GLMs and perform precomputation
  precomp_list <- lapply(X = Y_list, FUN = function(y) {
    if (is.null(theta)) {
      fit <- tryCatch({
        suppressWarnings(MASS::glm.nb(y ~ Z))
      }, error = function(e) {
        stats::glm(y ~ Z, family = poisson())
      })
    } else {
      fit <- stats::glm(y ~ Z, family = MASS::negative.binomial(theta))
    }
    if (resid_type != "all") {
      list(stats::setNames(stats::resid(fit, type = resid_type), NULL))
    } else {
      lapply(resid_types, function(curr_resid_type) {
        stats::resid(fit, type = curr_resid_type)
      }) |> stats::setNames(resid_types)
    }
  })
  if (resid_type == "all") {
    out <- lapply(X = resid_types, FUN = function(curr_resid_type) {
      curr_precomp_list <- lapply(precomp_list, FUN = function(item) item[curr_resid_type])
      curr_result <- run_adaptive_permutation_test(curr_precomp_list, x, side_code, h,
                                                   alpha, max_iterations, "compute_mean_over_treated_units") |>
        as.data.frame() |> setNames(c("p_value", "rejected", "n_losses", "stop_time")) |> dplyr::mutate(resid_type = curr_resid_type)

    }) |> data.table::rbindlist() |> as.data.frame()
  } else {
    out <- run_adaptive_permutation_test(precomp_list, x, side_code, h, alpha, max_iterations, "compute_mean_over_treated_units") |>
      as.data.frame() |> setNames(c("p_value", "rejected", "n_losses", "stop_time"))
  }
  return(out)
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
  data.frame(p_value = p_values, rejected = rejected)
}
