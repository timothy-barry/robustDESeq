compute_gaussian_p_value <- function(z, side) {
  switch(side,
         left = pnorm(z, lower.tail = TRUE),
         right = pnorm(z, lower.tail = FALSE),
         two_tailed = 2 * pnorm(abs(z), lower.tail = FALSE))
}

compute_t_dist_p_value <- function(t, deg_of_freedom, side) {
  switch(side,
         left = pt(q = t, df = deg_of_freedom, lower.tail = TRUE),
         right = pt(q = t, df = deg_of_freedom, lower.tail = FALSE),
         two_tailed = 2 * pt(q = abs(t), df = deg_of_freedom, lower.tail = FALSE))
}

get_side_code <- function(side) {
  if (!(side %in% c("left", "right", "two_tailed"))) stop("`side` not recognized.")
  switch(side, left = -1L, right = 1L, two_tailed = 0L)
}

get_result_metrics <- function(result, under_null) {
  rejected <- result$rejected
  n_true_discoveries <- sum(!under_null[rejected])
  n_false_discoveries <- sum(under_null[rejected])
  n_discoveries <- sum(rejected)
  n_alternatives <- sum(!under_null)
  fdp <- if (n_discoveries >= 1L) mean(under_null[rejected]) else 0
  return(list(n_true_discoveries = n_true_discoveries,
              n_false_discoveries = n_false_discoveries,
              n_discoveries = n_discoveries, fdp = fdp, tpr = n_true_discoveries/n_alternatives))
}

get_theta_from_fitted_glm <- function(fit) {
  w <- fit$weights
  mu <- fit$fitted.values
  theta <- mean((w * mu)/(mu - w))
  return(theta)
}


#' @export
make_qq_plot <- function(to_plot) {
  p_qq_trans <- ggplot2::ggplot(data = to_plot, mapping = ggplot2::aes(y = p_value)) +
    sceptre:::stat_qq_points(size = 0.8, ymin = 1e-8) +
    sceptre:::stat_qq_band() +
    ggplot2::scale_x_continuous(trans = sceptre:::revlog_trans(10)) +
    ggplot2::scale_y_continuous(trans = sceptre:::revlog_trans(10)) +
    ggplot2::labs(x = "Expected null p-value", y = "Observed p-value") +
    ggplot2::geom_abline(col = "black")
  p_qq_trans
}

generic_evaluation_function <- function(output, metric) {
  resid_types <- c("response", "deviance", "pearson")
  if (all(names(output) %in% resid_types)) {
    sapply(resid_types, function(resid_type) output[[resid_type]][[metric]])
  } else {
    output[[metric]]
  }
}
