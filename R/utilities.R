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

get_result_df <- function(p_values, rejected) {
  tibble::tibble(p_value = p_values,
                 rejected = rejected,
                 hyp_idx = seq_len(length(p_values))) |> dplyr::arrange(p_value)
}

get_side_code <- function(side) {
  if (!(side %in% c("left", "right", "two_tailed"))) stop("`side` not recognized.")
  switch(side, left = -1L, right = 1L, two_tailed = 0L)
}

get_result_metrics <- function(result, under_null) {
  n_true_discoveries <- sum(!under_null[result |> dplyr::filter(rejected) |> dplyr::pull(hyp_idx)])
  fdp <- mean(under_null[result |> dplyr::filter(rejected) |> dplyr::pull(hyp_idx)])
  return(list(n_true_discoveries = n_true_discoveries, fdp = fdp))
}

get_theta_from_fitted_glm <- function(fit) {
  w <- fit$weights
  mu <- fit$fitted.values
  theta <- mean((w * mu)/(mu - w))
  return(theta)
}
