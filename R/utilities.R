compute_gaussian_p_value <- function(z, side) {
  switch(side,
         left = pnorm(z, lower.tail = TRUE),
         right = pnorm(z, lower.tail = FALSE),
         two_tailed = 2 * pnorm(abs(z), lower.tail = FALSE))
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
