#' Title
#'
#' @param Y_list
#' @param x
#' @param Z
#' @param side
#' @param h
#' @param alpha
#' @param max_iterations
#'
#' @return
#' @export
#'
#' @examples
#' library(DESeq2)
#' dds <- readRDS("/Users/tib163/research_offsite/projects/camp/data/arrayed_crispr/deseq_object.rds")
#' design(dds) <- formula(~donor_id + stim_status + grna)
#' row_sums <- rowSums(assays(dds)$counts)
#' highly_expressed <- row_sums >= 6
#' dds <- dds[highly_expressed,]
run_robust_deseq <- function(dds, side = "two_tailed", h = 15L, alpha = 0.1, max_iterations = 200000L) {
  # get the side of the test
  side_code <- get_side_code(side)

  # modify the formula
  form <- design(dds)
  terms <- attr(terms(form), "term.labels")
  trt_covariate_name <- terms[length(terms)]
  new_form <- as.formula(paste0("~", paste0(terms[-length(terms)], collapse = "+")))
  design(dds) <- new_form
  trt_vector <- dds[[trt_covariate_name]]

  # convert x to binary variable; extract count matrix
  zero_label <- names(which.max(table(trt_vector)))
  x <- as.integer(trt_vector == zero_label)
  count_matrix <- assays(dds)$counts
  rownames(count_matrix) <- colnames(count_matrix) <- NULL

  # run the deseq regression under the null hypothesis
  dds <- DESeq2::DESeq(dds)

  # get mu and thetas
  thetas <- 1/dispersions(dds)
  Z_model <- stats::model.matrix(design(dds), colData(dds))
  beta_mat <- coef(dds)
  rownames(beta_mat) <- colnames(beta_mat) <- colnames(Z_model) <- rownames(Z_model) <- names(size_factors) <- NULL
  size_factors <- sizeFactors(dds)
  mu_mat <- size_factors * t(exp(Z_model %*% t(beta_mat)))

  # perform the precomputation
  precomp_list <- lapply(X = seq_len(nrow(count_matrix)), FUN = function(i) {
    compute_precomputation_pieces_deseq(count_matrix[i,], mu_mat[i,], Z_model, theta)
  })

  # run the permutation test
  result <- run_adaptive_permutation_test(precomp_list, x, side_code, h, alpha, max_iterations, "compute_score_stat")
  df <- as.data.frame(result) |> setNames(c("p_value", "rejected", "n_losses", "stop_time"))
}


compute_precomputation_pieces_deseq <- function(y, mu, Z_model, theta) {
  denom <- 1 + mu / theta
  w <- mu / denom
  a <- (y - mu) / denom
  wZ <- w * Z_model
  Zt_wZ <- t(Z_model) %*% wZ
  P_decomp <- eigen(Zt_wZ, symmetric = TRUE)
  U <- P_decomp$vectors
  Lambda_minus_half <- 1 / sqrt(P_decomp$values)
  D <- (Lambda_minus_half * t(U)) %*% t(wZ)
  D_list <- apply(D, 1L, function(row) row, simplify = FALSE)
  out <- list(a = a, w = w, D_list = D_list)
  return(out)
}
