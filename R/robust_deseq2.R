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
#' res <- run_robust_deseq(dds)
run_robust_deseq <- function(dds, side = "two_tailed", h = 15L, alpha = 0.1, max_iterations = 200000L) {
  # get the side of the test
  side_code <- get_side_code(side)

  # modify the formula
  form <- DESeq2::design(dds)
  terms <- attr(terms(form), "term.labels")
  trt_covariate_name <- terms[length(terms)]
  new_form <- as.formula(paste0("~", paste0(terms[-length(terms)], collapse = "+")))
  DESeq2::design(dds) <- new_form
  trt_vector <- dds[[trt_covariate_name]]

  # convert x to binary variable; extract count matrix
  one_label <- names(which.min(table(trt_vector)))
  x <- as.integer(trt_vector == one_label)
  count_matrix <- SummarizedExperiment::assays(dds)$counts
  rownames(count_matrix) <- colnames(count_matrix) <- NULL

  # run the deseq regression under the null hypothesis
  dds <- DESeq2::DESeq(dds)

  # get mu and thetas
  thetas <- 1/dispersions(dds)
  Z_model <- stats::model.matrix(design(dds), colData(dds))
  beta_mat <- coef(dds)
  size_factors <- sizeFactors(dds)
  rownames(beta_mat) <- colnames(beta_mat) <- colnames(Z_model) <- rownames(Z_model) <- names(size_factors) <- NULL
  mu_mat <- size_factors * t(exp(Z_model %*% t(beta_mat)))

  # perform the precomputation
  precomp_list <- lapply(X = seq_len(nrow(count_matrix)), FUN = function(i) {
    compute_precomputation_pieces_deseq(count_matrix[i,], mu_mat[i,], Z_model, thetas[i])
  })

  # run the permutation test
  print("Running permutations.")
  result <- run_adaptive_permutation_test(precomp_list, x, side_code, h, alpha, max_iterations, "compute_score_stat")
  as.data.frame(result) |> setNames(c("p_value", "rejected", "n_losses", "stop_time"))
}


#' Run robust DESeq (list interface)
#'
#' @inheritParams run_robust_nb_regression
#'
#' @return results data frame
#' @export
#'
#' @examples
#' n <- 100L
#' m <- 500L
#' theta_gt <- 5
#' Z <- MASS::mvrnorm(n = n, mu = c(-0.5, 0.5), Sigma = toeplitz(c(1, 0.5)))
# x <- rbinom(n = n, size = 1, prob = binomial()$linkinv(-1 + as.numeric(Z %*% c(0.8, 0.7))))
#' x <- rbinom(n = n, size = 1, prob = 0.3)
#" colnames(Z) <- c("1", "2")
#' family_object <- MASS::negative.binomial(theta_gt)
#' design_matrix <- cbind(x, Z)
#' colnames(design_matrix) <- c("x", "z1", "z2")
#' null_coefs <- log(c(25, 1.0, 0.8, 1.1))
#' alt_coefs <- log(c(25, 1.5, 0.8, 1.1))
#' under_null <- sample(c(rep(TRUE, 0.9 * m), rep(FALSE, 0.1 * m)))
#' Y_list <- sapply(X = seq_len(m), FUN = function(i) {
#'  generate_glm_data(design_matrix = design_matrix,
#'                    coefficients = if (under_null[i]) null_coefs else alt_coefs,
#'                    family_object = family_object,
#'                    add_intercept = TRUE)
#' }, simplify = FALSE)
#' h <- 15L; alpha <- 0.1; theta <- 30
#'
#' # mass implementation
#' res_robust_nb <- run_robust_nb_regression(Y_list, x, Z)
#' get_result_metrics(res_robust_nb, under_null)
#'
#' # deseq implementation
#' res_deseq <- run_deseq_list_interface(Y_list, x, Z, robust = FALSE)
#' get_result_metrics(res_deseq, under_null)
run_deseq_list_interface <- function(Y_list, x, Z, side = "two_tailed", h = 15L, alpha = 0.1, robust = TRUE, max_iterations = 200000L) {
  Y_mat <- sapply(X = Y_list, FUN = identity) |> t()
  x_fct <- factor(x = x, levels = c(0L, 1L), labels = c("untrt", "trt"))
  covariate_matrix <- data.frame(Z, x_fct)
  colnames(covariate_matrix) <- c(paste0("z", seq_len(ncol(Z))), "x")
  form <- paste0("~", paste0(colnames(covariate_matrix), collapse = "+")) |> as.formula()
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = Y_mat,
                                        colData = covariate_matrix,
                                        design = form)
  # pass to `run_robust_deseq`
  if (robust) {
    res <- run_robust_deseq(dds = dds, side = side, h = h, alpha = alpha, max_iterations = max_iterations)
  } else {
    dds <- DESeq2::DESeq(dds)
    res_orig <- DESeq2::results(dds, independentFiltering = FALSE, cooksCutoff = FALSE)
    rejected <- p.adjust(p = res_orig$pvalue, method = "BH") < alpha
    res <- data.frame(p_value = res_orig$pvalue, rejected = rejected)
  }
  return(res)
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
