#' Run robust DESeq
#'
#' @param dds a DESeq object
#' @inheritParams run_robust_nb_regression
#'
#' @return the result data frame
#' @export
#'
#' @examples
#' library(DESeq2)
#' dds <- readRDS("/Users/tib163/research_offsite/projects/camp/data/arrayed_crispr/deseq_object.rds")
#' design(dds) <- formula(~donor_id + stim_status + grna)
#' row_sums <- rowSums(assays(dds)$counts)
#' highly_expressed <- row_sums >= 5
#' dds <- dds[highly_expressed,]
#' res <- run_robust_deseq(dds)
run_robust_deseq <- function(dds, side = "two_tailed", h = 15L, alpha = 0.1, size_factors = NULL, dispersion_estimation = "local", max_iterations = 200000L) {
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

  # run deseq
  # 1. size factors
  if (is.null(size_factors)) {
    dds <- DESeq2::estimateSizeFactors(dds)
  } else {
    DESeq2::sizeFactors(dds) <- size_factors
  }
  # 2. dispersion
  if (dispersion_estimation == "raw") {
    dds <- DESeq2::estimateDispersions(dds, fitType = "mean")
    DESeq2::dispersions(dds) <- SummarizedExperiment::mcols(dds)$dispGeneEst
  } else {
    dds <- DESeq2::estimateDispersions(dds, fitType = dispersion_estimation)
  }
  # 3. NB GLM
  dds <- DESeq2::nbinomWaldTest(object = dds)

  # get mu, thetas, and Z_model
  thetas <- 1/DESeq2::dispersions(dds)
  Z_model <- stats::model.matrix(DESeq2::design(dds), SummarizedExperiment::colData(dds))
  mu_mat <- SummarizedExperiment::assays(dds)$mu

  # alternately, manually compute mu_mat
  # beta_mat <- coef(dds) * log(2)
  # size_factors <- sizeFactors(dds)
  # mu_mat_manual <- t(exp(Z_model %*% t(beta_mat) + log(size_factors)))

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
#' @export
#' @examples
#' n <- 100L
#' m <- 500L
#' theta_gt <- 5
#' offsets <- log(rpois(n = n, lambda = 50))
#' Z <- MASS::mvrnorm(n = n, mu = c(-0.5, 0.5), Sigma = toeplitz(c(1, 0.5)))
#' x <- rbinom(n = n, size = 1, prob = 0.3)
#' colnames(Z) <- c("1", "2")
#' family_object <- MASS::negative.binomial(theta_gt)
#' design_matrix <- cbind(x, Z)
#' colnames(design_matrix) <- c("x", "z1", "z2")
#' null_coefs <- log(c(5, 1.0, 0.8, 1.1))
#' alt_coefs <- log(c(5, 1.5, 0.8, 1.1))
#' under_null <- sample(c(rep(TRUE, 0.9 * m), rep(FALSE, 0.1 * m)))
#' Y_list <- sapply(X = seq_len(m), FUN = function(i) {
#'  generate_glm_data(design_matrix = design_matrix,
#'                    coefficients = if (under_null[i]) null_coefs else alt_coefs,
#'                    family_object = family_object,
#'                    add_intercept = TRUE,
#'                    offsets = offsets)
#' }, simplify = FALSE)
#'
#' res <- run_robust_deseq_list_interface(Y_list, x, Z)
run_robust_deseq_list_interface <- function(Y_list, x, Z, side = "two_tailed", h = 15L, alpha = 0.1, size_factors = NULL, dispersion_estimation = "local", max_iterations = 200000L) {
  dds <- make_deseq_object(Y_list, x, Z)
  out <- run_robust_deseq(dds = dds, side = side, h = h,
                          alpha = alpha, size_factors = size_factors,
                          dispersion_estimation = dispersion_estimation,
                          max_iterations = max_iterations)
  return(out)
}


#' Run standard DESeq (list interface)
#' @export
run_standard_deseq_list_interface <- function(Y_list, x, Z, side = "two_tailed", alpha = 0.1, dispersion_estimation = "local", size_factors = NULL) {
  dds <- make_deseq_object(Y_list, x, Z)
  if (is.null(size_factors)) {
    dds <- DESeq2::estimateSizeFactors(dds)
  } else {
    DESeq2::sizeFactors(dds) <- size_factors
  }
  dds <- DESeq2::estimateDispersions(dds, fitType = dispersion_estimation)
  dds <- DESeq2::nbinomWaldTest(dds)
  res <- DESeq2::results(dds, independentFiltering = FALSE, cooksCutoff = FALSE)
  rejected <- p.adjust(p = res$pvalue, method = "BH") < alpha
  out <- data.frame(p_value = res$pvalue, rejected = rejected)
  return(out)
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


make_deseq_object <- function(Y_list, x, Z) {
  Y_mat <- sapply(X = Y_list, FUN = identity) |> t()
  x_fct <- factor(x = x, levels = c(0L, 1L), labels = c("untrt", "trt"))
  covariate_matrix <- data.frame(Z, x_fct)
  colnames(covariate_matrix) <- c(paste0("z", seq_len(ncol(Z))), "x")
  form <- paste0("~", paste0(colnames(covariate_matrix), collapse = "+")) |> as.formula()
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = Y_mat,
                                        colData = covariate_matrix,
                                        design = form)
  return(dds)
}
