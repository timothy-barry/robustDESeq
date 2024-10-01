if (FALSE) {
  library(DESeq2)
  library(airway)
  data("airway")
  deseq_object <- DESeqDataSet(se = airway, design = formula(~cell))[1:100,]
  x <- as.integer(colData(deseq_object)$dex == "trt")
  res_orig <- DESeqDataSet(se = airway, design = formula(~cell + dex))[1:100,] |> DESeq() |> results()

  run_robust_deseq <- function(deseq_object, side = "two_tailed", h = 15L, alpha = 0.1) {
    side_code <- get_side_code(side)
    # run DESeq analysis under null hypothesis
    deseq_object <- deseq_object |>
      estimateSizeFactors() |>
      estimateDispersions() |>
      nbinomWaldTest()
    # extract pieces for perm test
    gene_all_zero <- rowData(deseq_object)$allZero
    expression_matrix <- assays(deseq_object)[["counts"]][!gene_all_zero,]
    mu_matrix <- assays(deseq_object)[["mu"]][!gene_all_zero,]
    theta <- 1/(dispersions(deseq_object)[!gene_all_zero])
    Z_model <- model.matrix(design(deseq_object), colData(deseq_object))
    precomp_list <- lapply(X = seq_len(nrow(expression_matrix)), FUN = function(i) {
      compute_precomputation_pieces_deseq(expression_matrix[i,], Z_model, mu_matrix[i,], theta[i])
    })
    result <- run_adaptive_permutation_test(precomp_list, x, side_code, h, alpha, "compute_score_stat")
  }

  compute_precomputation_pieces_deseq <- function(y, Z_model, mu, theta) {
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
}
