test_that("wor sample correcntess", {
  set.seed(10)
  n <- 20L
  n_trt <- 10L
  n_samples <- 100000L
  samp <- generate_wor_sample_test(n = n, n_trt = n_trt, n_samples = n_samples)
  # verify each sample consists of unique elements
  testthat::expect_true(all(sapply(samp, function(l) length(unique(l)) == length(l))))
  samp_unlist <- unlist(samp)
  testthat::expect_equal(min(samp_unlist), 0L)
  testthat::expect_equal(max(samp_unlist), n - 1L)
  # verify that the proportion of each element is about 1/n
  props <- prop.table(table(samp_unlist))
  all(props > 0.049 & props < 0.051)
})
