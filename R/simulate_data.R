#' Generate GLM data
#'
#' `generate_glm_data()` generates synthetic GLM data.
#'
#' @param design_matrix a design matrix
#' @param coefficients coefficients (possibly including an intercept)
#' @param family_object a family object
#' @param add_intercept logical indicating whether to add an intercept to the design matrix
#' @return a vector of synthetic data generated from the GLM
#' @export
#' @examples
#' ##########################
#' # GENERATE SINGLE DATASET
#' #########################
#' n <- 7000
#' x <- rbinom(n = n, size = 1, prob = 0.25)
#' Z <- MASS::mvrnorm(n = n, mu = c(-0.5, 0.5), Sigma = toeplitz(c(1, 0.5)))
#' colnames(Z) <- c("1", "2")
#' family_object <- MASS::negative.binomial(5)
#' design_matrix <- cbind(x, Z)
#' colnames(design_matrix) <- c("x", "z1", "z2")
#' coefs <- log(c(7, 1.3, 0.8, 1.1))
#' y <- generate_glm_data(
#'   design_matrix = design_matrix,
#'   coefficients = coefs,
#'   family_object = family_object,
#'   add_intercept = TRUE
#' )
#' # MASS with unknown size
#' fit_theta_est_mass <- MASS::glm.nb(formula = y ~ x + Z)
#' summary(fit_theta_est_mass)
#' # VGAM with unknown size
#' fit_theta_est_vgam <- VGAM::vglm(y ~ x + Z, family = VGAM::negbinomial())
#' VGAM::summary(fit_theta_est_vgam)
#' # GLM with known size
#' fit_theta_known <- glm(formula = y ~ x + Z, family = family_object)
#'
#' ############################
#' # GENERATE MULTIPLE DATASETS
#' ############################
#' Y_list <- replicate(n = 1000L,
#'   expr = {generate_glm_data(design_matrix, coefs, family_object, TRUE)},
#'   simplify = FALSE)
#' fits <- lapply(Y_list, function(y) {
#'   fit_theta_est_mass <- MASS::glm.nb(formula = y ~ x + Z)
#' })
generate_glm_data <- function(design_matrix, coefficients, family_object, add_intercept = TRUE) {
  if (!is(design_matrix, "matrix")) stop("`design_matrix` must be an object of class matrix.")
  family_object <- family_object |> augment_family_object()
  if (add_intercept) design_matrix <- cbind(1, design_matrix)
  eta <- as.numeric(design_matrix %*% coefficients)
  mus <- family_object$linkinv(eta)
  y <- family_object$generate_samples(length(mus), mus)
  return(y)
}


augment_family_object <- function(family_object) {
  if (is.null(family_object$augmented) || !family_object$augmented) {
    fam_str <- gsub(x = family_object$family, replacement = "", pattern = "\\((.*)")
    if (fam_str == "Negative Binomial") {
      family_object$theta <- get(x = ".Theta", envir = environment(family_object$variance))
    }
    generate_samples <- switch(EXPR = fam_str,
                               "poisson" = stats::rpois,
                               "binomial" = function(n, mu) stats::rbinom(n = n, size = 1, prob = mu),
                               "Negative Binomial" = function(n, mu) {
                                 MASS::rnegbin(n = n, mu = mu, theta = family_object$theta)
                               })
    family_object$generate_samples <- generate_samples
    family_object$augmented <- TRUE
  }
  return(family_object)
}
