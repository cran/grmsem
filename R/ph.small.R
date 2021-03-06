#' Phenotype data: small data set
#'
#' A small tri-variate data set was simulated assuming an underlying Cholesky model, 
#' with 100 observations per trait and low polygenicity (150 SNPs per genetic factor).
#' The data set is described in full in the vignette. Genetic and residual
#' covariances are assumed to be influenced by three independent genetic factors 
#' (A1, A2 and A3, based on a 150 SNPs each) and three independent residual factors (E1, E2 and E3), 
#' respectively.
#' @format A data frame with 100 rows and 3 columns
#' \describe{
#'   \item{Y1}{trait one}
#'   \item{Y2}{trait two}
#'   \item{Y3}{trait three}
#' }
"ph.small"
