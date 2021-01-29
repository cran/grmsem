#' Phenotype data: large data set
#'
#' A large quad-variate data set was simulated assuming an underlying Cholesky model, 
#' with 5000 observations per trait and high polygenicity (5,000 SNPs per genetic factor).
#' The data set is described in full in the vignette. Genetic and residual
#' covariances ar assumed to be influenced by four independent genetic factors
#' (A1, A2, A3 and A4, based on 5000 SNPs each) and four independent residual factors (E1, E2, E3 and E4), 
#' respectively.
#' @format A data frame with 5000 rows and 4 columns
#' \describe{
#'   \item{Y1}{trait one}
#'   \item{Y2}{trait two}
#'   \item{Y3}{trait three}
#'   \item{Y4}{trait four}
#' }
"ph.large"
