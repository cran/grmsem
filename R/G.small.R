#' Symmetric GRM data: small data set
#'
#' A genetic relationship matrix for a small tri-variate data set was simulated assuming an underlying 
#' Cholesky model, with 100 observations per trait and low polygenicity (150 SNPs per genetic factor).
#' Genetic trait variances were set to 0.30, 0.60, 0.60 and 0.70 respectively, and residual variances
#' to 0.70, 0.40, 0.40 and 0.30. The data set is described in full, including genetic and residual
#' covariances, in the vignette. The traits are influenced by three independent genetic factors 
#' (A1, A2 and A3, based on a 150 SNPs each) and three independent residual factors (E1, E2 and E3).
#' @format Symmetric GRM data frame with 100 rows and 100 columns (observations)
"G.small"
