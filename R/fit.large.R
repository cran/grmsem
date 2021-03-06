#' Prefitted model: large data set
#'
#' A large quad-variate data set was simulated assuming an underlying Cholesky model, 
#' with 5000 observations per trait and high polygenicity (20,000 SNPs per genetic factor).
#' Genetic trait variances were set to 0.30, 0.60, 0.60 and 0.70 respectively, and residual variances
#' to 0.70, 0.40, 0.40 and 0.30. The data set is described in full, including genetic and residual
#' covariances, in the vignette. A 4-variate Cholesky model was fitted to the data as described in 
#' the vignette and the output has been saved in fit.large.RData.
#' @format gsem.fit output (a list object)
#' \describe{
#'   \item{model.in}{input parameters}
#'   \item{formula}{fitted formula}
#'   \item{model.fit}{optimisation output}
#'   \item{model.out}{fitted gsem model}
#'   \item{VCOV}{variance covariance matrix}
#'   \item{k}{Number of phenotypes}
#'   \item{n}{Number of observations across all phenotypes}
#'   \item{n.obs}{Number of observations per phenotype}
#'   \item{n.ind}{Number of individuals with at least one phenotype}
#'   \item{model}{gsem model}
#'   \item{con}{constraint}
#'   \item{ph.nms}{phenotype names}   
#' }
"fit.large"
