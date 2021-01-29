#' grmsem bivariate heritability estimation function.
#'
#' This function estimates the bivariate heritability.
#' @details The \code{grmsem.biher} function estimates the bivariate heritability 
#' (DS, Cholesky, IP and IPC models) from the observed phenotype data and a \code{grmsem.var} object. 
#' All standard errors are derived with the Delta method.
#' @param grmsem.var.out A grmsem.var object with unstandardised parameters (factor loadings). Default NULL.
#' @param ph Phenotype file as R dataframe (columns: >=2 phenotypes, rows: ni individuals in the same order as G). No default.
#' @keywords grmsem
#' @export
#' @return \code{grmsem.biher} returns a list object consisting of the following matrices:
#' \item{VPO}{observed phenotypic variance/covariance matrix}
#' \item{VA}{estimated genetic variance}
#' \item{BIHER}{estimated bivariate heritability (off-diagonals): VA / VPO}
#' \item{BIHER.se}{standard error of estimated bivariate heritability}
#' \item{BIHER.Z}{Z (wald) of estimated bivariate heritability}
#' \item{BIHER.p}{p (Wald) of estimated bivariate heritability}
#' @examples
#' #(runtime should be less than one minute)
#' \donttest{
#' out <- grmsem.fit(ph.small, G.small, LogL = TRUE, estSE = TRUE)
#' var.out <- grmsem.var(out)
#' grmsem.biher(ph.small, var.out)}
grmsem.biher <- function(ph, grmsem.var.out = NULL){
  local
  ################################################################################
  #Main body
  ################################################################################
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Full length
  if (class(grmsem.var.out) != "grmsem.var") {
    stop('The function requires a grmsem.var object.')
  }
  if (is.null(grmsem.var.out$VA)) {
    stop('There is no grmsem.var.out information on genetic variances available.')
  }
  VPO<- stats::cov(ph,use="pairwise.complete.obs") #Observed phenotypic correlations
  colnames(VPO) <- seq_len(ncol(VPO))
  if (ncol(VPO) < 2){
    stop('The estimation of bivariate heritabilities requires at least 2 phenotypes.')
  }
  VA <- grmsem.var.out$VA
  colnames(VA) <- seq_len(ncol(VPO))
  rownames(VA) <- rownames(VPO)
  if (any(is.na(VA))) {
    stop('Genetic variance/covariance matrix A contains missing values.')
  }
  VA.se <- grmsem.var.out$VA.se
  if (any(is.na(VA.se))) {
    stop('SE matrix of genetic variance/covariance matrix A contains missing values.')
  }
  if (ncol(VA) != ncol(VPO)) {
    m.msg <- paste("The number of phenotypes in the observed phenotypic covariance matrix VPO (", ncol(VPO), ") does not match those in VA (", ncol(VA), ").", sep = "")
    stop(m.msg)
  }
  BIHER <- VA / VPO
  colnames(BIHER) <- seq_len(ncol(VPO))
  rownames(BIHER) <- rownames(VPO)
  #Approximation of BIHER SE estimation, 
  #as observed variance covariance matrix (VPO) is estimated with high precision
  BIHER.se <- VA.se / VPO
  BIHER.se <- as.matrix(BIHER.se)
  colnames(BIHER.se) <- seq_len(ncol(VPO))
  rownames(BIHER.se) <- rownames(VPO)
  BIHER.Z <- BIHER / BIHER.se
  colnames(BIHER.Z) <- seq_len(ncol(VPO))
  rownames(BIHER.Z) <- rownames(VPO)
  BIHER.p <- 2*stats::pnorm(-abs(BIHER.Z))
  colnames(BIHER.p) <- seq_len(ncol(VPO))
  rownames(BIHER.p) <- rownames(VPO)
  out <- (list(VPO = VPO, VA = VA, BIHER = BIHER, BIHER.se = BIHER.se, BIHER.Z = BIHER.Z, BIHER.p = BIHER.p))
  class(out) <- ("grmsem.biher")
  return(out)
}
