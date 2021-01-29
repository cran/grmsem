#' grmsem standardised parameters function
#'
#' This function estimates standardised parameters for a grmsem.fit object.
#' @details \code{grmsem.stpar} standardises grmsem.fit estimates so that derived or estimated A and E variance components will add up to phenotypic unit variance. 
#' The SEs of standardised parameters are derived using the Jacobian matrix. 
#' @param grmsem.out grmsem.fit object as provided by the grmsem.fit function. Default NULL.
#' @keywords grmsem
#' @export
#' @return \code{grmsem.stpar} returns a list object consisting of:
#' \item{model.in}{list of input parameters}
#' \item{stand.model.out}{dataframe of fitted grmsem model with standardised parameters and SEs}
#' \item{stVCOV}{standardised variance/covariance matrix}
#' \item{k}{number of phenotypes}
#' \item{n}{total number of observations across all phenotypes}
#' \item{model}{type of grmsem model}
#' \item{ph.nms}{vector of phenotype names}
#' 
#' \code{model.in} list of input parameters:
#' \item{part}{a - genetic, e - residual parameters}
#' \item{label}{parameter label}
#' \item{value}{starting values}
#' \item{freepar}{free parameters}
#' 
#' \code{stand.model.out} data.frame of fitted grmsem model with standardised parameters and SEs:
#' \item{label}{parameter label}
#' \item{estimates}{standardised estimated parameters}
#' \item{se}{SE}
#' \item{Z}{Z (Wald)}
#' \item{p}{p (Wald)}
#' @examples
#' #(runtime should be less than one minute)
#' \donttest{
#' out <- grmsem.fit(ph.small, G.small, LogL = TRUE, estSE = TRUE)
#' stout <- grmsem.stpar(out)
#' print(stout)}
grmsem.stpar <- function(grmsem.out = NULL){
  local
  #################################################################################
  #Subfunctions
  #################################################################################
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Generate standardised coefficients
  grmsem.ds.standcf <- function(estimates) {
    #A part
    VA <- diag(k)
    VA[lower.tri(VA, diag = TRUE)] <- estimates[a.sel]
    VA[upper.tri(VA)] <- t(VA[lower.tri(VA)])
    #E part
    VE <- diag(k)
    VE[lower.tri(VE, diag = TRUE)] <- estimates[e.sel]
    VE[upper.tri(VE)] <- t(VE[lower.tri(VE)])
    VP <- VA + VE
    iSD <- try(solve((diag(k) * VP)), silent = TRUE)
    #Standardised coefficients
    stVA <- iSD %*% VA
    stVE <- iSD %*% VE
    stpar <- c(as.vector(stVA[lower.tri(stVA, diag = TRUE)]), as.vector(stVE[lower.tri(stVE, diag = TRUE)]))
    return(stpar)
  }
  grmsem.chol.standcf <- function(estimates) {
    #A part
    A <- diag(k)
    A[lower.tri(A, diag = TRUE)] <- estimates[a.sel]
    VA <- A %*% t(A)
    #E part
    E <- diag(k)
    E[lower.tri(E, diag = TRUE)] <- estimates[e.sel]
    VE <- E %*% t(E)
    VP <- VA + VE
    iSD <- try(solve(sqrt(diag(k) * VP)), silent = TRUE)
    #Standardised coefficients
    stA <- iSD %*% A
    stE <- iSD %*% E
    stpar <- c(as.vector(stA[lower.tri(stA, diag = TRUE)]), as.vector(stE[lower.tri(stE, diag = TRUE)]))
    return(stpar)
  }
  #Generate standardised coefficients
  grmsem.ind.standcf <- function(estimates) {
    #A part
    A.c <- length(k)
    A.c <- estimates[a.sel][1:k]
    A.s <- diag(k)
    diag(A.s) <- estimates[a.sel][(k + 1):(2 * k)]
    VA <- A.c %*% t(A.c) + A.s %*% t(A.s)
    #E part
    E.c <- length(k)
    E.c <- estimates[e.sel][1:k]
    E.s <- diag(k)
    diag(E.s) <- estimates[e.sel][(k + 1):(2 * k)]
    VE <- E.c %*% t(E.c) + E.s %*% t(E.s)
    VP <- VA + VE
    iSD <- try(solve(sqrt(diag(k) * VP)), silent = TRUE)
    #Standardised coefficients
    stA.c <- iSD %*% A.c
    stA.s <- iSD %*% A.s
    stE.c <- iSD %*% E.c
    stE.s <- iSD %*% E.s
    stpar <- c(stA.c, diag(stA.s), stE.c, diag(stE.s))
    return(stpar)
  }
  grmsem.indchol.standcf <- function(estimates) {
    #A part
    A.c <- length(k)
    A.c <- estimates[a.sel][1:k]
    A.s <- diag(k)
    diag(A.s) <- estimates[a.sel][(k + 1):(2 * k)]
    VA <- A.c %*% t(A.c) + A.s %*% t(A.s)
    #E part
    E <- diag(k)
    E[lower.tri(E, diag = TRUE)] <- estimates[e.sel]
    VE <- E %*% t(E)
    VP <- VA + VE
    iSD <- try(solve(sqrt(diag(k) * VP)), silent = TRUE)
    #Standardised coefficients
    stA.c <- iSD %*% A.c
    stA.s <- iSD %*% A.s
    stE <- iSD %*% E
    stpar <- c(stA.c, diag(stA.s), as.vector(stE[lower.tri(stE, diag = TRUE)]))
    return(stpar)
  }
  ################################################################################
  #Main body
  ################################################################################
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Full length
  if (class(grmsem.out) != "grmsem.fit") {
    stop('The function requires a grmsem.fit object.')
  }
  if (is.null(grmsem.out$model.out)) {
    stop('There is no model.out information containing parameter estimates.')
  }
  notfree <- which(grmsem.out$model.in$freepar == 0)
  free <- which(grmsem.out$model.in$freepar == 1)
  label <- grmsem.out$model.in$label
  estimates <- grmsem.out$model.out$estimates
  #Adjust estimates (full vector)
  estimates[notfree] <- 0
  if (any(is.na(estimates))) {
    stop('Some parameter estimates in model.out contain missing values.')
  }
  model.in <- grmsem.out$model.in
  k <- grmsem.out$k
  n <- grmsem.out$n #total number of observations
  model <- grmsem.out$model
  ph.nms <- grmsem.out$ph.nms
  a.sel <- which(grmsem.out$model.in$part == "a")
  e.sel <- which(grmsem.out$model.in$part == "e")
  standcf <- NULL
  jac.standcf <- NULL
  if (model == "DS") {
    standcf <- grmsem.ds.standcf(estimates)
    jac.standcf <- numDeriv::jacobian(grmsem.ds.standcf, x = standcf)
  } else if (model == "Cholesky") {
    standcf <- grmsem.chol.standcf(estimates)
    jac.standcf <- numDeriv::jacobian(grmsem.chol.standcf, x = standcf)
  } else if (model == "IP") {
    standcf <- grmsem.ind.standcf(estimates)
    jac.standcf <- numDeriv::jacobian(grmsem.ind.standcf, x = standcf)
  } else if (model == "IPC") {
    standcf <- grmsem.indchol.standcf(estimates)
    jac.standcf <- numDeriv::jacobian(grmsem.indchol.standcf, x = standcf)
  } 
  #Free parameters only
  stand.model.out <- NULL
  stVCOV <- NULL
  stVCOV <- try((jac.standcf)[free, free] %*% grmsem.out$VCOV %*% t(jac.standcf)[free, free], silent = TRUE)
  standse <- rep(NA, length(label))
  model.standse <- sqrt(diag(stVCOV))
  standse[free] <- model.standse
  Z <- rep(NA, length(label))
  pZ <- rep(NA, length(label))
  Z[free] <- standcf[free]/standse[free]
  pZ[free] <- 2 * stats::pnorm(-abs(Z[free]))
  standcf[notfree] <- NA
  colnames(stVCOV) <- label[free]
  rownames(stVCOV) <- label[free]
  stand.model.out <- data.frame(label = label, stand.estimates = standcf, stand.se = standse, Z = Z, p = pZ)
  stout <- list(model.in = model.in, stand.model.out = stand.model.out, stVCOV = stVCOV, k = k, n = n, model = model, ph.nms = ph.nms) #EDIT
  class(stout) <- ("grmsem.stpar")
  return(stout)
}
