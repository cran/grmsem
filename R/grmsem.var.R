#' grmsem variance estimation function
#'
#' This function estimates genetic and residual variances, and genetic correlations.
#' @details The \code{grmsem.var} function can be used to estimate genetic and residual covariance and correlations for DS, Cholesky, IP and IPC models, based on 
#' \code{grmsem.fit} or \code{grmsem.stpar} objects. For the latter, the diagonal elements of the VA output matrix 
#' detail the heritabilities. Except for directly estimated variance components using the DS model, all standard errors 
#' are derived with the Delta method. 
#' @param grmsem.out A grmsem.fit or grmsem.stpar object. Default NULL.
#' @keywords grmsem
#' @export 
#' @return \code{grmsem.var} returns a list object consisting of the following matrices:
#' \item{VA}{estimated genetic variance}
#' \item{VA.se}{standard error of estimated genetic variance}
#' \item{VE}{estimated residual variance}
#' \item{VE.se}{standard error of estimated residual variance}
#' \item{VP}{estimated total phenotypic variance}
#' \item{RG}{genetic correlation}
#' \item{RG.se}{standard error genetic correlation}
#' \item{RE}{residual correlation}
#' \item{RG.se}{standard error residual correlation}
#' @examples
#' #(runtime should be less than one minute)
#' \donttest{
#' out <- grmsem.fit(ph.small, G.small, LogL = TRUE, estSE = TRUE)
#' var.out <- grmsem.var(out)
#' print(var.out)}
grmsem.var <- function(grmsem.out = NULL){
  local
  ################################################################################
  #Main body
  ################################################################################
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Full length
  if (!((class(grmsem.out) %in% c("grmsem.fit","grmsem.stpar")))) {
    stop('The function requires a grmsem.fit or grmsem.stpar object.')
  }
  if (class(grmsem.out) == "grmsem.fit" & is.null(grmsem.out$model.out)) {
    stop('There is no model.out information available containing parameter estimates.')
  }
  if (class(grmsem.out) == "grmsem.stpar" & is.null(grmsem.out$stand.model.out)) {
    stop('There is no stand.model.out information available containing parameter estimates.')
  }
  notfree <- which(grmsem.out$model.in$freepar == 0)
  free <- which(grmsem.out$model.in$freepar == 1)
  label <- grmsem.out$model.in$label
  k <- grmsem.out$k
  ph.nms <- grmsem.out$ph.nms
  model <- grmsem.out$model
  model.vector <- c("Cholesky", "IP","IPC", "DS")
  if (!is.element(model, model.vector)) {
    m.msg <- paste("Variances are only reported for Cholesky, IP, IPC and DS models.", sep = "")
    stop(m.msg)
  }
  a.sel <- which(grmsem.out$model.in$part == "a")
  e.sel <- which(grmsem.out$model.in$part == "e")
  if(class(grmsem.out) == "grmsem.fit") {
    estimates <- grmsem.out$model.out$estimates
    VCOV <- grmsem.out$VCOV
    grmsem.out$model.out<-grmsem.out$model.out
    if(model == "DS") {
      estimates.se <- grmsem.out$model.out$se
    }
  }else{
    estimates <- grmsem.out$stand.model.out$stand.estimates
    VCOV <- grmsem.out$stVCOV
    grmsem.out$model.out<-grmsem.out$stand.model.out
    if(model == "DS") {
      estimates.se <- grmsem.out$stand.model.out$stand.se
    }
  }
  #Adjust estimates (wrt full vector)
  estimates[notfree] <- 0
  if (any(is.na(estimates))) {
    stop('Parameter estimates contain missing values.')
  }
  if (any(is.na(VCOV))) {
    stop('VCOV contains missing values.')
  }
  #Estimated variance components
  VA <- NULL # Genetic var/covar
  VE <- NULL # Residual var/covar
  VP <- NULL # Derived phenotypic var/covar
  VA.se <- NULL # Genetic var/covar se
  VE.se <- NULL # Residual var/covar se
  RG <- NULL # Genetic correlation
  RG.se <- NULL # Genetic correlation se
  RE <- NULL # Residual correlation
  RE.se <- NULL # Residual correlation se
  if (model == "DS") {
    #A part
    VA <- diag(k)
    VA[lower.tri(VA, diag = TRUE)] <- estimates[a.sel]
    VA[upper.tri(VA)] <- t(VA[lower.tri(VA)])
    VA.se <- diag(k)
    VA.se[lower.tri(VA.se, diag = TRUE)] <- estimates.se[a.sel]
    VA.se[upper.tri(VA.se)] <- t(VA.se[lower.tri(VA.se)])
    #E part
    VE <- diag(k)
    VE[lower.tri(VE, diag = TRUE)] <- estimates[e.sel]
    VE[upper.tri(VE)] <- t(VE[lower.tri(VE)])
    VE.se <- diag(k)
    VE.se[lower.tri(VE.se, diag = TRUE)] <- estimates.se[e.sel]
    VE.se[upper.tri(VE.se)] <- t(VE.se[lower.tri(VE.se)])
    VP <- VA + VE
    #Calculate standard errors
    #Free parameters
    iSDVA <- try(solve(sqrt(diag(k) * VA)), silent = TRUE)
    iSDVE <- try(solve(sqrt(diag(k) * VE)), silent = TRUE)
    #Genetic correlation: RG
    RG <- iSDVA %*% VA %*% iSDVA
    #labels
    xlabel <- paste(rep("x", length(label[free])), seq_along(label[free]), sep = "")
    #Full label
    full.xlabel <- rep(0, length(label))
    full.xlabel[free] <- xlabel
    XA <- diag(k)
    diag(XA) <- 0
    XA[lower.tri(XA, diag = TRUE)] <- full.xlabel[a.sel]
    XA[upper.tri(XA)] <- t(XA[lower.tri(XA)])
    RG.se.formula <- matrix(nrow = k, ncol = k)
    #RG.se.formula <- as.data.frame(RG.se.formula)
    RG.se <- RG.se.formula
    #RG se formula
    RG.tmp.formula <- rep(NA, length(k))
    for (i in 1:k) {
      for (j in 1:k) {
        if (i == j) {
          RG.tmp.formula[i] <- paste("(", XA[i,j], ")", sep = "")
        }
      }
    }
    for (i in 1:k) {
      for (j in 1:k) {
        RG.se.formula[i, j] <- paste("msm::deltamethod(~(", XA[i,j], ")/sqrt(", paste(RG.tmp.formula[i], RG.tmp.formula[j], sep = " * "), "), estimates[free], VCOV)", sep = "")
        RG.se[i, j] <- eval(parse(text = RG.se.formula[i, j]))
      }
    }
    #print(RG.se.formula)
    #Residual correlation: RE
    RE <- iSDVE %*% VE %*% iSDVE
    #Res correlation se: RG.se
    XE <- diag(k)
    diag(XE) <- 0
    XE[lower.tri(XE, diag = TRUE)] <- full.xlabel[e.sel]
    XE[upper.tri(XE)] <- t(XE[lower.tri(XE)])
    RE.se.formula <- matrix(nrow = k, ncol = k)
    #RE.se.formula <- as.data.frame(RE.se.formula)
    RE.se <- RE.se.formula
    RE.tmp.formula <- rep(NA, length(k))
    for (i in 1:k) {
      for (j in 1:k) {
        if (i == j) {
          RE.tmp.formula[i] <- paste("(", XE[i,j], ")", sep = "")
        }
      }
    }
    for (i in 1:k) {
      for (j in 1:k) {
        RE.se.formula[i, j] <- paste("msm::deltamethod(~(", XE[i,j], ")/sqrt(", paste(RE.tmp.formula[i], RE.tmp.formula[j], sep = " * "), "), estimates[free], VCOV)", sep = "")
        RE.se[i, j] <- eval(parse(text = RE.se.formula[i, j]))
      }
    }
  #print(RE.se.formula)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  } else if (model == "Cholesky") {
    #A part
    A <- diag(k)
    A[lower.tri(A, diag = TRUE)] <- estimates[a.sel]
    VA <- A %*% t(A)
    #E part
    E <- diag(k)
    E[lower.tri(E, diag = TRUE)] <- estimates[e.sel]
    VE <- E %*% t(E)
    VP <- VA + VE
    #Calculate standard errors
    #Free parameters
    xlabel <- paste(rep("x", length(label[free])), seq_along(label[free]), sep = "")
    #Full label
    full.xlabel <- rep(0, length(label))
    full.xlabel[free] <- xlabel
    XA <- diag(k)
    diag(XA) <- 0
    XAT <- XA
    XA[lower.tri(XA, diag = TRUE)] <- full.xlabel[a.sel]
    XAT <- t(XA)
    VA.se.formula <- matrix(nrow = k, ncol = k)
    #VA.se.formula <- as.data.frame(VA.se.formula)
    VA.se <- VA.se.formula
    for (i in 1:k) {
      for (j in 1:k) {
        VA.se.formula[i, j] <- paste("msm::deltamethod(~", paste(paste(XA[i, ], XAT[, j], sep = " * "), collapse = " + "), ", estimates[free], VCOV)", sep = "")
        VA.se[i, j] <- eval(parse(text = VA.se.formula[i, j]))
      }
    }
    XE <- diag(k)
    diag(XE) <- 0
    XET <- XE
    XE[lower.tri(XE, diag = TRUE)] <- full.xlabel[e.sel]
    XET <- t(XE)
    VE.se.formula <- matrix(nrow = k, ncol = k)
    #VE.se.formula <- as.data.frame(VE.se.formula)
    VE.se <- VE.se.formula
    for (i in 1:k) {
      for (j in 1:k) {
        VE.se.formula[i, j] <- paste("msm::deltamethod(~", paste(paste(XE[i, ], XET[, j], sep = " * "), collapse = " + "), ", estimates[free], VCOV)", sep = "")
        VE.se[i, j] <- eval(parse(text = VE.se.formula[i, j]))
      }
    }
    iSDVA <- try(solve(sqrt(diag(k) * VA)), silent = TRUE)
    iSDVE <- try(solve(sqrt(diag(k) * VE)), silent = TRUE)
    #Genetic correlation: RG
    RG <- iSDVA %*% VA %*% iSDVA
    RG.se.formula <- matrix(nrow = k, ncol = k)
    #RG.se.formula <- as.data.frame(RG.se.formula)
    RG.se <- RG.se.formula
    RG.tmp.formula <- rep(NA, length(k))
    for (i in 1:k) {
      for (j in 1:k) {
        if (i == j) {
          RG.tmp.formula[i] <- paste("(", paste(paste(XA[i, ], XAT[, j], sep = " * "), collapse = " + "), ")", sep = "")
        }
      }
    }
    for (i in 1:k) {
      for (j in 1:k) {
        RG.se.formula[i, j] <- paste("msm::deltamethod(~(", paste(paste(XA[i, ], XAT[, j], sep = " * "), collapse = " + "), ")/sqrt(", paste(RG.tmp.formula[i], RG.tmp.formula[j], sep = " * "), "), estimates[free], VCOV)", sep = "")
        RG.se[i, j] <- eval(parse(text = RG.se.formula[i, j]))
      }
    }
    #Residual correlation: RE
    RE <- iSDVE %*% VE %*% iSDVE
    RE.se.formula <- matrix(nrow = k, ncol = k)
    #RE.se.formula <- as.data.frame(RE.se.formula)
    RE.se <- RE.se.formula
    RE.tmp.formula <- rep(NA, length(k))
    for (i in 1:k) {
      for (j in 1:k) {
        if (i == j) {
          RE.tmp.formula[i] <- paste("(", paste(paste(XE[i, ], XET[, j], sep = " * "), collapse = " + "), ")", sep = "")
        }
      }
    }
    for (i in 1:k) {
      for (j in 1:k) {
        RE.se.formula[i, j] <- paste("msm::deltamethod(~(", paste(paste(XE[i, ], XET[, j], sep = " * "), collapse = " + "), ")/sqrt(", paste(RE.tmp.formula[i], RE.tmp.formula[j], sep = " * "), "), estimates[free], VCOV)", sep = "")
        RE.se[i, j] <- eval(parse(text = RE.se.formula[i, j]))
      }
    }
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  } else if (model == "IP") {
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
    #Calculate standard errors
    #Free parameters
    xlabel <- paste(rep("x", length(label[free])), seq_along(label[free]), sep = "")
    #Full label
    full.xlabel <- rep(0, length(label))
    full.xlabel[free] <- xlabel
    XA.c <- length(k)
    XA.c <- full.xlabel[a.sel][1:k]
    XA.s <- diag(k)
    diag(XA.s) <- full.xlabel[a.sel][(k + 1):(2 * k)]
    VA.se.formula <- matrix(nrow = k, ncol = k)
    #VA.se.formula <- as.data.frame(VA.se.formula)
    VA.se.formula.c <- VA.se.formula
    VA.se.formula.s <- VA.se.formula
    VA.se <- VA.se.formula
    for (i in 1:k) {
      for (j in 1:k) {
        VA.se.formula.c[i, j] <- paste(XA.c[i], t(XA.c[j]), sep = " * ")
        VA.se.formula.s[i, j] <- 0
        if (i == j) {
          VA.se.formula.s[i, j] <- paste(XA.s[i, j], XA.s[i, j], sep = " * ")
        }
        VA.se.formula[i, j] <- paste("msm::deltamethod(~", paste(VA.se.formula.c[i, j], VA.se.formula.s[i, j], sep = " + "), ", estimates[free], VCOV)", sep = "")
        VA.se[i, j] <- eval(parse(text = VA.se.formula[i, j]))
      }
    }
    #print(VA.se.formula)
    XE.c <- length(k)
    XE.c <- full.xlabel[e.sel][1:k]
    XE.s <- diag(k)
    diag(XE.s) <- full.xlabel[e.sel][(k + 1):(2 * k)]
    VE.se.formula <- matrix(nrow = k, ncol = k)
    #VE.se.formula <- as.data.frame(VE.se.formula)
    VE.se.formula.c <- VE.se.formula
    VE.se.formula.s <- VE.se.formula
    VE.se <- VE.se.formula
    for (i in 1:k) {
      for (j in 1:k) {
        VE.se.formula.c[i, j] <- paste(XE.c[i], t(XE.c[j]), sep = " * ")
        VE.se.formula.s[i, j] <- 0
        if (i == j) {
          VE.se.formula.s[i, j] <- paste(XE.s[i, j], XE.s[i, j], sep = " * ")
        }
        VE.se.formula[i, j] <- paste("msm::deltamethod(~", paste(VE.se.formula.c[i, j], VE.se.formula.s[i, j], sep = " + "), ", estimates[free], VCOV)", sep = "")
        VE.se[i, j] <- eval(parse(text = VE.se.formula[i, j]))
      }
    }
    #print(VE.se.formula)
    iSDVA <- try(solve(sqrt(diag(k) * VA)), silent = TRUE)
    iSDVE <- try(solve(sqrt(diag(k) * VE)), silent = TRUE)
    #Genetic correlation: RG
    RG <- iSDVA %*% VA %*% iSDVA
    RG.se.formula <- matrix(nrow = k, ncol = k)
    #RG.se.formula <- as.data.frame(RG.se.formula)
    RG.se <- RG.se.formula
    RG.tmp.formula.num <- RG.se.formula
    RG.tmp.formula.denom <- rep(NA, length(k))
    for (i in 1:k) {
      for (j in 1:k) {
        RG.tmp.formula.num[i, j] <- paste("(", paste(VA.se.formula.c[i, j], VA.se.formula.s[i, j], sep = " + "), ")", sep = "")
        if (i == j) {
          RG.tmp.formula.denom[i] <- paste("(", paste(VA.se.formula.c[i, j], VA.se.formula.s[i, j], sep = " + "), ")", sep = "")
        }
      }
    }
    for (i in 1:k) {
      for (j in 1:k) {
        RG.se.formula[i, j] <- paste("msm::deltamethod(~(", RG.tmp.formula.num[i, j], ")/sqrt(", paste(RG.tmp.formula.denom[i], RG.tmp.formula.denom[j], sep = " * "), "), estimates[free], VCOV)", sep = "")
        RG.se[i, j] <- eval(parse(text = RG.se.formula[i, j]))
      }
    }
    #Residual correlation: RE
    RE <- iSDVE %*% VE %*% iSDVE
    RE.se.formula <- matrix(nrow = k, ncol = k)
    #RE.se.formula <- as.data.frame(RE.se.formula)
    RE.se <- RE.se.formula
    RE.tmp.formula.num <- RE.se.formula
    RE.tmp.formula.denom <- rep(NA, length(k))
    for (i in 1:k) {
      for (j in 1:k) {
        RE.tmp.formula.num[i, j] <- paste("(", paste(VE.se.formula.c[i, j], VE.se.formula.s[i, j], sep = " + "), ")", sep = "")
        if (i == j) {
          RE.tmp.formula.denom[i] <- paste("(", paste(VE.se.formula.c[i, j], VE.se.formula.s[i, j], sep = " + "), ")", sep = "")
        }
      }
    }
    for (i in 1:k) {
      for (j in 1:k) {
        RE.se.formula[i, j] <- paste("msm::deltamethod(~(", RE.tmp.formula.num[i, j], ")/sqrt(", paste(RE.tmp.formula.denom[i], RE.tmp.formula.denom[j], sep = " * "), "), estimates[free], VCOV)", sep = "")
        RE.se[i, j] <- eval(parse(text = RE.se.formula[i, j]))
      }
    }
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  } else if (model == "IPC") { # "IPC")
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
    #Free parameters
    xlabel <- paste(rep("x", length(label[free])), seq_along(label[free]), sep = "")
    #Full label
    full.xlabel <- rep(0, length(label))
    full.xlabel[free] <- xlabel
    XA.c <- length(k)
    XA.c <- full.xlabel[a.sel][1:k]
    XA.s <- diag(k)
    diag(XA.s) <- full.xlabel[a.sel][(k + 1):(2 * k)]
    VA.se.formula <- matrix(nrow = k, ncol = k)
    #VA.se.formula <- as.data.frame(VA.se.formula)
    VA.se.formula.c <- VA.se.formula
    VA.se.formula.s <- VA.se.formula
    VA.se <- VA.se.formula
    for (i in 1:k) {
      for (j in 1:k) {
        VA.se.formula.c[i, j] <- paste(XA.c[i], t(XA.c[j]), sep = " * ")
        VA.se.formula.s[i, j] <- 0
        if (i == j) {
          VA.se.formula.s[i, j] <- paste(XA.s[i, j], XA.s[i, j], sep = " * ")
        }
        VA.se.formula[i, j] <- paste("msm::deltamethod(~", paste(VA.se.formula.c[i, j], VA.se.formula.s[i, j], sep = " + "), ", estimates[free], VCOV)", sep = "")
        VA.se[i, j] <- eval(parse(text = VA.se.formula[i, j]))
      }
    }
    #print(VA.se.formula)
    XE <- diag(k)
    diag(XE) <- 0
    XET <- XE
    XE[lower.tri(XE, diag = TRUE)] <- full.xlabel[e.sel]
    XET <- t(XE)
    VE.se.formula <- matrix(nrow = k, ncol = k)
    #VE.se.formula <- as.data.frame(VE.se.formula)
    VE.se <- VE.se.formula
    for (i in 1:k) {
      for (j in 1:k) {
        VE.se.formula[i, j] <- paste("msm::deltamethod(~", paste(paste(XE[i, ], XET[, j], sep = " * "), collapse = " + "), ", estimates[free], VCOV)", sep = "")
        VE.se[i, j] <- eval(parse(text = VE.se.formula[i, j]))
      }
    }
    #print(VE.se.formula)
    iSDVA <- try(solve(sqrt(diag(k) * VA)), silent = TRUE)
    iSDVE <- try(solve(sqrt(diag(k) * VE)), silent = TRUE)
    #Genetic correlation: RG
    RG <- iSDVA %*% VA %*% iSDVA
    RG.se.formula <- matrix(nrow = k, ncol = k)
    #RG.se.formula <- as.data.frame(RG.se.formula)
    RG.se <- RG.se.formula
    RG.tmp.formula.num <- RG.se.formula
    RG.tmp.formula.denom <- rep(NA, length(k))
    for (i in 1:k) {
      for (j in 1:k) {
        RG.tmp.formula.num[i, j] <- paste("(", paste(VA.se.formula.c[i, j], VA.se.formula.s[i, j], sep = " + "), ")", sep = "")
        if (i == j) {
          RG.tmp.formula.denom[i] <- paste("(", paste(VA.se.formula.c[i, j], VA.se.formula.s[i, j], sep = " + "), ")", sep = "")
        }
      }
    }
    for (i in 1:k) {
      for (j in 1:k) {
        RG.se.formula[i, j] <- paste("msm::deltamethod(~(", RG.tmp.formula.num[i, j], ")/sqrt(", paste(RG.tmp.formula.denom[i], RG.tmp.formula.denom[j], sep = " * "), "), estimates[free], VCOV)", sep = "")
        RG.se[i, j] <- eval(parse(text = RG.se.formula[i, j]))
      }
    }
    #Residual correlation: RE
    RE <- iSDVE %*% VE %*% iSDVE
    RE.se.formula <- matrix(nrow = k, ncol = k)
    #RE.se.formula <- as.data.frame(RE.se.formula)
    RE.se <- RE.se.formula
    RE.tmp.formula <- rep(NA, length(k))
    for (i in 1:k) {
      for (j in 1:k) {
        if (i == j) {
          RE.tmp.formula[i] <- paste("(", paste(paste(XE[i, ], XET[, j], sep = " * "), collapse = " + "), ")", sep = "")
        }
      }
    }
    for (i in 1:k) {
      for (j in 1:k) {
        RE.se.formula[i, j] <- paste("msm::deltamethod(~(", paste(paste(XE[i, ], XET[, j], sep = " * "), collapse = " + "), ")/sqrt(", paste(RE.tmp.formula[i], RE.tmp.formula[j], sep = " * "), "), estimates[free], VCOV)", sep = "")
        RE.se[i, j] <- eval(parse(text = RE.se.formula[i, j]))
      }
    }	
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  }
  
  VA.se <- as.matrix(VA.se)
  VE.se <- as.matrix(VE.se)
  colnames(VA) <- seq(1:k)
  colnames(VA.se) <- seq(1:k)
  colnames(VE) <- seq(1:k)
  colnames(VE.se) <- seq(1:k)
  colnames(VP) <- seq(1:k)
  colnames(RG) <- seq(1:k)
  colnames(RG.se) <- seq(1:k)
  colnames(RE) <- seq(1:k)
  colnames(RE.se) <- seq(1:k)
  rownames(VA) <- ph.nms
  rownames(VA.se) <- ph.nms
  rownames(VE) <- ph.nms
  rownames(VE.se) <- ph.nms
  rownames(VP) <- ph.nms
  rownames(RG) <- ph.nms
  rownames(RG.se) <- ph.nms
  rownames(RE) <- ph.nms
  rownames(RE.se) <- ph.nms
  out <- (list(VA = VA, VA.se = VA.se, VE = VE, VE.se = VE.se, VP = VP, RG = RG, RG.se = RG.se, RE = RE, RE.se = RE.se))
  class(out) <- ("grmsem.var")
  return(out)
}
