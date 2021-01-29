#' grmsem factorial co-heritability and co-environmentality estimation function.
#'
#' This function estimates factorial co-heritabilities and factorial co-environmentalities.
#' @details The \code{grmsem.fcoher} function can be used to estimate factorial co-heritabilities and factorial co-environmentalities 
#' for models estimating latent variables (Cholesky, IP or IPC models), based on \code{grmsem.fit} or \code{grmsem.stpar} objects.
#' The factorial co-heritability of a genetic factor m for trait t is the ratio of the genetic variance explained by factor m (A_mt) to the total genetic variance (A_t): A_mt / A_t. 
#' The factorial co-environmentality of a residual factor n for trait t is the ratio of the residual variance explained by factor n (E_nt) to the total residual variance (E_t):  E_nt / E_t.
#' All standard errors are derived with the Delta method.
#' @param grmsem.out grmsem.fit or grmsem.stpar object. Default NULL.
#' @keywords grmsem
#' @export
#' @return \code{grmsem.fcoher} returns an extended model.out dataframe, fcoher.model.out, with the following columns:
#' \item{label}{parameter label}
#' \item{estimates}{estimated parameters}
#' \item{gradient}{gradient}
#' \item{se}{SE}
#' \item{Z}{Z (Wald) of factor loading}
#' \item{p}{p (Wald) of factor loading}
#' \item{Vi}{squared factor loading, explained phenotypic variation by the factor}
#' \item{Vi.se}{SE of squared factor loading}
#' \item{FCOHER}{factorial co-heritability}
#' \item{FCOHER.se}{SE of factorial co-heritability}
#' \item{FCOHER.Z}{Z (wald) of factorial co-heritability}
#' \item{FCOHER.se}{p (Wald) of factorial co-heritability}
#' \item{FCOENV}{factorial co-environmentality}
#' \item{FCOENV.se}{SE of factorial co-environmentality}
#' \item{FCOENV.Z}{Z (wald) of factorial co-environmentality}
#' \item{FCOENV.se}{p (Wald) of factorial co-environmentality}
#' @examples
#' #(runtime should be less than one minute)
#' \donttest{
#' out <- grmsem.fit(ph.small, G.small, LogL = TRUE, estSE = TRUE)
#' grmsem.fcoher(out)}
grmsem.fcoher <- function(grmsem.out = NULL){
  local
  ################################################################################
  #Main body
  ################################################################################
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
  if(class(grmsem.out) == "grmsem.fit") {
    estimates <- grmsem.out$model.out$estimates
    VCOV <- grmsem.out$VCOV
    grmsem.out$model.out<-grmsem.out$model.out
  }else{
    estimates <- grmsem.out$stand.model.out$stand.estimates
    VCOV <- grmsem.out$stVCOV
    grmsem.out$model.out<-grmsem.out$stand.model.out
  }
  #Adjust estimates (wrt full vector)
  estimates[notfree] <- 0
  if (any(is.na(estimates))) {
    stop('Parameter estimates contain missing values.')
  }
  k <- grmsem.out$k
  if (k<2) {
    stop('The estimation of factorial co-heritabilities requires at least 2 phenotypes.')
  }
  if (any(is.na(VCOV))) {
    stop('VCOV contains missing values.')
  }
  model <- grmsem.out$model
  model.vector <- c("Cholesky", "IP", "IPC")
  if (!is.element(model, model.vector)) {
    m.msg <- paste("Please provide a grmsem.fit object for a Cholesky, Independent Pathway or IPC model (factor loadings).", sep = "")
    stop(m.msg)
  }
  a.sel <- which(grmsem.out$model.in$part == "a")
  e.sel <- which(grmsem.out$model.in$part == "e")
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (model == "Cholesky") {
    #Derived parameters
    FCOHER.se.formula <- matrix(nrow = k, ncol = k)
    FCOHER.se.formula <- as.data.frame(FCOHER.se.formula)
    AA.se <- AA.se.formula <- AA.se.tmp.formula <- FCOHER <- FCOHER.se <- FCOHER.se.tmp.formula <- FCOHER.se.formula
    FCOENV <- FCOHER
    # Factorial coheritability
    A <- diag(k)
    A[lower.tri(A, diag = TRUE)] <- estimates[a.sel]
    VA <- A %*% t(A)
    AA <- A*A
    FCOHER <- AA/diag(VA)
    #E part  - genetic path coefficients
    E <- diag(k)
    E[lower.tri(E, diag = TRUE)] <- estimates[e.sel]
    VE <- E %*% t(E)
    EE <- E*E
    FCOENV <- EE/diag(VE)
    #Calculate standard errors
    #Free parameters
    xlabel <- paste(rep("x", length(label[free])), seq_along(label[free]), sep = "")
    #Full label
    full.xlabel <- rep(0, length(label))
    full.xlabel[free] <- xlabel
    XA <- diag(k)
    diag(XA) <- 0
    XA[lower.tri(XA, diag = TRUE)] <- full.xlabel[a.sel]
    for (i in 1:k) {
      for (j in 1:k) {
        AA.se.tmp.formula[i, j] <- (paste(XA[i,j ], XA[i, j], sep = " * "))
      }
    }
    for (i in 1:k) {
      for (j in 1:k) {
        FCOHER.se.tmp.formula[i, j] <- (paste(AA.se.tmp.formula[i, ], collapse = " + "))
      }
    }
    for (i in 1:k) {
      for (j in 1:k) {
        AA.se.formula[i, j] <- paste("msm::deltamethod(~", AA.se.tmp.formula[i, j], ", estimates[free], VCOV)", sep = "")
        AA.se[i, j] <- eval(parse(text = AA.se.formula[i, j]))
        FCOHER.se.formula[i, j] <- paste("msm::deltamethod(~", paste("(", AA.se.tmp.formula[i, j], ")/(", FCOHER.se.tmp.formula[i, j], ")", sep=""), ", estimates[free], VCOV)", sep = "")
        FCOHER.se[i, j] <- eval(parse(text = FCOHER.se.formula[i, j]))
      }
    }
    # Factorial coenvironmentality
    XE <- diag(k)
    diag(XE) <- 0
    XE[lower.tri(XE, diag = TRUE)] <- full.xlabel[e.sel]
    FCOENV.se.formula <- matrix(nrow = k, ncol = k)
    FCOENV.se.formula <- as.data.frame(FCOENV.se.formula)
    EE.se <- EE.se.formula <- EE.se.tmp.formula <- FCOENV.se <- FCOENV.se.tmp.formula  <- FCOENV.se.formula
    for (i in 1:k) {
      for (j in 1:k) {
        EE.se.tmp.formula[i, j] <- (paste(XE[i,j ], XE[i, j], sep = " * "))
      }
    }
    for (i in 1:k) {
      for (j in 1:k) {
        FCOENV.se.tmp.formula[i, j] <- (paste(EE.se.tmp.formula[i, ], collapse = " + "))
      }
    }
    for (i in 1:k) {
      for (j in 1:k) {
        EE.se.formula[i, j] <- paste("msm::deltamethod(~", EE.se.tmp.formula[i, j], ", estimates[free], VCOV)", sep = "")
        EE.se[i, j] <- eval(parse(text = EE.se.formula[i, j]))
        FCOENV.se.formula[i, j] <- paste("msm::deltamethod(~", paste("(", EE.se.tmp.formula[i, j], ")/(", FCOENV.se.tmp.formula[i, j], ")", sep=""), ", estimates[free], VCOV)", sep = "")
        FCOENV.se[i, j] <- eval(parse(text = FCOENV.se.formula[i, j]))
      }
    }
    #Updated model.out dataframe accounting for missing parameters
    grmsem.out$model.out$Vi <- NA
    grmsem.out$model.out$Vi[a.sel]<-AA[lower.tri(AA, diag = TRUE)]
    grmsem.out$model.out$Vi[e.sel]<-EE[lower.tri(EE, diag = TRUE)]
    grmsem.out$model.out$Vi.se <- NA
    grmsem.out$model.out$Vi.se[a.sel]<-AA.se[lower.tri(AA.se, diag = TRUE)]
    grmsem.out$model.out$Vi.se[e.sel]<-EE.se[lower.tri(EE.se, diag = TRUE)]
    grmsem.out$model.out$FCOHER <- NA
    grmsem.out$model.out$FCOHER.se <- NA
    grmsem.out$model.out$FCOHER[a.sel] <- FCOHER[lower.tri(FCOHER, diag = TRUE)]
    grmsem.out$model.out$FCOHER.se[a.sel]<-FCOHER.se[lower.tri(FCOHER.se, diag = TRUE)]
    grmsem.out$model.out$FCOHER[notfree]<-NA
    grmsem.out$model.out$FCOHER.se[notfree]<-NA
    grmsem.out$model.out$FCOHER.Z<-grmsem.out$model.out$FCOHER/grmsem.out$model.out$FCOHER.se
    grmsem.out$model.out$FCOHER.p<-2*stats::pnorm(-abs(grmsem.out$model.out$FCOHER.Z))
    grmsem.out$model.out$FCOENV <- NA
    grmsem.out$model.out$FCOENV.se <- NA
    grmsem.out$model.out$FCOENV[e.sel] <- FCOENV[lower.tri(FCOENV, diag = TRUE)]
    grmsem.out$model.out$FCOENV.se[e.sel]<-FCOENV.se[lower.tri(FCOENV.se, diag = TRUE)]
    grmsem.out$model.out$FCOENV[notfree]<-NA
    grmsem.out$model.out$FCOENV.se[notfree]<-NA
    grmsem.out$model.out$FCOENV.Z<-grmsem.out$model.out$FCOENV/grmsem.out$model.out$FCOENV.se
    grmsem.out$model.out$FCOENV.p<-2*stats::pnorm(-abs(grmsem.out$model.out$FCOENV.Z))
    grmsem.out$model.out<-as.data.frame(grmsem.out$model.out)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  }else if (model == "IP") {
    #Derived parameters - vectors
    FCOHER.c <- NULL
    FCOHER.s <- NULL
    FCOHER.c.se <- NULL
    FCOHER.s.se <- NULL
    FCOHER.c.se.formula <- NULL
    FCOHER.s.se.formula <- NULL
    FCOENV.c <- NULL
    FCOENV.s <- NULL
    FCOENV.c.se <- NULL
    FCOENV.s.se <- NULL
    FCOENV.c.se.formula <- NULL
    FCOENV.s.se.formula <- NULL
    AA.c.se.formula <- NULL
    AA.s.se.formula <- NULL
    EE.c.se.formula <- NULL
    EE.s.se.formula <- NULL
    AA.c.se <- NULL
    AA.s.se <- NULL
    EE.c.se <- NULL
    EE.s.se <- NULL
    #Factorial co-heritabilities
    A.c <- length(k)
    A.c <- estimates[a.sel][1:k]
    A.s <- diag(k)
    diag(A.s) <- estimates[a.sel][(k + 1):(2 * k)]
    VA <- A.c %*% t(A.c) + A.s %*% t(A.s)
    AA.c <- A.c*A.c
    AA.s <- A.s*A.s
    FCOHER.c <- AA.c/diag(VA)
    FCOHER.s <- AA.s/diag(VA)
    #E part
    E.c <- length(k)
    E.c <- estimates[e.sel][1:k]
    E.s <- diag(k)
    diag(E.s) <- estimates[e.sel][(k + 1):(2 * k)]
    VE <- E.c %*% t(E.c) + E.s %*% t(E.s)
    EE.c <- E.c*E.c
    EE.s <- E.s*E.s
    FCOENV.c <- EE.c/diag(VE)
    FCOENV.s <- EE.s/diag(VE)
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
    VA.se.formula <- as.data.frame(VA.se.formula)
    VA.se.formula.c <- VA.se.formula.s <- VA.se.formula
    for (i in 1:k) {
      for (j in 1:k) {
        VA.se.formula.c[i, j] <- paste(XA.c[i], t(XA.c[j]), sep = " * ")
        VA.se.formula.s[i, j] <- 0
        if (i == j) {
          VA.se.formula.s[i, j] <- paste(XA.s[i, j], XA.s[i, j], sep = " * ")
        }
        VA.se.formula[i, j] <- paste(VA.se.formula.c[i, j], VA.se.formula.s[i, j], sep = " + ")
      }
    }
    for (i in 1:k) {
      AA.c.se.formula[i] <- paste("msm::deltamethod(~", paste("(", VA.se.formula.c[i,i],")",sep=""), ", estimates[free], VCOV)", sep = "")
      AA.s.se.formula[i] <- paste("msm::deltamethod(~", paste("(", VA.se.formula.s[i,i],")",sep=""), ", estimates[free], VCOV)", sep = "")
      AA.c.se[i] <- eval(parse(text = AA.c.se.formula[i]))
      AA.s.se[i] <- eval(parse(text = AA.s.se.formula[i]))
      FCOHER.c.se.formula[i] <-   paste("msm::deltamethod(~", paste("(", VA.se.formula.c[i,i],")/(", VA.se.formula[i,i], ")",sep=""), ", estimates[free], VCOV)", sep = "")
      FCOHER.c.se[i] <- eval(parse(text = FCOHER.c.se.formula[i]))
      FCOHER.s.se.formula[i] <-   paste("msm::deltamethod(~", paste("(", VA.se.formula.s[i,i],")/(", VA.se.formula[i,i], ")",sep=""), ", estimates[free], VCOV)", sep = "")
      FCOHER.s.se[i] <- eval(parse(text = FCOHER.s.se.formula[i]))
    }
    #Factorial co-environmentalities
    XE.c <- length(k)
    XE.c <- full.xlabel[e.sel][1:k]
    XE.s <- diag(k)
    diag(XE.s) <- full.xlabel[e.sel][(k + 1):(2 * k)]
    VE.se.formula <- matrix(nrow = k, ncol = k)
    VE.se.formula <- as.data.frame(VE.se.formula)
    VE.se.formula.c <- VE.se.formula.s <- VE.se.formula
    for (i in 1:k) {
      for (j in 1:k) {
        VE.se.formula.c[i, j] <- paste(XE.c[i], t(XE.c[j]), sep = " * ")
        VE.se.formula.s[i, j] <- 0
        if (i == j) {
          VE.se.formula.s[i, j] <- paste(XE.s[i, j], XE.s[i, j], sep = " * ")
        }
        VE.se.formula[i, j] <- paste(VE.se.formula.c[i, j], VE.se.formula.s[i, j], sep = " + ")
      }
    }
    for (i in 1:k) {
      EE.c.se.formula[i] <- paste("msm::deltamethod(~", paste("(", VE.se.formula.c[i,i],")",sep=""), ", estimates[free], VCOV)", sep = "")
      EE.s.se.formula[i] <- paste("msm::deltamethod(~", paste("(", VE.se.formula.s[i,i],")",sep=""), ", estimates[free], VCOV)", sep = "")
      EE.c.se[i] <- eval(parse(text = EE.c.se.formula[i]))
      EE.s.se[i] <- eval(parse(text = EE.s.se.formula[i]))
      FCOENV.c.se.formula[i] <-   paste("msm::deltamethod(~", paste("(", VE.se.formula.c[i,i],")/(", VE.se.formula[i,i], ")",sep=""), ", estimates[free], VCOV)", sep = "")
      FCOENV.c.se[i] <- eval(parse(text = FCOENV.c.se.formula[i]))
      FCOENV.s.se.formula[i] <-   paste("msm::deltamethod(~", paste("(", VE.se.formula.s[i,i],")/(", VE.se.formula[i,i], ")",sep=""), ", estimates[free], VCOV)", sep = "")
      FCOENV.s.se[i] <- eval(parse(text = FCOENV.s.se.formula[i]))
    }
    #Updated model.out dataframe accounting for missing parameters
    grmsem.out$model.out$Vi <- NA
    grmsem.out$model.out$Vi[a.sel][1:k]<-AA.c
    grmsem.out$model.out$Vi[a.sel][(k + 1):(2 * k)]<-diag(AA.s)
    grmsem.out$model.out$Vi[e.sel][1:k]<-EE.c
    grmsem.out$model.out$Vi[e.sel][(k + 1):(2 * k)]<-diag(EE.s)
    grmsem.out$model.out$Vi.se <- NA
    grmsem.out$model.out$Vi.se[a.sel][1:k]<-AA.c.se
    grmsem.out$model.out$Vi.se[a.sel][(k + 1):(2 * k)]<-AA.s.se
    grmsem.out$model.out$Vi.se[e.sel][1:k]<-EE.c.se
    grmsem.out$model.out$Vi.se[e.sel][(k + 1):(2 * k)]<-EE.s.se
    grmsem.out$model.out$FCOHER <- NA
    grmsem.out$model.out$FCOHER.se <- NA
    grmsem.out$model.out$FCOHER[a.sel][1:k] <- FCOHER.c
    grmsem.out$model.out$FCOHER[a.sel][(k + 1):(2 * k)] <- diag(FCOHER.s)
    grmsem.out$model.out$FCOHER.se[a.sel][1:k]<-FCOHER.c.se
    grmsem.out$model.out$FCOHER.se[a.sel][(k + 1):(2 * k)]<-FCOHER.s.se
    grmsem.out$model.out$FCOHER[notfree]<-NA
    grmsem.out$model.out$FCOHER.se[notfree]<-NA
    grmsem.out$model.out$FCOHER.Z<-grmsem.out$model.out$FCOHER/grmsem.out$model.out$FCOHER.se
    grmsem.out$model.out$FCOHER.p<-2*stats::pnorm(-abs(grmsem.out$model.out$FCOHER.Z))
    grmsem.out$model.out$FCOENV <- NA
    grmsem.out$model.out$FCOENV.se <- NA
    grmsem.out$model.out$FCOENV[e.sel][1:k] <- FCOENV.c
    grmsem.out$model.out$FCOENV[e.sel][(k + 1):(2 * k)] <- diag(FCOENV.s)
    grmsem.out$model.out$FCOENV.se[e.sel][1:k]<-FCOENV.c.se
    grmsem.out$model.out$FCOENV.se[e.sel][(k + 1):(2 * k)]<-FCOENV.s.se
    grmsem.out$model.out$FCOENV[notfree]<-NA
    grmsem.out$model.out$FCOENV.se[notfree]<-NA
    grmsem.out$model.out$FCOENV.Z<-grmsem.out$model.out$FCOENV/grmsem.out$model.out$FCOENV.se
    grmsem.out$model.out$FCOENV.p<-2*stats::pnorm(-abs(grmsem.out$model.out$FCOENV.Z))
    grmsem.out$model.out<-as.data.frame(grmsem.out$model.out)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  } else { #"IPC")
    #Derived parameters - vectors
    #Factorial co-heritabilities estimates
    FCOHER.c <- NULL
    FCOHER.s <- NULL
    FCOHER.c.se <- NULL
    FCOHER.s.se <- NULL
    FCOHER.c.se.formula <- NULL
    FCOHER.s.se.formula <- NULL
    AA.c.se.formula <- NULL
    AA.s.se.formula <- NULL
    AA.c.se <- NULL
    AA.s.se <- NULL
    A.c <- length(k)
    A.c <- estimates[a.sel][1:k]
    A.s <- diag(k)
    diag(A.s) <- estimates[a.sel][(k + 1):(2 * k)]
    VA <- A.c %*% t(A.c) + A.s %*% t(A.s)
    AA.c <- A.c*A.c
    AA.s <- A.s*A.s
    FCOHER.c <- AA.c/diag(VA)
    FCOHER.s <- AA.s/diag(VA)
    #Factorial co-environmentalities estimates
    FCOENV <- matrix(nrow = k, ncol = k)
    E <- diag(k)
    E[lower.tri(E, diag = TRUE)] <- estimates[e.sel]
    VE <- E %*% t(E)
    EE <- E*E
    FCOENV <- EE/diag(VE)
    #SE estimation
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
    VA.se.formula <- as.data.frame(VA.se.formula)
    VA.se.formula.c <- VA.se.formula.s <- VA.se.formula
    for (i in 1:k) {
      for (j in 1:k) {
        VA.se.formula.c[i, j] <- paste(XA.c[i], t(XA.c[j]), sep = " * ")
        VA.se.formula.s[i, j] <- 0
        if (i == j) {
          VA.se.formula.s[i, j] <- paste(XA.s[i, j], XA.s[i, j], sep = " * ")
        }
        VA.se.formula[i, j] <- paste(VA.se.formula.c[i, j], VA.se.formula.s[i, j], sep = " + ")
      }
    }
    for (i in 1:k) {
      AA.c.se.formula[i] <- paste("msm::deltamethod(~", paste("(", VA.se.formula.c[i,i],")",sep=""), ", estimates[free], VCOV)", sep = "")
      AA.s.se.formula[i] <- paste("msm::deltamethod(~", paste("(", VA.se.formula.s[i,i],")",sep=""), ", estimates[free], VCOV)", sep = "")
      AA.c.se[i] <- eval(parse(text = AA.c.se.formula[i]))
      AA.s.se[i] <- eval(parse(text = AA.s.se.formula[i]))
      FCOHER.c.se.formula[i] <-   paste("msm::deltamethod(~", paste("(", VA.se.formula.c[i,i],")/(", VA.se.formula[i,i], ")",sep=""), ", estimates[free], VCOV)", sep = "")
      FCOHER.c.se[i] <- eval(parse(text = FCOHER.c.se.formula[i]))
      FCOHER.s.se.formula[i] <-   paste("msm::deltamethod(~", paste("(", VA.se.formula.s[i,i],")/(", VA.se.formula[i,i], ")",sep=""), ", estimates[free], VCOV)", sep = "")
      FCOHER.s.se[i] <- eval(parse(text = FCOHER.s.se.formula[i]))
    }
    # Factorial coenvironmentality
    XE <- diag(k)
    diag(XE) <- 0
    XE[lower.tri(XE, diag = TRUE)] <- full.xlabel[e.sel]
    FCOENV.se.formula <- matrix(nrow = k, ncol = k)
    FCOENV.se.formula <- as.data.frame(FCOENV.se.formula)
    EE.se <- EE.se.formula <- EE.se.tmp.formula <- FCOENV.se <- FCOENV.se.tmp.formula  <- FCOENV.se.formula
    for (i in 1:k) {
      for (j in 1:k) {
        EE.se.tmp.formula[i, j] <- (paste(XE[i,j ], XE[i, j], sep = " * "))
      }
    }
    for (i in 1:k) {
      for (j in 1:k) {
        FCOENV.se.tmp.formula[i, j] <- (paste(EE.se.tmp.formula[i, ], collapse = " + "))
      }
    }
    for (i in 1:k) {
      for (j in 1:k) {
        EE.se.formula[i, j] <- paste("msm::deltamethod(~", EE.se.tmp.formula[i, j], ", estimates[free], VCOV)", sep = "")
        EE.se[i, j] <- eval(parse(text = EE.se.formula[i, j]))
        FCOENV.se.formula[i, j] <- paste("msm::deltamethod(~", paste("(", EE.se.tmp.formula[i, j], ")/(", FCOENV.se.tmp.formula[i, j], ")", sep=""), ", estimates[free], VCOV)", sep = "")
        FCOENV.se[i, j] <- eval(parse(text = FCOENV.se.formula[i, j]))
      }
    }
    #Updated model.out dataframe accounting for missing parameters
    grmsem.out$model.out$Vi <- NA
    grmsem.out$model.out$Vi[a.sel][1:k]<-AA.c
    grmsem.out$model.out$Vi[a.sel][(k + 1):(2 * k)]<-diag(AA.s)
    grmsem.out$model.out$Vi[e.sel]<-EE[lower.tri(EE, diag = TRUE)]
    grmsem.out$model.out$Vi.se <- NA
    grmsem.out$model.out$Vi.se[a.sel][1:k]<-AA.c.se
    grmsem.out$model.out$Vi.se[a.sel][(k + 1):(2 * k)]<-AA.s.se
    grmsem.out$model.out$Vi.se[e.sel]<-EE.se[lower.tri(EE.se, diag = TRUE)]
    grmsem.out$model.out$FCOHER <- NA
    grmsem.out$model.out$FCOHER.se <- NA
    grmsem.out$model.out$FCOHER[a.sel][1:k] <- FCOHER.c
    grmsem.out$model.out$FCOHER[a.sel][(k + 1):(2 * k)] <- diag(FCOHER.s)
    grmsem.out$model.out$FCOHER.se[a.sel][1:k]<-FCOHER.c.se
    grmsem.out$model.out$FCOHER.se[a.sel][(k + 1):(2 * k)]<-FCOHER.s.se
    grmsem.out$model.out$FCOHER[notfree]<-NA
    grmsem.out$model.out$FCOHER.se[notfree]<-NA
    grmsem.out$model.out$FCOHER.Z<-grmsem.out$model.out$FCOHER/grmsem.out$model.out$FCOHER.se
    grmsem.out$model.out$FCOHER.p<-2*stats::pnorm(-abs(grmsem.out$model.out$FCOHER.Z))
    grmsem.out$model.out$FCOENV <- NA
    grmsem.out$model.out$FCOENV.se <- NA
    grmsem.out$model.out$FCOENV[e.sel] <- FCOENV[lower.tri(FCOENV, diag = TRUE)]
    grmsem.out$model.out$FCOENV.se[e.sel]<-FCOENV.se[lower.tri(FCOENV.se, diag = TRUE)]
    grmsem.out$model.out$FCOENV[notfree]<-NA
    grmsem.out$model.out$FCOENV.se[notfree]<-NA
    grmsem.out$model.out$FCOENV.Z<-grmsem.out$model.out$FCOENV/grmsem.out$model.out$FCOENV.se
    grmsem.out$model.out$FCOENV.p<-2*stats::pnorm(-abs(grmsem.out$model.out$FCOENV.Z))
  }
  out <- list(fcoher.model.out = grmsem.out$model.out)
  class(out) <- ("grmsem.fcoher")
  return(out)
}

