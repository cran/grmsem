#' grmsem model fitting function
#'
#' This function fits a grmsem model.
#' 
#' @details grmsem models estimate genetic (A) and residual (E) variance/covariance of quantitative traits (AE model), where E in GRM-based methods can capture both, shared and unique residual influences. 
#' The estimation of parameters and their SEs is performed with the function \code{grmsem.fit}. 
#' Specifically, the loglikelihood is estimated with \code{stats::optim} and the BFGS (Broyden-Fletcher-Goldfarb-Shannon) approach and the variance/covariance matrix of estimated parameters 
#' with \code{numDeriv::hessian}. The statistical significance of estimated parameters is assessed using a Wald test, assuming multivariate normality. 
#' 
#' \code{grmsem.fit} allows fitting different models describing the underlying multivariate genetic architecture of quantitative traits, 
#' as captured by a genetic-relationship-matrix (GRM), using structural equation modelling (SEM) techniques and a maximum likelihood approach. 
#' The user can fit multiple predefined model structures to the data. A Cholesky decomposition, Independent Pathway, 
#' and hybrid Independent Pathway/Cholesky models can be fitted by setting the \code{model} option to \code{Cholesky}, \code{IP} or \code{IPC}, respectively. 
#' In addition, the Cholesky model can be re-parametrised as Direct Symmetric model, 
#' estimating genetic and residual covariances directly, using the \code{model} option \code{DS}. Each model can be adapted by the user by setting free parameters (\code{A.v.free} and \code{E.v.free} options) and starting values (\code{A.v.startval} and  \code{E.v.startval} 
#' options). 
#' 
#' Input parameters are returned as \code{model.in} list object. Output from the maximum likelihood estimation is also given as 
#' list \code{model.fit} and the fitted grmsem model with estimated parameters and SEs is returned as dataframe \code{model.out}. 
#' The returned \code{grmsem.fit} object can be used to estimate genetic and residual covariance and correlations (\code{grmsem.var} function), bivariate heritabilities (\code{grmsem.biher} function), 
#' and factorial co-heritabilities and co-environmentalities (\code{grmsem.fcoher} function). All estimated parameters of the fitted grmsem model can also be standardised using the function \code{grmsem.stpar}.
#' 
#' Listwise complete observations can be selected with the option \code{compl.ph}=\code{TRUE}. Otherwise, \code{grmsem.fit} fits, like GREML, all available data to the model 
#' with the default option \code{compl.ph} \code{FALSE}. Using the option \code{LogL}=\code{FALSE}, the user can check the model input parameters and formula without a 
#' maximum likelihood estimation. Using the option \code{estSE}=\code{FALSE}, the user can carry out a maximum likelihood estimation without the estimation of SEs for 
#' estimated parameters that require calculating the Hessian. Note that \code{grmsem.fit} should preferably be run in parallel, by setting the \code{cores} option to the required number of cores.
#' 
#' When \code{grmsem.fit} is called with \code{LogL} \code{TRUE}, the user will see also the iterations of the \code{stats::optim} loglikelihood estimation, which are not included in the exported grmsem.fit list object.
#' 
#' @param ph phenotype file as R dataframe, even for single vectors (columns: k phenotypes, rows: ni individuals in same order as G). No default.
#' @param G GRM matrix as provided by the grm.input or grm.bin.input.R function. Use the same order of individuals as in ph. No default.
#' @param A.v.free vector of free parameters for genetic factor loadings (free:1, not-free:0). Default NULL, all parameters are estimated.
#' @param E.v.free vector of free parameters for residual factor loadings (free:1, not-free:0). Default NULL, all parameters are estimated.
#' @param A.v.startval vector of starting values for genetic factor loadings. Default NULL, all starting values are generated.
#' @param E.v.startval vector of starting values for residual factor loadings. Default NULL, all starting values are generated.
#' @param LogL estimation of the loglikelihood using optim BFGS (TRUE/FALSE). Default FALSE.
#' @param estSE estimation of standard errors by recalculating the Hessian matrix. Default FALSE.
#' @param cores number of cores for multi-threaded calculations (numeric). Default 1.
#' @param model grmsem model selection. Options: "Cholesky", "IP", "IPC", "DS". Default "Cholesky".
#' @param compl.ph listwise complete observations across all phenotypes (all NA are excluded). Default FALSE.
#' @param printest print output of the model.fit function including estimates (printest.txt) that can be used as starting values. Default FALSE.
#' @param cluster cluster type. Options: "PSOCK", "FORK". Default "PSOCK".
#' @param optim optimisation function from stats or optimParallel libraries. Options: "optim", "optimParallel". Default "optim".
#' @param verbose additional model fit information: (i) phenotype vector, (ii) n of GRM and corresponding I matrix when data are missing, (iii) Hessian if \code{estSE} \code{TRUE}. Default FALSE.
#' @keywords grmsem
#' @export
#' @return \code{grmsem.fit} returns a grmsem.fit list object consisting of:
#' \item{model.in}{list of input parameters}
#' \item{formula}{matrix of the model specification}
#' \item{model.fit}{list output of the maximum likelihood estimation, if \code{LogL} \code{TRUE}}
#' \item{model.out}{dataframe of fitted grmsem model with estimated parameters and SEs, if \code{estSE} \code{TRUE}}
#' \item{VCOV}{variance/covariance matrix}
#' \item{k}{number of phenotypes}
#' \item{n}{total number of observations across all phenotypes}
#' \item{n.obs}{number of observations per phenotype}
#' \item{n.ind}{number of individuals with at least one phenotype}
#' \item{model}{type of grmsem model}
#' \item{ph.nms}{vector of phenotype names}
#' 
#' \code{model.in} list of input parameters:
#' \item{part}{a - genetic, e - residual parameters}
#' \item{label}{parameter label}
#' \item{value}{starting values}
#' \item{freepar}{free parameters}
#' 
#' \code{model.fit} list output of the maximum likelihood estimation:
#' \item{optimisation}{output via optim()}
#' \item{estimates}{estimated parameters: factor loadings for `Cholesky`, `IP` and `IPC` models, but variance components for `DS`}
#' \item{LL}{loglikelihood}
#' \item{calls}{optim() calls}
#' \item{convergence}{optim() convergence}
#' \item{message}{optim() message}
#' 
#' \code{model.out} data.frame of fitted grmsem model:
#' \item{label}{parameter label}
#' \item{estimates}{estimated parameters}
#' \item{gradient}{gradient}
#' \item{se}{SE}
#' \item{Z}{Z (Wald)}
#' \item{p}{p (Wald)}
#' @examples
#' #Set up a Cholesky model: Model formula and total number of parameters
#' #ph.small is a trivariate phenotype file for 100 individuals in the same order as G.small
#' #nrow = 100, ncol = 3
#' #(runtime should be less than one minute)
#' \donttest{out <- grmsem.fit(ph.small, G.small, LogL = FALSE, estSE = FALSE)}
#' 
#' #Run a Cholesky model:
#' \donttest{out <- grmsem.fit(ph.small, G.small, LogL = TRUE, estSE = TRUE)}
grmsem.fit <- function(ph, G, A.v.free = NULL, E.v.free = NULL, A.v.startval = NULL, E.v.startval = NULL, 
                     LogL = FALSE, estSE = FALSE, cores = 1, model = "Cholesky", 
                     compl.ph = FALSE, printest = FALSE, cluster = "PSOCK", optim = "optim", verbose = FALSE) {
  local
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Subfunctions
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Generate matrix of parameter labels for genetic factor loadings (Cholesky)
  grmsem.chol.a.par <- function(k) {
    #Number of genetic factor loadings in Cholesky model
    Ng <- k * (k + 1)/2
    #Create matrices of genetic factor loading labels
    A.label <- diag(k)
    diag(A.label) <- 0
    A.label <- as.data.frame(A.label)
    A.label.T <- A.label
    for (i in 1:k) {
      for (j in 1:i) {
        A.label[i, j] <- paste("a", i, j, sep = "")
        A.label.T[j, i] <- paste("a", i, j, sep = "")
      }
    }
    # Set full number of genetic path labels
    A.v.label <- NULL
    tmp <- diag(k)
    A.v.label <- A.label[lower.tri(tmp, diag = TRUE)]
    print("Vector of genetic factor loadings:")
    print(A.v.label)
    if (Ng != length(A.v.label)) {
      m.msg <- "The number of genetic parameters 'a' does not match the number of parameter labels"
      stop(m.msg)
    }
    A.formula <- matrix(nrow = k, ncol = k)
    #A.formula <- as.data.frame(A.formula)
    for (i in 1:k) {
      for (j in 1:k) {
        if (compl.ph){
          A.formula[i, j] <- paste("(", paste(paste(A.label[i, ], A.label.T[, j], sep = " * "), collapse = " + "), ") * G", sep = "")
        }else{
          A.formula[i, j] <- paste("(", paste(paste(A.label[i, ], A.label.T[, j], sep = " * "), collapse = " + "), ") * ", G.label [i, j], sep = "")
        }
      }
    }
    return(list(A.v.label = A.v.label, Ng = Ng, A.formula = A.formula))
  }
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Generate matrix of parameter labels for residual factor loadings (Cholesky)
  grmsem.chol.e.par <- function(k) {
    #Number of residual factor loadings in Cholesky model
    Ne <- k * (k + 1)/2
    # Create matrices of residual factor loading labels
    E.label <- diag(k)
    E.label <- as.data.frame(E.label)
    E.label.T <- E.label
    for (i in 1:k) {
      for (j in 1:i) {
        E.label[i, j] <- paste("e", i, j, sep = "")
        E.label.T[j, i] <- paste("e", i, j, sep = "")
      }
    }
    # Set full number of genetic path labels
    E.v.label <- NULL
    tmp <- diag(k)
    E.v.label <- E.label[lower.tri(tmp, diag = TRUE)]
    print("Vector of residual factor loadings:")
    print(E.v.label)
    if (Ne != length( E.v.label)) {
      m.msg <- "Number of residual parameters 'e' does not match the number of parameter labels"
      stop(m.msg)
    }
    # Generate matrix of parameter combinations
    E.formula <- matrix(nrow = k, ncol = k)
    #E.formula <- as.data.frame(E.formula)
    for (i in 1:k) {
      for (j in 1:k) {
        if (compl.ph){
          E.formula[i, j] <- paste("(", paste(paste(E.label[i, ], E.label.T[, j], sep = " * "), collapse = " + "), ") * I", sep = "")
        }else{
          E.formula[i, j] <- paste("(", paste(paste(E.label[i, ], E.label.T[, j], sep = " * "), collapse = " + "), ") * ", I.label [i, j], sep = "")
        }
      }
    }
    return(list(E.v.label = E.v.label, Ne = Ne, E.formula = E.formula))
  }
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Generate matrix of parameter labels for genetic factor loadings (Independent pathway model)
  grmsem.ind.a.par <- function(k) {
    #Number of genetic factor loadings
    Ng <- 2 * k
    # Common and specific genetic factors
    A.common <- paste("ac", seq(1:k), sep = "")
    A.specific <- paste("as", seq(1:k), sep = "")
    A.v.label <- c(A.common, A.specific)
    print("Vector of genetic factor loadings:")
    print(A.v.label)
    if (Ng != length(A.v.label)) {
      g.msg <- "Number of genetic parameters 'a' does not match the number of parameter labels"
      stop(g.msg)
    }
    A.formula.c <- matrix(nrow = k, ncol = k)
    #A.formula.c <- as.data.frame(A.formula.c)
    for (i in 1:k) {
      for (j in 1:k) {
        A.formula.c[i, j] <- paste(A.common[i], t(A.common[j]), sep = " * ")
      }
    }
    # Specific factors
    A.formula.s <- diag(k)
    #A.formula.s <- as.data.frame(A.formula.s)
    for (i in 1:k) {
      for (j in 1:k) {
        if (i == j) {
          A.formula.s[i, j] <- paste(A.specific[i], A.specific[i], sep = " * ")
        }
      }
    }
    # Final A formula
    A.formula <- diag(k)
    #A.formula <- as.data.frame(A.formula)
    for (i in 1:k) {
      for (j in 1:k) {
        if(compl.ph){
          A.formula[i, j] <- paste("(", paste(A.formula.c[i, j], A.formula.s[i, j], sep = " + "), ") * G", sep = "")
        }else{
          A.formula[i, j] <- paste("(", paste(A.formula.c[i, j], A.formula.s[i, j], sep = " + "), ") * ", G.label [i, j], sep = "")
        }
      }
    }
    #print(A.formula)
    return(list(A.v.label = A.v.label, Ng = Ng, A.formula = A.formula))
  }
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Generate matrix of parameter labels for residual factor loadings (Independent pathway model)
  grmsem.ind.e.par <- function(k) {
    # Number of residual factor loadings
    Ne <- 2 * k
    # Common and specific residual factors
    E.common <- paste("ec", seq(1:k), sep = "")
    E.specific <- paste("es", seq(1:k), sep = "")
    E.v.label <- c(E.common, E.specific)
    print("Vector of residual factor loadings:")
    print(E.v.label)
    if (Ne != length(E.v.label)) {
      m.msg <- "Number of residual parameters 'e' does not match the number of parameter labels"
      stop(m.msg)
    }
    E.formula.c <- matrix(nrow = k, ncol = k)
    #E.formula.c <- as.data.frame(E.formula.c)
    for (i in 1:k) {
      for (j in 1:k) {
        E.formula.c[i, j] <- paste(E.common[i], t(E.common[j]), sep = " * ")
      }
    }
    # Specific factors
    E.formula.s <- diag(k)
    #E.formula.s <- as.data.frame(E.formula.s)
    for (i in 1:k) {
      for (j in 1:k) {
        if (i == j) {
          E.formula.s[i, j] <- paste(E.specific[i], E.specific[i], sep = " * ")
        }
      }
    }
    # Final E formula
    E.formula <- diag(k)
    #E.formula <- as.data.frame(E.formula)
    for (i in 1:k) {
      for (j in 1:k) {
        if(compl.ph){
          E.formula[i, j] <- paste("(", paste(E.formula.c[i, j], E.formula.s[i, j], sep = " + "), ") * I", sep = "")
        }else{
          E.formula[i, j] <- paste("(", paste(E.formula.c[i, j], E.formula.s[i, j], sep = " + "), ") * ", I.label [i, j], sep = "")
        }
      }
    }
    #print(E.formula)
    return(list(E.v.label = E.v.label, Ne = Ne, E.formula = E.formula))
  }
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Generate matrix of parameter labels for genetic variances (Direct symmetric model)
  grmsem.ds.a.par <- function(k) {
    #Number of genetic factor loadings in DS model
    Ng <- k * (k + 1)/2
    #Create matrices of genetic factor loading labels
    A.label <- diag(k)
    diag(A.label) <- 0
    A.label <- as.data.frame(A.label)
    for (i in 1:k) {
      for (j in 1:i) {
        A.label[i, j] <- paste("A", i, j, sep = "")
        A.label[j, i] <- paste("A", i, j, sep = "")
      }
    }
    # Set full number of genetic path labels
    A.v.label <- NULL
    tmp <- diag(k)
    A.v.label <- A.label[lower.tri(tmp, diag = TRUE)]
    print("Vector of genetic variances:")
    print(A.v.label)
    if (Ng != length(A.v.label)) {
      m.msg <- "The number of genetic variances 'A' does not match the number of parameter labels"
      stop(m.msg)
    }
    A.formula <- matrix(nrow = k, ncol = k)
    #A.formula <- as.data.frame(A.formula)
    for (i in 1:k) {
      for (j in 1:k) {
        if (compl.ph){
          A.formula[i, j] <- paste("(", A.label[i,j], ") * G", sep = "")
        }else{
          A.formula[i, j] <- paste("(", A.label[i,j], ") * ", G.label [i, j], sep = "")
        }
      }
    }
    return(list(A.v.label = A.v.label, Ng = Ng, A.formula = A.formula))
  }
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Generate matrix of parameter labels for residual factor loadings (Direct symmetric model)
  grmsem.ds.e.par <- function(k) {
    #Number of residual factor loadings in DS model
    Ne <- k * (k + 1)/2
    # Create matrices of residual factor loading labels
    E.label <- diag(k)
    E.label <- as.data.frame(E.label)
    for (i in 1:k) {
      for (j in 1:i) {
        E.label[i, j] <- paste("E", i, j, sep = "")
        E.label[j, i] <- paste("E", i, j, sep = "")
      }
    }
    # Set full number of genetic path labels
    E.v.label <- NULL
    tmp <- diag(k)
    E.v.label <- E.label[lower.tri(tmp, diag = TRUE)]
    print("Vector of residual factor loadings:")
    print(E.v.label)
    if (Ne != length( E.v.label)) {
      m.msg <- "Number of residual variances 'E' does not match the number of parameter labels"
      stop(m.msg)
    }
    # Generate matrix of parameter combinations
    E.formula <- matrix(nrow = k, ncol = k)
    #E.formula <- as.data.frame(E.formula)
    for (i in 1:k) {
      for (j in 1:k) {
        if (compl.ph){
          E.formula[i, j] <- paste("(", E.label[i, j], ") * I", sep = "")
        }else{
          E.formula[i, j] <- paste("(", E.label[i, j], ") * ", I.label [i, j], sep = "")
        }
      }
    }
    return(list(E.v.label = E.v.label, Ne = Ne, E.formula = E.formula))
  }
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Loglikelihood estimation
  grmsem.ll <- function(x) {
    estimates[free] <- x
    names(estimates) <- label
    # mnemonic assignments
    for(i in seq_along(estimates)){
      assign(names(estimates)[i], estimates[i], envir = as.environment(-1), inherits=FALSE)
      #assign(names(estimates)[i], estimates[i], inherits=TRUE)
      
    }
    # Set up expected var/covar matrix
    V <- matrix(0,  n,  n)
    for (i in 1:k) {
      for (j in i:k) { #recommended by A. Klassmann
        x.block.start <- start.obs[i]
        x.block.stop <- stop.obs[i]
        y.block.start <- start.obs[j]
        y.block.stop <- stop.obs[j]
        #V[x.block.start:x.block.stop, y.block.start:y.block.stop] <- as.formula(all.formula[i, j]) #too long
        V[x.block.start:x.block.stop, y.block.start:y.block.stop] <- eval(parse(text = all.formula[i, j]))
        #       if (i != j) {
        #          V[y.block.start:y.block.stop, x.block.start:x.block.stop] <- t(V[x.block.start:x.block.stop, y.block.start:y.block.stop])
        #        }
      }
    }
    # log likelihood
    cholV <- try(chol(V), silent = TRUE)
    if (!is.matrix(cholV)) return(NA)
    logSqrtDetV <- sum(log(diag(cholV)))
    temp <- crossprod(pvec, backsolve(cholV, diag(length(pvec))))
    q <- tcrossprod(temp)
    return(-logSqrtDetV - 0.5 * q)
  }
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Estimation of standard errors
  # For BFGS, the update of the Hessian at convergence of the parameters might not yet be a good approximation of the true Hessian,
  # therefore we calculate the Hession using the following steps:
  grmsem.se <- function(model.estimates, free, label) {
    #Full length vectors
    estimates <- rep(NA, length(label))
    gradient <- rep(NA, length(label))
    se <- rep(NA, length(label))
    Z <- rep(NA, length(label))
    pZ <- rep(NA, length(label))
    #Model-specific vector of length free
    model.se <- rep(NA, length(free))
    #Gradient evaluated at the model parameter estimates
    model.gradient <- numDeriv::grad(grmsem.ll, model.estimates)
    #Hessian matrix evaluated at the model parameter estimates
    model.hess <- numDeriv::hessian(grmsem.ll, model.estimates)
    if (verbose) {
      print ("Hessian")
      print(model.hess)
    }
    # Estimated variance-covariance matrix of the parameter estimates
    VCOV <- try(solve(-model.hess), silent = TRUE)
    if (is.matrix(VCOV)) {
      model.se <- sqrt(diag(VCOV))
      estimates[free] <- model.estimates
      gradient[free] <- model.gradient
      se[free] <- model.se
      Z[free] <- model.estimates/model.se
      pZ[free] <- 2 * stats::pnorm(-abs(model.estimates/model.se))
      colnames(VCOV) <- label[free]
      rownames(VCOV) <- label[free]
    }
    model.out <- data.frame(label = label, estimates = estimates, gradient = gradient, se = se, Z = Z, p = pZ)
    return(list(model.out = model.out, VCOV = VCOV))
  }
  ################################################################################
  #Main body
  ################################################################################
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Check input parameters
  ph.nms <- colnames(ph)
  ph <- as.matrix(ph)
  if (!is.numeric(ph)) {
    stop('ph is not numeric.')
  }
  if (!is.numeric(G)) {
    stop('G is not numeric.')
  }
  if (!is.numeric(cores)) {
    stop('The number of cores has to be numeric.')
  }
  #if (any(is.na(ph))) {
  #  stop('ph contains missing data.')
  #}
  if (any(is.na(G))) {
    stop('G contains missing data.')
  }
  model.vector <- c("Cholesky", "IP","IPC","DS")
  if (!is.element(model, model.vector)) {
    m.msg <- paste("The ", model,  " model is not specified.", sep = "")
    stop(m.msg)
  }
  cluster.vector <- c("PSOCK", "FORK")
  if (!is.element( cluster,  cluster.vector)) {
    m.msg <- paste("The ", cluster,  " cluster is not specified.", sep = "")
    stop(m.msg)
  }
  optim.vector <- c("optim", "optimParallel")
  if (!is.element( optim,  optim.vector)) {
    m.msg <- paste("The ", optim,  " optimisation function is not specified.", sep = "")
    stop(m.msg)
  }
  # Number of phenotypes or repeated assessments
  k <- ncol(ph)
  # Number of individuals within the phenotype data file  
  nph <- nrow(ph)
  if (nph < k) {
    m.msg <- paste("Model misspecification: The number of phenotypes (", k, ") is larger than the number of individuals per phenotype (", nph , ").", sep = "")
    stop(m.msg)
  }
  if (nph != nrow(G)) {
    m.msg <- paste("The number of phenotype rows (", nrow(ph), ") does not match G (", nrow(G), ").", sep = "")
    stop(m.msg)
  }
  if (k == 1 & (model == "IP" | model == "IPC" )) {
    m.msg <- paste("Single phenotypes should be analysed using a Cholesky model or direct symmetric (DS) model.", sep = "")
    stop(m.msg)
  }
  # Implement complete case analysis option
  if (compl.ph) {
    missing <- stats::complete.cases(ph)
    G <- G[missing,missing]
    ph <- as.data.frame(ph[missing,])
    I<-diag(nrow(G))
  }
  # Generate vector of missing phenotypes
  miss.ph <- matrix(ncol = ncol(ph), nrow = nrow(ph))
  for (i in 1:k){
    miss.ph[, i] <- !is.na(ph[, i])
  }
  # N individuals with at least one phenotype
  row.has.all.na <- apply(ph, 1, function(x){all(is.na(x))})
  n.ind <- dim(ph[!row.has.all.na,])[1]
  #Determine phenotypic missingness and create phenotype vector
  pvec<-NULL
  n.obs<-length(k)
  start.obs<-length(k)
  stop.obs<-length(k)
  count<-0
  for (i in 1:k){
    pvec <- as.vector(c(pvec,ph[, i][miss.ph[, i]]))
    n.obs[i] <- length(ph[, i][miss.ph[, i]])
    count <- count + 1
    start.obs[i] <- count
    count <- count + n.obs[i]-1
    stop.obs[i]<-count
  }
  # Total number of observations
  n<-length(pvec)
  if (verbose) {
    print("Total number of observations")
    print(n)
    print("Number of observations per phenotype")
    print(n.obs)
    print("Phenotype vector")
    print(start.obs)
    print(stop.obs)
  }
  # Select observations for incomplete phenotypes
  if (!compl.ph) {
    # G label
    G.label <- diag(k)
    G.label <- as.data.frame(G.label)
    for (i in 1:k) {
      for (j in 1:i) {
        G.label [i, j] <- paste("G", i, j, sep = "")
        G.label [j, i] <- paste("G", j, i, sep = "")
      }
    }
    # I label
    I.label <- diag(k)
    I.label <- as.data.frame(I.label)
    for (i in 1:k) {
      for (j in 1:i) {
        I.label [i, j] <- paste("I", i, j, sep = "")
        I.label [j, i] <- paste("I", j, i, sep = "")
      }
    }
    # Create G matrices
    for (i in 1:k) {
      for (j in 1:i) {
        Gtmp <- G[ miss.ph[,i], miss.ph[,j]]
        #assign(G.label[i,j], Gtmp, inherits=TRUE)
        assign(G.label[i,j], Gtmp, envir = as.environment(-1), inherits=FALSE)
        #assign(G.label[j,i], t(Gtmp), inherits=TRUE)
        assign(G.label[j,i], t(Gtmp), envir = as.environment(-1), inherits=FALSE)
        if (verbose) {
          print("GRM with missing data")
          print(c(i,j))
          print(dim(Gtmp))
          print(c(j,i))
          print(dim(t(Gtmp)))
        }
      }
    }
    # Create I matrices
    I<-diag(nrow(G))
    for (i in 1:k) {
      for (j in 1:i) {
        Itmp <- I[ miss.ph[,i], miss.ph[,j]]
        #assign(I.label[i,j], Itmp, inherits=TRUE)
        assign(I.label[i,j], Itmp, envir = as.environment(-1), inherits=FALSE)
        #(I.label[j,i], t(Itmp), inherits=TRUE)
        assign(I.label[j,i], t(Itmp), envir = as.environment(-1), inherits=FALSE)
        if (verbose) {
          print("I with missing data")
          print(c(i,j))
          print(dim(Itmp))
          print(c(j,i))
          print(dim(t(Itmp)))
        }
      }
    }
  }
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Model selection
  a.par <- NULL
  e.par <- NULL
  if (model == "Cholesky") {
    print("Parameters for a Cholesky model:")
    a.par <- grmsem.chol.a.par(k)
    e.par <- grmsem.chol.e.par(k)
  }
  else if (model == "IP") {
    print("Parameters for an Independent Pathway model:")
    a.par <- grmsem.ind.a.par(k)
    e.par <- grmsem.ind.e.par(k)
  }
  else if (model == "IPC") {
    print("Parameters for a Mixed Independent Pathway model (genetic part) and Cholesky model (residual part):")
    a.par <- grmsem.ind.a.par(k)
    e.par <- grmsem.chol.e.par(k)
  }
  else if (model == "DS") {
    print("Parameters for a direct symmetric model (variance components):")
    a.par <- grmsem.ds.a.par(k)
    e.par <- grmsem.ds.e.par(k)
  }
  A.v.label <- a.par$A.v.label
  A.formula <- a.par$A.formula
  Ng <- a.par$Ng
  E.v.label <- e.par$E.v.label
  Ne <- e.par$Ne
  E.formula <- e.par$E.formula
  all.formula <- matrix(nrow = k, ncol = k)
  #all.formula <- as.data.frame(all.formula)
  for (i in 1:k) {
    for (j in 1:k) {
      all.formula[i, j] <- paste(A.formula[i, j], E.formula[i, j], sep = " + ")
    }
  }
  print("Model formula")
  print(all.formula)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Generate free parameter vectors if required, and check whether the existing parameters are numeric(0, 1)
  if (is.null(A.v.free)) {
    A.v.free <- rep(1, Ng)
  }
  if ( any(A.v.free != 1 & A.v.free != 0)) {
    m.msg <- paste("Free parameters in", A.v.free, " must be given using 0 and 1.", sep = "")
    stop(m.msg)
  }
  if (is.null(E.v.free)) {
    E.v.free <- rep(1, Ne)
  }
  if ( any(E.v.free != 1 & E.v.free != 0)) {
    m.msg <- paste("Free parameters in", E.v.free, " must be given using 0 and 1.", sep = "")
    stop(m.msg)
  }
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Generate starting values if required, and check whether the existing parameters are numeric
  if (model == "Cholesky") {
    if (is.null(A.v.startval)) {
      ph.var <- stats::var(ph, na.rm=TRUE)
      lambda <- matrix(data=NA, nrow=k,ncol=k)
      lambda <- try(t(chol(ph.var)))
      A.v.startval <- lambda[lower.tri(lambda, diag = TRUE)]
    }
    if (is.null(E.v.startval)) {
      ph.var <- stats::var(ph, na.rm=TRUE)
      lambda <- matrix(data=NA, nrow=k,ncol=k)
      lambda <- try(t(chol(ph.var)))
      E.v.startval <- lambda[lower.tri(lambda, diag = TRUE)]
    }
    if (!is.numeric(A.v.startval)) {
      stop('Automatic A.v.startval is not numeric, please supply starting values.')
    }
    if (!is.numeric(E.v.startval)) {
      stop('Automatic E.v.startval is not numeric, please supply starting values')
    }
  } else  if (model == "IPC") {
    if (is.null(A.v.startval)) {
      #A.v.startval <- stats::runif(Ng, min = -1, max = 1)
      #Factor solution for 1 factor
      fact.out<-try(stats::factanal(stats::na.omit(ph), facto=1))
      #common factors loadings: fact.out$loadings
      #specific factor loadings:  sqrt(fact.out$uniquenesses)
      fact.start <- c(fact.out$loading,sqrt(fact.out$uniquenesses))
      names(fact.start)<-NULL
      A.v.startval <- fact.start
    }
    if (is.null(E.v.startval)) {
      ph.var <- stats::var(ph, na.rm=TRUE)
      lambda <- matrix(data=NA, nrow=k,ncol=k)
      lambda <- try(t(chol(ph.var)))
      E.v.startval <- lambda[lower.tri(lambda, diag = TRUE)]
    }
    if (!is.numeric(A.v.startval)) {
      stop('Automatic A.v.startval is not numeric, please supply starting values.')
    }
    if (!is.numeric(E.v.startval)) {
      stop('Automatic E.v.startval is not numeric, please supply starting values.')
    }
  } else if (model == "IP") {
    if (is.null(A.v.startval)) {
      #A.v.startval <- stats::runif(Ng, min = -1, max = 1)
      #Factor solution for 1 factor
      fact.out<-try(stats::factanal(stats::na.omit(ph), facto=1))
      #common factors loadings: fact.out$loadings
      #specific factor loadings:  sqrt(fact.out$uniquenesses)
      fact.start <- c(fact.out$loading,sqrt(fact.out$uniquenesses))
      names(fact.start)<-NULL
      A.v.startval<-fact.start
    }
    if (!is.numeric(A.v.startval)) {
      stop('Automatic A.v.startval is not numeric, please supply starting values.')
    }
    if (is.null(E.v.startval)) {
      #A.v.startval <- stats::runif(Ng, min = -1, max = 1)
      #Factor solution for 1 factor
      fact.out<-try(stats::factanal(stats::na.omit(ph), facto=1))
      #common factors loadings: fact.out$loadings
      #specific factor loadings:  sqrt(fact.out$uniquenesses)
      fact.start <- c(fact.out$loading,sqrt(fact.out$uniquenesses))
      names(fact.start)<-NULL
      E.v.startval <- fact.start
    }
    if (!is.numeric(E.v.startval)) {
      stop('Automatic E.v.startval is not numeric, please supply starting values.')
    }
  }else if (model == "DS") {
    if (is.null(A.v.startval)) {
      ph.var<-stats::var(ph, na.rm=TRUE)
      A.v.startval <- ph.var[lower.tri(ph.var, diag = TRUE)]
    }
    if (is.null(E.v.startval)) {
      ph.var<-stats::var(ph, na.rm=TRUE)
      E.v.startval <- ph.var[lower.tri(ph.var, diag = TRUE)]
    }
    if (!is.numeric(A.v.startval)) {
      stop('Automatic A.v.startval is not numeric, please supply starting values.')
    }
    #if (is.null(E.v.startval)) {
    #  E.v.startval <- stats::runif(Ne, min = 0.1, max = 0.9)
    #}
    if (!is.numeric(E.v.startval)) {
      stop('Automatic E.v.startval is not numeric, please supply starting values.')
    }
  } 
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Check whether the vector for parameter values, free values and starting values match each other
  if (length(A.v.label) != length(A.v.free) || length(A.v.label) != length(A.v.startval)) {
    m.msg <- paste("The number of values for genetic factor loadings (", length(A.v.label), "), free parameters (", length(A.v.free), ") and starting values (", length(A.v.startval), ") does not match!", sep = "")
    stop(m.msg)
  }
  if (length(E.v.label) != length(E.v.free) || length(E.v.label) != length(E.v.startval)) {
    m.msg <- paste("The number of values for residual factor loadings (", length(E.v.label), "), free parameters (", length(E.v.free), ") and starting values (", length(E.v.startval), ") does not match!", sep = "")
    stop(m.msg)
  }
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Model parameter summary
  part <- c(rep("a", length(A.v.label)), rep("e", length(E.v.label)))
  label <- c(A.v.label, E.v.label)
  value <- c(A.v.startval, E.v.startval)
  freepar <- c(A.v.free, E.v.free)
  model.in <- data.frame(part = part, label = label, value = value, freepar = freepar)
  free <- which(model.in$freepar == 1)
  x <-value[free]
  estimates <- value
  print("Input parameters")
  print(model.in) # EDIT
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Fit loglikelihood
  model.fit <- NULL
  model.out <- NULL
  VCOV <- NULL
  if (LogL) {
    print("Fitting the model...")
    cl <- parallel::makeCluster(cores, type = cluster)
    parallel::setDefaultCluster(cl = cl)
    if(compl.ph){
      parallel::clusterExport(cl, c("grmsem.ll","n","k", "G", "I", "pvec", "x", "estimates", "free", "label", "all.formula","start.obs", "stop.obs"), envir = environment())
    }else{
      export.cl<-c("grmsem.ll", "n", "k", unname(unlist(G.label)), unname(unlist(I.label)), "pvec", "x", "estimates", "free", "label", "all.formula", "start.obs", "stop.obs")
      parallel::clusterExport(cl, export.cl, envir = environment())
    }
    if (optim == "optim") {
      model.fit <- stats::optim(x, grmsem.ll, method = "BFGS", control = list(fnscale = -1, trace = 10, maxit = 100000, reltol = 1e-9, REPORT = 1))
    }else{
      model.fit <- optimParallel::optimParallel(x, grmsem.ll, method = "BFGS", control = list(fnscale = -1, trace = 10, maxit = 100000, reltol = 1e-9, REPORT = 1))
    }
    parallel::stopCluster(cl)
    print(model.fit)
    names(model.fit) <- c("estimates", "LL", "calls","convergence", "message")
    if (printest == TRUE) {
      sink("temp.estimates.txt") #modelfit estimates that can be used as starting values (system crashes)
      print(model.fit)
      sink()
    }
    if (estSE & all(!is.na(model.fit$estimates))) {
      grmsem.se.out <- grmsem.se(model.fit$estimates, free, label)
      model.out <- grmsem.se.out$model.out
      VCOV <- grmsem.se.out$VCOV
    } else {
      model.out <- NULL
      VCOV <- NULL
    }
    print("Done.")
    
  }
  out <- list(model.in = model.in, formula = all.formula, model.fit = model.fit, model.out = model.out, VCOV = VCOV, k = k, n = n, n.obs = n.obs, n.ind = n.ind, model = model, ph.nms = ph.nms)
  class(out) <- ("grmsem.fit")
  #rm(list=estimates)
  return(out)
}
