#' grm import function
#'
#' This function imports genetic relationship matrices (GRMs) in gz format as e.g.
#' stored by the GCTA software command \code{gcta64 --grm test --make-grm-gz --out test}
#' @param file the name of the gz compressed grm file. No default.
#' @keywords grmsem
#' @importFrom utils read.table
#' @export
#' @return \code{grmsem.input} imports a zipped GCTA GRM (.gz) file and transforms it into a symmetric matrix
grm.input <- function(file) {
  if (!file.exists(file)) {
    stop(paste("File ", file, " not found."))
  }
  t <- read.table(file, header = FALSE)
  if (ncol(t) != 4) {
    stop("Table has to contain 4 columns")
  }
  ids <- NULL
  filebase <- gsub(".gz", "", file)
  idfile <- paste0(filebase, ".id")
  if (file.exists(idfile)) {
    ids <- unique(read.table(idfile)$V2)
  }
  return(vector2grm(t$V4, ids))
}

vector2grm <- function(v, ids = NULL) {
  n <- (sqrt(1 + 8 * length(v)) - 1) / 2 # number individuals
  A <- matrix(0, nrow = n, ncol = n)
  #if (!is.null(ids)) {
  #  rownames(A) <- ids
  #  colnames(A) <- ids
  #}
  A[upper.tri(A, diag = TRUE)] <- v
  A <- A + t(A)
  diag(A) <- .5 * diag(A)
  return(A)
}
