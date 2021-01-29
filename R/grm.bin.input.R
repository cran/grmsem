#' grm import function
#'
#' This function imports genetic relationship matrices (GRMs)
#' in binary format as e.g. stored by the GCTA software command \code{gcta64 --grm test --make-grm-bin --out test}
#' @param file name of the binary grm file. No default.
#' @param size byte size used (typically gcta uses 4). See eponymous parameter of function \code{\link{readBin}}.
#' @return object of type matrix
#' @export
#' @return \code{grmsem.bin.input} imports a binary GCTA GRM (.bin) file and transforms it into a symmetric matrix
grm.bin.input <-
  function(file, size = 4)
  {
    if (!file.exists(file)) {
      stop(paste("File ", file, " not found."))
    }
    ids <- NULL
    filebase <- gsub(".bin", "", file)
    idfile <- paste0(filebase, ".id")
    if (file.exists(idfile)) {
      ids <- unique(read.table(idfile)$V2)
    }
    vlength <- file.size(file) / size  # length of vector
    n <- (sqrt(1 + 8 * vlength) - 1) / 2      # number individuals
    if (as.integer(n) != n) {
      stop("Cannot determine number of entries in the grm.")
    }
    BinFile <-
      file(file, "rb") # matrix entries (serialized triangle)
    v <- readBin(BinFile,
                 n = n * (n + 1) / 2,
                 what = numeric(0),
                 size = size)
    close(BinFile)
    return(vector2grm(v, ids))
  }
