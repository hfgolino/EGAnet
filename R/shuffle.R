#'  Estimating the number of dimensions for n datasets with shuffled variables.
#'
#' \code{shuffle} Apply EGA to n datasets with m shuffled variables. The number of
#' variables shuffled is defined as the square root of the total number of variables in the data.
#'
#' @param data A data.frame object.
#' @param n Numeric. Number of estimates.
#' @param ncores Numeric. Number of cores to use in parallel.
#' @author Hudson F. Golino <hfg9s at virginia.edu>
#' @examples
#' \dontrun{
#' wmt.surrogate <- shuffle(data = wmt2[,7:24], n = 100, ncores = 4)
#' }
#' @seealso \code{\link{subsamples}} to estimate the number of dimensions via EGA in n random subsamples of the original data and \code{\link{surrogate}}
#' to apply a surrogate method for EGA.
#' @export

## Surrogate EGA:
shuffle <- function(data, n, ncores){
  require(compiler)
  require(foreach)
  require(doParallel)
  vars <- sqrt(ncol(data))
  N <- ncol(data)
  sample.vars <- vector("list", n)
  sample.data <- vector("list", n)
  ega.surrugate <- vector("list", n)
  sample.data <- vector("list", n)
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  vector
  all.ega.surrugate <- foreach(i = 1:n, .combine=rbind) %dopar% {
    sample.vars[[i]] <- sort(sample(seq_len(N), vars))
    sample.data[[i]] <- data[, sample.vars[[i]], drop = FALSE]
    sample.data[[i]] <- as.data.frame(apply(sample.data[[i]], 2, sample))
    sample.data[[i]] <- cbind(data[,-(sample.vars[[i]])], sample.data[[i]])
    ega.surrugate[[i]] <- ndim(sample.data[[i]])
  }
  stopCluster(cl)
  rownames(all.ega.surrugate) <- NULL
  return(all.ega.surrugate)
}

