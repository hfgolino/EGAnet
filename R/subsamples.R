#'  Applying EGA to n random subsamples.
#'
#' \code{subsamples} Applies EGA to n random subsamples of the original data.
#'
#' @param data A data.frame object.
#' @param n Numeric. Number of random subsamples
#' @param ncores Numeric. Number of cores to use in parallel.
#' @author Hudson F. Golino <hfg9s at virginia.edu>
#' @examples
#' \dontrun{
#' wmt.subsample <- subsamples(data = wmt2[,7:24], n = 100, ncores = 4)
#' }
#' \dontrun{
#' EGA.subsamples()
#' }
#' @seealso \code{\link{EGA}} to estimate the number of dimensions of an instrument using EGA, \code{\link{CFA}} to
#' verify the fit of the structure suggested by EGA using confirmatory factor analysis, \code{\link{shuffle}} to generate n
#' estimations of the number of dimensions in surrogate samples of the original dataset.
#' @export

## EGA in multiple subsamples:
subsamples <- function(data, n, ncores){
  require(compiler)
  require(foreach)
  require(doParallel)
  sample.data <- list(n)
  ega.subsamples <- list(n)
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  vector
  all.ega.subsamples <-
    foreach(i = 1:n, .combine=rbind) %dopar% {
    sample.data[[i]] <- data[sample(1:nrow(data), nrow(data)/2, replace=FALSE),]
    ega.subsamples[[i]] <- ndim(sample.data[[i]])
    }
  stopCluster(cl)
  rownames(all.ega.subsamples) <- NULL
  return(all.ega.subsamples)
}
