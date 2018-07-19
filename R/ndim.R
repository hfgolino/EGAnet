#'  Estimates the number of dimensions (for simulation studies, resampling or surrogate analysis).
#'
#' \code{ndim} Applies EGA to estimate the number of dimensions for simulation studies,
#' resampling or surrogate analysis.
#'
#' @param data A data.frame object.
#' @examples
#' wmt.dim <- ndim(data = wmt2[,7:24])
#'
#' \dontrun{
#' ndim(data = wmt2[,7:24])
#' }
#' @seealso \code{\link{EGA}} to estimate the number of dimensions of an instrument using EGA, \code{\link{CFA}} to
#' verify the fit of the structure suggested by EGA using confirmatory factor analysis, \code{\link{shuffle}} to generate n
#' estimations of the number of dimensions in shuffled versions of the original dataset and \code{\link{subsamples}} to
#' apply EGA to n random subsamples of the original data.
#' @export

## Estimating the numnber of latent dimensions:

ndim <- function (data) {
  require(qgraph)
  require(igraph)
  data <- as.data.frame(data)
  cor.data <- cor_auto(data)
  glasso.ebic <- EBICglasso(S = cor.data, n = nrow(data),
                            lambda.min.ratio = 0.1)
  graph.glasso <- as.igraph(qgraph(abs(glasso.ebic), DoNotPlot = TRUE), attributes = FALSE)
  wc <- walktrap.community(graph.glasso)
  n.dim <- max(wc$membership)
  return(n.dim)
}
