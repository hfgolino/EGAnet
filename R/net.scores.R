#' Network Scores
#'
#' @description This function computes network scores computed based on
#' each node's \code{\link[NetworkToolbox]{strength}} within each
#' community (i.e., factor) in the network (see \code{\link[EGAnet]{net.loads}}).
#' These values are used as network "factor loadings" for the weights of each item.
#' Notably, network analysis allows nodes to contribution to more than one community.
#' These loadings are considered in the network scores. In addition,
#' if the construct is a hierarchy (e.g., personality questionnaire;
#' items in facet scales in a trait domain), then an overall
#' score can be computed (see argument \code{global}). An important difference
#' is that the network scores account for cross-loadings in their
#' estimation of scores
#'
#' @param data Matrix or data frame.
#' Must be a dataset
#'
#' @param A Matrix, data frame, or \code{\link[EGAnet]{EGA}} object.
#' An adjacency matrix of network data
#'
#' @param wc Numeric.
#' A vector of community assignments.
#' Not necessary if an \code{\link[EGAnet]{EGA}} object
#' is input for argument \code{A}
#'
#' @param global Boolean.
#' Should general network loadings be computed in scores?
#' Defaults to \code{FALSE}.
#' If there is more than one dimension and there is theoretically
#' one global dimension, then general loadings of the dimensions
#' onto the global dimension can be included in the weighted
#' scores
#'
#' @param impute Character.
#' In the presence of missing data, imputation can be implemented. Currently,
#' three options are available:
#'
#' \itemize{
#'
#' \item{\strong{\code{none}}}
#' {No imputation is performed. This is the default.}
#'
#' \item{\strong{\code{mean}}}
#' {The "mean" value of the columns are used to replace the missing data.}
#'
#' \item{\strong{\code{median}}}
#' {The "median" value of the columns are used to replace the missing data.}
#'
#' }
#'
#' @param ... Additional arguments for \code{\link[EGAnet]{EGA}}
#'
#' @return Returns a list containing:
#'
#' \item{unstd.scores}{The unstandardized network scores for each participant
#' and community (including the overall score)}
#'
#' \item{std.scores}{The standardized network scores for each participant
#' and community (including the overall score)}
#'
#' \item{commCor}{Partial correlations between the specified or identified communities}
#'
#' \item{loads}{Standardized network loadings for each item in each dimension
#' (computed using \code{\link[EGAnet]{net.loads}})}
#'
#' @details For more details, type \code{vignette("Network_Scores")}
#'
#' @examples
#' # Load data
#' wmt <- wmt2[,7:24]
#'
#' \dontrun{
#'  # Estimate EGA
#'  ega.wmt <- EGA(wmt)
#' }
#'
#' # Network scores
#' net.scores(data = wmt, A = ega.wmt)
#'
#' @references
#' Christensen, A. P., & Golino, H. (2021).
#' On the equivalency of factor and network loadings.
#' \emph{Behavior Research Methods}.
#' \doi{10.3758/s13428-020-01500-6}
#'
#' Christensen, A. P., Golino, H., & Silvia, P. J. (in press).
#' A psychometric network perspective on the measurement and assessment of personality traits.
#' \emph{European Journal of Personality}.
#' \doi{10.1002/per.2265}
#'
#' Golino, H., Christensen, A. P., Moulder, R., Kim, S., & Boker, S. M. (under review).
#' Modeling latent topics in social media using Dynamic Exploratory Graph Analysis: The case of the right-wing and left-wing trolls in the 2016 US elections.
#' \emph{PsyArXiv}.
#' \doi{10.31234/osf.io/tfs7c}
#'
#' @author Alexander P. Christensen <alexpaulchristensen@gmail.com> and Hudson F. Golino <hfg9s at virginia.edu>
#'
#' @export
#'
#Network Scores
#Updated: 27.12.2020
net.scores <- function (data, A, wc, global = FALSE, impute, ...)
{
  if (missing(data)) {
    stop("Argument 'data' is required for analysis")
  }
  if (any(class(A) == "EGA")) {
    wc <- A$wc
    A <- A$network
  }
  else if (missing(A)) {
    stop("Adjacency matrix is required for analysis")
  }
  if (missing(impute)) {
    impute <- "none"
    warning("Argument 'impute' is missing. No imputation will be used.")
  }
  else if (missing(wc)) {
    wc <- rep(1, ncol(data))
  }
  P <- net.loads(A = A, wc = wc, pos.manifold = TRUE)$std
  nfacts <- length(unique(wc))
  if (nfacts > 1) {
    if (global) {
      fact.res <- as.data.frame(matrix(0, nrow = nrow(data),
                                       ncol = (nfacts + 1)))
    }
    else {
      fact.res <- as.data.frame(matrix(0, nrow = nrow(data),
                                       ncol = nfacts))
    }
  }
  else {
    fact.res <- as.data.frame(matrix(0, nrow = nrow(data),
                                     ncol = nfacts))
  }

  missing <- rowSums(is.na(data))
  if (impute != "none") {
    data <- data.matrix(data)
    miss <- which(is.na(data), arr.ind = TRUE)
    if (impute == "mean") {
      item.means <- colMeans(data, na.rm = TRUE)
      data[miss] <- item.means[miss[, 2]]
    }
    else {
      item.med <- apply(data, 2, median, na.rm = TRUE)
      data[miss] <- item.med[miss[, 2]]
    }
  }

  net.score.fxn <- function(loads, data) {
    net.sco <- matrix(0, nrow = nrow(data), ncol = ncol(loads))
    for (i in 1:ncol(loads)) {
      f.load <- loads[which(loads[, i] != 0), i]
      names(f.load) <- row.names(loads)[which(loads[, i] !=
                                                0)]
      dat <- data[, names(f.load)]
      f.sds <- apply(dat, 2, sd, na.rm = TRUE)
      rel <- f.load/f.sds
      rel.wei <- rel/sum(rel)
      net.sco[, i] <- as.vector(rowSums(t(t(dat) * rel.wei), na.rm = TRUE))
    }
    colnames(net.sco) <- colnames(loads)
    return(net.sco)
  }
  net.sco <- net.score.fxn(P, data)
  fact.res[, 1:nfacts] <- net.sco
  if (nfacts > 1) {
    colnames(fact.res)[1:nfacts] <- colnames(P)
  }
  else {
    colnames(fact.res) <- "1"
  }
  res <- list()

  if (nfacts > 1) {
    if (global) {
      ega.gen <- suppressMessages(EGA(net.sco, plot.EGA = FALSE,
                                      model = "glasso", ...))
      C <- ega.gen$network
      res$commCor <- C
      nl.gen <- net.loads(A = ega.gen$network, wc = rep(1,
                                                        ncol(net.sco)), min.load = 0, pos.manifold = TRUE)$std
      sds.gen <- apply(net.sco, 2, sd, na.rm = TRUE)
      rel.gen <- nl.gen/sds.gen
      rel.wei.gen <- (rel.gen/sum(rel.gen)) + 1
      G <- as.vector(rowSums(t(t(net.sco) * as.vector(as.matrix(rel.wei.gen)))))
      fact.res[, (nfacts + 1)] <- G
      colnames(fact.res)[nfacts + 1] <- "Overall"
    }
  }
  res$unstd.scores <- as.data.frame(round(fact.res, 3))
  res$std.scores <- as.data.frame(round(apply(fact.res, 2,
                                              scale), 3))
  res$loads <- P
  
  # Class
  class(res) <- "NetScores"
  
  return(res)
}
#----
