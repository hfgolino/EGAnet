#' Network Scores
#'
#' @description This function computes network scores computed based on
#' each node's \code{strength} within each
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
#' @param rotation Character.
#' A rotation to use, like factor loadings, to obtain
#' a simple structure. For a list of rotations,
#' see \link{GPArotation}
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
#' @param na.rm Boolean.
#' Should missing data be removed from network score computation?
#' Defaults to \code{TRUE}
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
#' # Estimate EGA
#' ega.wmt <- EGA(
#'   data = wmt,
#'   plot.EGA = FALSE # No plot for CRAN checks
#' )}
#'
#' # Network scores
#' net.scores(data = wmt, A = ega.wmt)
#' 
#' \dontrun{
#' # Produce Methods section
#' methods.section(
#'   ega.wmt,
#'   stats = "net.scores"
#' )}
#'
#' @references
#' Christensen, A. P., & Golino, H. (2021).
#' On the equivalency of factor and network loadings.
#' \emph{Behavior Research Methods}, \emph{53}, 1563-1580.
#' 
#' Christensen, A. P., Golino, H., & Silvia, P. J. (2020).
#' A psychometric network perspective on the validity and validation of personality trait questionnaires.
#' \emph{European Journal of Personality}, \emph{34}, 1095-1108.
#'
#' Golino, H., Christensen, A. P., Moulder, R., Kim, S., & Boker, S. M. (2021).
#' Modeling latent topics in social media using Dynamic Exploratory Graph Analysis: The case of the right-wing and left-wing trolls in the 2016 US elections.
#' \emph{Psychometrika}.
#'
#' @author Alexander P. Christensen <alexpaulchristensen@gmail.com> and Hudson F. Golino <hfg9s at virginia.edu>
#'
#' @export
#'
# Network Scores
# Updated: 07.01.2022
# Add rotation: 20.10.2022
net.scores <- function (
    data, A, wc, rotation = "oblimin",
    impute, na.rm = TRUE, ...
)
{
  
  # Depracate global argument
  if("global" %in% names(list(...))){
    message("Argument 'global' has been deprecated. Use `hierEGA` object instead.")
  }
  
  # Missing arguments handling
  if(missing(data)){
    stop("Argument 'data' is required for analysis")
  }
  
  # Set network and memberships
  if(is(A, "EGA")){
    wc <- A$wc
    A <- A$network
  }else if(is(A, "dynEGA")){
    wc <- A$dynEGA$wc
    A <- A$dynEGA$network
  }else if(missing(A)){
    stop("Adjacency matrix is required for analysis")
  }else if(missing(wc)){
    wc <- rep(1, ncol(data))
  }
  
  # Check for imputation
  if(missing(impute)){
    impute <- "none"
    warning("Argument 'impute' is missing. No imputation will be used.")
  }
  
  # Compute network loadings
  P <- net.loads(A = A, wc = wc, ...)
  nfacts <- length(unique(wc))
  fact.res <- as.data.frame(matrix(0, nrow = nrow(data),
                                   ncol = nfacts))

  missing <- rowSums(is.na(data))
  if (impute != "none") {
    data <- data.matrix(data)
    miss <- which(is.na(data), arr.ind = na.rm)
    if (impute == "mean") {
      item.means <- colMeans(data, na.rm = na.rm)
      data[miss] <- item.means[miss[, 2]]
    }else {
      item.med <- apply(data, 2, median, na.rm = na.rm)
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
      f.sds <- apply(dat, 2, sd, na.rm = na.rm)
      rel <- f.load/f.sds
      rel.wei <- rel/sum(rel)
      net.sco[, i] <- as.vector(rowSums(t(t(dat) * rel.wei), na.rm = na.rm))
    }
    colnames(net.sco) <- colnames(loads)
    return(net.sco)
  }
  net.sco <- net.score.fxn(P$std, data)
  fact.res[, 1:nfacts] <- net.sco
  if (nfacts > 1) {
    colnames(fact.res)[1:nfacts] <- colnames(P$std)
  }else{
    colnames(fact.res) <- "1"
  }
  
  net.sco.rotated <- net.score.fxn(P$rotated$loadings, data)
  fact.res.rotated <- fact.res
  fact.res.rotated[, 1:nfacts] <- net.sco.rotated
  if (nfacts > 1) {
    colnames(fact.res.rotated)[1:nfacts] <- colnames(P$rotated$loadings)
  }else{
    colnames(fact.res.rotated) <- "1"
  }
  
  res <- list(
    scores = list(
      unstd.scores = as.data.frame(round(fact.res, 3)),
      std.scores = as.data.frame(round(apply(fact.res, 2, scale), 3)),
      rot.scores = as.data.frame(round(apply(fact.res.rotated, 2, scale), 3))
    ),
    loadings = P
    
  )

  # Class
  class(res) <- "NetScores"
  
  return(res)
}
#----
