#' Conditional Mutual Information
#'
#' Computes the conditional mutual information metric using a modification of the matrix of partial correlations (see Zhao, Zhou,Zhang, & Chen, 2016).
#' If the raw data is provided, the correlation matrix will be computed using the \code{\link[qgraph]{cor_auto}} function of the \code{\link[qgraph]{qgraph}} package.
#'
#' @param data A dataframe with the variables to be used in the analysis or a correlation matrix.
#'
#' @param network Logical.
#' If TRUE, returns a plot of the conditional mutual information network.
#' Defaults to FALSE.
#'
#' @param EGA Logical.
#' If TRUE, exploratory graph analysis is performed using the conditional mutual information network.
#'
#' @param steps Number of steps to be used in \code{\link[igraph]{cluster_walktrap}} algorithm (necessary only if the EGA argument is set to TRUE).
#' Defaults to 4.
#'
#' @param scale Logical.
#' If TRUE, the conditional mutual information will be scaled using the following formula: trace(X)/p, where X is the matrix with the conditional
#' mutual information and p is the number of variables.
#'
#' @author Hudson F. Golino <hfg9s at virginia.edu>
#'
#' @examples
#'
#' \donttest{
#' #estimate EGA using CMI:
#'
#' cmi1 <- cmi(data = depression[,48:68], network = FALSE, EGA = TRUE, steps = 4, scale = TRUE)
#'
#' @seealso \code{\link{bootEGA}} to investigate the stability of EGA's estimation via bootstrap
#' and \code{\link{EGA}} to apply the exploratory graph analysis techinique.
#'
#' @references
#' Zhao, J., Zhou, Y., Zhang, X., & Chen, L. (2016).
#' Part mutual information for quantifying direct associations in networks.
#' \emph{Proceedings of the National Academy of Sciences}, \emph{113(18)}, 5130-5135.
#'
#' @importFrom Matrix nearPD
#' @importFrom stats na.omit cov2cor
#' @importFrom qgraph cor_auto
#'
#' @export
#'
#'
cmi <- function(data, network = FALSE, EGA = TRUE, steps = 4, scale = TRUE){
  if(nrow(data)!=ncol(data)){
    cor.data <- qgraph::cor_auto(data)
  } else{
    cor.data <- data
  }
  pcor <- solve(cor.data)
  pcor <- -stats::cov2cor(pcor)
  cmi <- -(1/2)*log10(1-pcor)

  if(scale == TRUE){
    cmi <- cmi/(sum(diag(cmi))/ncol(cmi))
    cmi <- -cmi
    diag(cmi) <- -diag(cmi)
  }

  if(network==TRUE){
    qgraph::qgraph(cmi, layout = "spring", vsize = 6)
  }

  if(EGA==TRUE){

    graph <- NetworkToolbox::convert2igraph(abs(cmi))
    wc <- igraph::walktrap.community(graph, steps = steps)
    names(wc$membership) <- colnames(data)
    n.dim <- max(wc$membership)
    wc <- wc$membership

    a <- list()
    # Returning only communities that have at least two items:
    if(length(unique(wc))>1){
      indices <- seq_along(wc)
      indices2 <- indices[wc %in% wc[duplicated(wc)]]
      wc[indices[-indices2]] <- NA
      a$n.dim <- length(unique(na.omit(wc)))
    }else{
      a$n.dim <- length(unique(wc))
    }
    a$cmi <- cmi
    a$wc <- wc
    dim.variables <- data.frame(items = colnames(data), dimension = a$wc)
    dim.variables <- dim.variables[order(dim.variables[, 2]),                                ]
    a$dim.variables <- dim.variables

    # Plot:
    plot.ega <- qgraph::qgraph(a$cmi, layout = "spring",
                                   vsize = 6, groups = as.factor(a$wc), label.prop = 1, legend = TRUE)
    } else{

    a <- list()
    a$cmi <- cmi

    }
  class(a) <- "CMI"
  return(a)
}
#----

