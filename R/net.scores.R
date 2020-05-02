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
#' Christensen, A. P., & Golino, H. (under review).
#' Statistical equivalency of factor and network loadings.
#' \emph{PsyArXiv}.
#' doi: \href{https://doi.org/10.31234/osf.io/xakez}{10.31234/osf.io/xakez}
#' 
#' Christensen, A. P., Golino, H. F., & Silvia, P. J. (in press).
#' A psychometric network perspective on the measurement and assessment of personality traits.
#' \emph{European Journal of Personality}.
#' doi: \href{https://doi.org/10.1002/per.2265}{10.1002/per.2265}
#' 
#' Golino, H., Christensen, A. P., Moulder, R., Kim, S., & Boker, S. M. (under review).
#' Modeling latent topics in social media using Dynamic Exploratory Graph Analysis: The case of the right-wing and left-wing trolls in the 2016 US elections.
#' \emph{PsyArXiv}.
#' doi: \href{https://doi.org/10.31234/osf.io/tfs7c}{10.31234/osf.io/tfs7c}
#'
#' @author Alexander P. Christensen <alexpaulchristensen@gmail.com> and Hudson F. Golino <hfg9s at virginia.edu>
#'
#' @export
#'
#Network Scores
#Updated: 02.05.2020
net.scores <- function (data, A, wc, global = FALSE, ...)
{
  ####Missing arguments checks####
  if(missing(data))
  {stop("Argument 'data' is required for analysis")}

  # Detect if input is an 'EGA' object
  if(class(A) == "EGA")
  {
    # Grab communities
    wc <- A$wc

    # Replace 'A' with 'EGA' network
    A <- A$network

  }else if(missing(A))
  {stop("Adjacency matrix is required for analysis")
  }else if(missing(wc)) #Default to single  variable
  {wc <- rep(1,ncol(data))}
  ####Missing arguments checks####
  
  #Compute network loadings
  P <- net.loads(A = A, wc = wc, pos.manifold = TRUE)$std

  #Number of factors
  nfacts <- length(unique(wc))

  #Initialize factor result matrix
  if(nfacts > 1)
  {
    if(global)
    {fact.res <- as.data.frame(matrix(0, nrow = nrow(data), ncol = (nfacts + 1)))
    }else{fact.res <- as.data.frame(matrix(0, nrow = nrow(data), ncol = nfacts))}
  }else{fact.res <- as.data.frame(matrix(0, nrow = nrow(data), ncol = nfacts))}

  ####NETWORK SCORE FUNCTION####
  net.score.fxn <- function(loads, data)
  {
    #Initialize participant  scores
    net.sco <- matrix(0, nrow = nrow(data), ncol = ncol(loads))

    #Compute  factor scores (ML)
    for(i in 1:ncol(loads))
    {
      #Network loadings for each factor
      f.load <- loads[which(loads[,i]!=0),i]
      names(f.load) <- row.names(loads)[which(loads[,i]!=0)]

      #Grab items associated with factor
      dat <- data[,names(f.load)]

      #Grab std dev of items associated with factor
      f.sds <- apply(dat,2,sd,na.rm = TRUE)

      #Obtain relative weights
      rel <- f.load / f.sds
      rel.wei <- rel / sum(rel)

      #Compute scores
      net.sco[,i] <- as.vector(rowSums(t(t(dat) * rel.wei)))
    }

    colnames(net.sco) <- colnames(loads)

    return(net.sco)
  }
  ####NETWORK SCORE FUNCTION####

  #Populate factor result matrix
  net.sco <- net.score.fxn(P, data)
  fact.res[,1:nfacts] <- net.sco

  if(nfacts > 1)
  {colnames(fact.res)[1:nfacts] <- colnames(P)
  }else{colnames(fact.res) <- "1"}

  #Initialize results list
  res <- list()

  #Global network loadings
  if(nfacts > 1)
  {
    if(global)
    {
      # Compute EGA
      ega.gen <- suppressMessages(EGA(net.sco, plot.EGA = FALSE, model = "glasso", ...))
      
      #Regularized correlations between factors
      C <- ega.gen$network
      
      #Return partial correlations between factors
      res$commCor <- C
      
      # Compute network loadings
      nl.gen <- net.loads(A = ega.gen$network, wc = rep(1, ncol(net.sco)), min.load = 0, pos.manifold = TRUE)$std
      
      #Grab std dev of items associated with factor
      sds.gen <- apply(net.sco , 2, sd, na.rm = TRUE)
      
      #Obtain relative weights
      rel.gen <- nl.gen / sds.gen
      rel.wei.gen <- (rel.gen / sum(rel.gen)) +1
      
      #Compute overall score
      G <- as.vector(rowSums(t(t(net.sco) * as.vector(as.matrix(rel.wei.gen)))))
      fact.res[,(nfacts + 1)] <- G
      colnames(fact.res)[nfacts + 1] <- "Overall"

    }
  }

  res$unstd.scores <- as.data.frame(round(fact.res,3))
  res$std.scores <- as.data.frame(round(apply(fact.res,2,scale),3))
  res$loads <- P

  return(res)
}
#----
