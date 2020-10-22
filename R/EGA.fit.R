#' \code{\link[EGAnet]{EGA}} Optimal Model Fit using the Total Entropy Fit Index  (\code{\link[EGAnet]{tefi}})
#'
#' @description Estimates the best fitting model using \code{\link[EGAnet]{EGA}}.
#' The number of steps in the \code{\link[igraph]{cluster_walktrap}} detection
#' algorithm is varied and unique community solutions are compared using
#' \code{\link[EGAnet]{tefi}}.
#'
#' @param data Matrix or data frame.
#' Dataset or correlation matrix
#' 
#' @param n Integer.
#' Sample size (if the data provided is a correlation matrix)
#'
#' @param model Character.
#' A string indicating the method to use.
#' Defaults to \code{"glasso"}
#'
#' Current options are:
#'
#' \itemize{
#'
#' \item{\strong{\code{"glasso"}}}
#' {Estimates the Gaussian graphical model using graphical LASSO with
#' extended Bayesian information criterion to select optimal regularization parameter.
#' See \code{\link[EGAnet]{EBICglasso.qgraph}}}
#'
#' \item{\strong{\code{"TMFG"}}}
#' {Estimates a Triangulated Maximally Filtered Graph.
#' See \code{\link[NetworkToolbox]{TMFG}}}
#'
#' }
#'
#' @param steps Numeric vector.
#' Range of steps to be used in the model selection.
#' Defaults from 3 to 8 steps (based on Pons & Latapy, 2006)
#'
#' @return Returns a list containing:
#'
#' \item{EGA}{The \code{\link[EGAnet]{EGA}} output for the best fitting model}
#'
#' \item{steps}{The number of steps used in the best fitting model from
#' the \code{\link[igraph]{cluster_walktrap}} algorithm}
#'
#' \item{EntropyFit}{The \code{\link[EGAnet]{tefi}} Index for the unique solutions given the range of steps
#' (vector names represent the number of steps)}
#'
#' \item{Lowest.EntropyFit}{The lowest value for the \code{\link[EGAnet]{tefi}} Index}
#'
#' @examples
#'
#' # Load data
#' wmt <- wmt2[,7:24]
#'
#' \donttest{
#' # Estimate normal EGAtmfg
#' tmfg <- EGA(data = wmt, model = "TMFG", plot.EGA = FALSE)
#'
#' # Estimate optimal EGAtmfg
#' tmfg.opt <- EGA.fit(data = wmt, model = "TMFG")
#'
#' # Compare with CFA
#' cfa.tmfg <- CFA(tmfg, estimator = "WLSMV", data = wmt)
#' cfa.opt <- CFA(tmfg.opt$EGA, estimator = "WLSMV", data = wmt)
#'
#' lavaan::lavTestLRT(cfa.tmfg$fit, cfa.opt$fit, method = "satorra.bentler.2001")
#' }
#'
#' @references
#' # Entropy fit measures \cr
#' Golino, H., Moulder, R. G., Shi, D., Christensen, A. P., Garrido, L. E., Neito, M. D., Nesselroade, J., Sadana, R., Thiyagarajan, J. A., & Boker, S. M. (in press).
#' Entropy fit indices: New fit measures for assessing the structure and dimensionality of multiple latent variables.
#' \emph{Multivariate Behavioral Research}.
#' doi: \href{https://doi.org/10.31234/osf.io/mtka2}{10.31234/osf.io/mtka2}
#' 
#' # Original implementation of EGA.fit \cr
#' Golino, H., Thiyagarajan, J. A., Sadana, M., Teles, M., Christensen, A. P., & Boker, S. M. (under review).
#' Investigating the broad domains of intrinsic capacity, functional ability, and environment:
#' An exploratory graph analysis approach for improving analytical methodologies for measuring healthy aging.
#' \emph{PsyArXiv}.
#' doi: \href{https://doi.org/10.31234/osf.io/hj5mc}{10.31234/osf.io/hj5mc}
#' 
#' # Walktrap algorithm \cr
#' Pons, P., & Latapy, M. (2006).
#' Computing communities in large networks using random walks.
#' \emph{Journal of Graph Algorithms and Applications}, \emph{10}, 191-218.
#' doi:\href{https://doi.org/10.7155/jgaa.00185}{10.7155/jgaa.00185}
#'
#' @seealso \code{\link[EGAnet]{bootEGA}} to investigate the stability of EGA's estimation via bootstrap,
#' \code{\link[EGAnet]{EGA}} to estimate the number of dimensions of an instrument using EGA,
#' and \code{\link[EGAnet]{CFA}} to verify the fit of the structure suggested by EGA using confirmatory factor analysis.
#'
#' @author Hudson Golino <hfg9s at virginia.edu> and Alexander P. Christensen <alexpaulchristensen@gmail.com>
#'
#' @export
# EGA fit
# Updated 20.10.2020
EGA.fit <- function (data, model = c("glasso","TMFG"),
                     steps = c(3,4,5,6,7,8), n = NULL)
{
  if(missing(model))
  {model <- "glasso"
  }else{model <- match.arg(model)}

  if(missing(steps))
  {steps <- c(3,4,5,6,7,8)
  }else{steps <- steps}

  best.fit <- list()

      num <- length(steps)

      mods <- list()
      dims <- matrix(NA, nrow = ncol(data), ncol = num)

      #Generate walktrap models
      for(i in 1:num)
      {
        message(paste("Estimating EGA -- Walktrap model",i,"of",num,sep=" "))
        mods[[as.character(steps[i])]] <- EGA(data = data,
                                              n = n,
                                              model = model,
                                              model.args = list(steps = steps[i]),
                                              algorithm = "walktrap",
                                              plot.EGA = FALSE)

        dims[,i] <- mods[[as.character(steps[i])]]$wc
      }

      colnames(dims) <- as.character(steps)
      
      #remove solutions with missing dimensions
      rm.cols <- which(apply(apply(dims, 2, is.na), 2, any))
      
      if(length(rm.cols) != 0)
      {
        dims <- dims[,-rm.cols]
        steps <- steps[-rm.cols]
      }

      #check for unique number of dimensions
      step <- as.numeric(colnames(dims)[which(!duplicated(homogenize.membership(dims[,1], dims), MARGIN = 2))])
      
      len <- length(step)

      #if all models are the same
      if(len==1)
      {
        best.fit$EGA <- mods[[1]]
        best.fit$steps <- 4
        message("All EGA models are identical.")
      }else{

        ent.vec <- vector("numeric",length=len)

        for(i in 1:len)
        {ent.vec[i] <- tefi(abs(mods[[as.character(step[i])]]$correlation), mods[[as.character(step[i])]]$wc)$VN.Entropy.Fit}

        names(ent.vec) <- step

        best.fit$EGA <- mods[as.character(step[which(ent.vec==min(ent.vec))])]
        best.fit$steps <- step[which(ent.vec==min(ent.vec))]
        best.fit$EntropyFit <- ent.vec
        best.fit$Lowest.EntropyFit <- ent.vec[which(ent.vec==min(ent.vec))]
      }
      
  # Get information for EGA Methods section
  args <- list()
  
  args$model <- model
  args$algorithm <- "walktrap"
  args$steps <- range(steps)
  args$entropy <- best.fit$Lowest.EntropyFit
  args$solutions <- best.fit$EntropyFit
  
  best.fit$Methods <- args
      
  class(best.fit) <- "EGA.fit"
      
  return(best.fit)
}
#----
