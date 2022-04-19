#' Hierarchical \code{\link[EGAnet]{EGA}}
#'
#' Estimates EGA using the lower-order solution of \code{\link[igraph]{cluster_louvain}}
#' to identify the lower-order dimensions and then uses factor or network loadings to estimate factor or
#' network scores to estimate the higher-order dimensions
#'
#' @param data Matrix or data frame.
#' Variables (down columns) only.
#' Does not accept correlation matrices
#' 
#' @param scores Character.
#' How should scores for the lower-order structure be estimated?
#' Defaults to \code{"network"} for network scores computed using
#' the \code{\link[EGAnet]{net.scores}} function.
#' Set to \code{"factor"} for factor scores computed using
#' \code{\link[psych]{fa}}. Factors are assumed to be correlated
#' using the \code{"oblimin"} rotation
#' 
#' @param uni.method Character.
#' What unidimensionality method should be used? 
#' Defaults to \code{"LE"}.
#' Current options are:
#' 
#' \itemize{
#'
#' \item{\strong{\code{expand}}}
#' {Expands the correlation matrix with four variables correlated .50.
#' If number of dimension returns 2 or less in check, then the data 
#' are unidimensional; otherwise, regular EGA with no matrix
#' expansion is used. This is the method used in the Golino et al. (2020)
#' \emph{Psychological Methods} simulation.}
#'
#' \item{\strong{\code{LE}}}
#' {Applies the leading eigenvalue algorithm (\code{\link[igraph]{cluster_leading_eigen}})
#' on the empirical correlation matrix. If the number of dimensions is 1,
#' then the leading eigenvalue solution is used; otherwise, regular EGA
#' is used. This is the final method used in the Christensen, Garrido,
#' and Golino (2021) simulation.}
#' 
#' }
#' 
#' @param corr Type of correlation matrix to compute. The default uses \code{\link[qgraph]{cor_auto}}.
#' Current options are:
#'
#' \itemize{
#'
#' \item{\strong{\code{cor_auto}}}
#' {Computes the correlation matrix using the \code{\link[qgraph]{cor_auto}} function from
#' \code{\link[qgraph]{qgraph}}}.
#'
#' \item{\strong{\code{pearson}}}
#' {Computes Pearson's correlation coefficient using the pairwise complete observations via
#' the \code{\link[stats]{cor}}} function.
#'
#' \item{\strong{\code{spearman}}}
#' {Computes Spearman's correlation coefficient using the pairwise complete observations via
#' the \code{\link[stats]{cor}}} function.
#' }
#'
#' @param model Character.
#' A string indicating the method to use.
#' Defaults to \code{"glasso"}.
#' Current options are:
#'
#' \itemize{
#'
#' \item{\strong{\code{glasso}}}
#' {Estimates the Gaussian graphical model using graphical LASSO with
#' extended Bayesian information criterion to select optimal regularization parameter}
#'
#' \item{\strong{\code{TMFG}}}
#' {Estimates a Triangulated Maximally Filtered Graph}
#'
#' }
#'
#' @param model.args List.
#' A list of additional arguments for \code{\link[EGAnet]{EBICglasso.qgraph}}
#' or \code{\link[NetworkToolbox]{TMFG}}
#'
#' @param algorithm A string indicating the algorithm to use or a function from \code{\link{igraph}}
#' Defaults to \code{"walktrap"}.
#' Current options are:
#'
#' \itemize{
#'
#' \item{\strong{\code{walktrap}}}
#' {Computes the Walktrap algorithm using \code{\link[igraph]{cluster_walktrap}}}
#'
#' \item{\strong{\code{louvain}}}
#' {Computes the Louvain algorithm using \code{\link[igraph]{cluster_louvain}}}
#'
#' }
#'
#' @param algorithm.args List.
#' A list of additional arguments for \code{\link[igraph]{cluster_walktrap}}, \code{\link[igraph]{cluster_louvain}},
#' or some other community detection algorithm function (see examples)
#'
#' @param plot.EGA Boolean.
#' If \code{TRUE}, returns a plot of the network and its estimated dimensions.
#' Defaults to \code{TRUE}
#'
#' @param plot.type Character.
#' Plot system to use.
#' Current options are \code{\link[qgraph]{qgraph}} and \code{\link{GGally}}.
#' Defaults to \code{"GGally"}
#'
#' @param plot.args List.
#' A list of additional arguments for the network plot.
#' For \code{plot.type = "qgraph"}:
#'
#' \itemize{
#'
#' \item{\strong{\code{vsize}}}
#' {Size of the nodes. Defaults to 6.}
#'
#'}
#' For \code{plot.type = "GGally"} (see \code{\link[GGally]{ggnet2}} for
#' full list of arguments):
#'
#' \itemize{
#'
#' \item{\strong{\code{vsize}}}
#' {Size of the nodes. Defaults to 6.}
#'
#' \item{\strong{\code{label.size}}}
#' {Size of the labels. Defaults to 5.}
#'
#' \item{\strong{\code{alpha}}}
#' {The level of transparency of the nodes, which might be a single value or a vector of values. Defaults to 0.7.}
#'
#' \item{\strong{\code{edge.alpha}}}
#' {The level of transparency of the edges, which might be a single value or a vector of values. Defaults to 0.4.}
#'
#'  \item{\strong{\code{legend.names}}}
#' {A vector with names for each dimension}
#'
#' \item{\strong{\code{color.palette}}}
#' {The color palette for the nodes. For custom colors,
#' enter HEX codes for each dimension in a vector.
#' See \code{\link[EGAnet]{color_palette_EGA}} for
#' more details and examples}
#'
#' }
#' 
#' @param verbose Boolean.
#' Should network estimation parameters be printed?
#' Defaults to \code{TRUE}.
#' Set to \code{FALSE} for no print out
#' 
#' @return Returns a list containing:
#' 
#' \item{lower}{Lower-order results from \code{\link[EGAnet]{EGA}}}
#' 
#' \item{lower_scores}{Estimated scores for the lower dimensions}
#' 
#' \item{lower_loadings}{Estaimted loadings for the lower dimensions}
#' 
#' \item{score_type}{Type of scores estimated and used in the higher-order \code{\link[EGAnet]{EGA}}}
#' 
#' \item{higher}{Higher-order results from \code{\link[EGAnet]{EGA}}}
#'
#' @author Luis E. Garrido <garrido.luiseduardo@gmail.com>,
#' Alexander P. Christensen <alexpaulchristensen@gmail.com>, and
#' Hudson F. Golino <hfg9s at virginia.edu>
#'
#' @examples
#' # Obtain example data
#' data <- optimism
#' 
#' \dontrun{
#' # hierEGA example
#' opt.res <- hierEGA(data = optimism, algorithm = "louvain")
#' }
#'
#'
#' @export
#' 
# Hierarchical EGA
# Updated 19.04.2022
hierEGA <- function(
    data, scores = c("factor", "network"),
    uni.method = c("expand", "LE"),
    corr = c("cor_auto", "pearson", "spearman"),
    model = c("glasso", "TMFG"), model.args = list(),
    algorithm = c("walktrap", "louvain"), algorithm.args = list(),
    plot.EGA = TRUE, plot.type = c("GGally", "qgraph"),
    plot.args = list(), verbose = TRUE
)
{
 
  #### ARGUMENTS HANDLING ####
  
  if(missing(scores)){
    scores <- "network"
  }else{scores <- match.arg(scores)}
  
  if(missing(uni.method)){
    uni.method <- "LE"
  }else{uni.method <- match.arg(uni.method)}
  
  if(missing(corr)){
    corr <- "cor_auto"
  }else{corr <- match.arg(corr)}
  
  if(missing(model)){
    model <- "glasso"
  }else{model <- match.arg(model)}
  
  if(missing(algorithm)){
    algorithm <- "walktrap"
  }else if(!is.function(algorithm)){
    algorithm <- tolower(match.arg(algorithm))
  }
  
  if(missing(plot.type)){
    plot.type <- "GGally"
  }else{plot.type <- match.arg(plot.type)}
  
  #### ARGUMENTS HANDLING ####
  
  # Ensure data is a matrix
  data <- as.matrix(data)
  
  # Check for symmetric matrix
  if(isSymmetric(data)){
    stop(
      paste(
        "Symmetric matrix detected. Data with variables down the columns and cases across the rows are required to compute",
        scores,
        "scores for this analysis."
      )
    )
  }
  
  # Obtain EGA defaults
  ega_defaults <- formals(EGA)
  
  # Remove "..."
  ega_defaults <- ega_defaults[-length(ega_defaults)]
  
  # Start lower-order ----
  
  # Set lower-order defaults
  lower_order_defaults <- ega_defaults
  lower_order_defaults$data <- data
  lower_order_defaults$uni.method <- uni.method
  lower_order_defaults$corr <- corr
  lower_order_defaults$model <- model
  lower_order_defaults$model.args <- model.args
  lower_order_defaults$algorithm <- "louvain" # for lower order communities
  lower_order_defaults$algorithm.args <- algorithm.args
  lower_order_defaults$plot.lower_order <- FALSE # do not plot
  lower_order_defaults$plot.type <- plot.type
  lower_order_defaults$plot.args <- plot.args
  lower_order_defaults$verbose <- FALSE # quiet
  lower_order_defaults$lower.louvain <- TRUE # provides lower order Louvain
  
  # Get EGA
  lower_order_result <- suppressMessages(
    do.call(
      EGA.estimate, lower_order_defaults
    )
  )
  
  # Get S3 print information
  if(is.null(colnames(data))){
    dim.variables <- data.frame(
      items = paste("V", 1:ncol(data), sep = ""),
      dimension = lower_order_result$wc
    )
  }else{
    dim.variables <- data.frame(
      items = colnames(data), 
      dimension = lower_order_result$wc
    )
  }
  dim.variables <- dim.variables[order(dim.variables[, 2]),]
  lower_order_result$dim.variables <- dim.variables
  
  # Set class for lower-order
  class(lower_order_result) <- "EGA"
  
  # Estimate scores
  if(scores == "factor"){
    
    # Determine appropriate correlations
    ## Categories
    categories <- apply(data, 2, function(x){
      length(na.omit(unique(x)))
    })
    
    ## Check categories
    if(any(categories > 2)){
      
      if(any(categories <= 2) & any(categories >= 3)){ # Mixed
        psych_cor <- "mixed"
      }else if(any(categories < 6)){ # Polychoric
        psych_cor <- "poly"
      }else{
        psych_cor <- "cor"
      }
      
    }else{# Tetrachoric
      psych_cor <- "tet"
    }
    
    # Estimate factor model
    fm <- psych::fa(
      r = data, # needs data for scores
      nfactors = length(na.omit(unique(lower_order_result$wc))), # number of factors
      rotate = "oblimin",
      cor = psych_cor
    )
    
    # Score estimates
    score_est <- fm$scores
    
    # Lower-order loadings
    lower_loads <- fm$loadings[,1:length(na.omit(unique(lower_order_result$wc)))]
    
  }else if(scores == "network"){
    
    # Compute network scores
    nt <- suppressWarnings(
      net.scores(
        data = data,
        A = lower_order_result$network,
        wc = lower_order_result$wc
      )
    )
    
    # Score estimates
    score_est <- nt$std.scores
    
    # Lower-order loadings
    lower_loads <- nt$loads
  
  }
  
  # End lower-order ----
  
  # Start higher-order ----
  
  # Set the rest of the arguments
  ega_defaults$data <- score_est # set data as score estimates
  ega_defaults$uni.method <- uni.method
  ega_defaults$corr <- corr
  ega_defaults$model <- model
  ega_defaults$model.args <- model.args
  ega_defaults$algorithm <- algorithm
  ega_defaults$algorithm.args <- algorithm.args
  ega_defaults$plot.EGA <- plot.EGA
  ega_defaults$plot.type <- plot.type
  ega_defaults$plot.args <- plot.args
  ega_defaults$verbose <- verbose
  
  # Get EGA
  ega_result <- do.call(
    EGA, ega_defaults
  )
  
  # Set up plots
  if(isTRUE(plot.EGA)){
    
    # Set up output
    ggpubr::ggarrange(
      plot(lower_order_result, produce = FALSE), # plot lower-order
      plot(ega_result, produce = FALSE), # plot higher-order
      labels = c("Lower-order", "Higher-order")
    )
    
    
    
  }
  
  # Return results
  results <- list()
  results$lower <- lower_order_result
  results$lower_scores <- score_est
  results$lower_loadings <- lower_loads
  results$score_type <- scores
  results$higher <- ega_result
  
  # Make class "hierEGA"
  class(results) <- "hierEGA"
  
  return(results)
}
