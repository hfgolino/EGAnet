#' Applies the Exploratory Graph Analysis technique
#'
#' Estimates the number of dimensions of a given dataset or correlation matrix
#' using the graphical lasso (\code{\link{EBICglasso.qgraph}}) or the
#' Triangulated Maximally Filtered Graph (\code{\link[EGAnet]{TMFG}})
#' network estimation methods.
#'
#' Two community detection algorithms, Walktrap (Pons & Latapy, 2006) and
#' Louvain (Blondel et al., 2008), are pre-programmed because of their
#' superior performance in simulation studies on psychological
#' data generated from factor models (Christensen & Golino; 2020; Golino et al., 2020).
#' Notably, any community detection algorithm from the \code{\link{igraph}}
#' can be used to estimate the number of communities (see examples).
#'
#' @param data Matrix or data frame.
#' Variables (down columns) or correlation matrix.
#' If the input is a correlation matrix,
#' then argument \code{n} (number of cases) is \strong{required}
#'
#' @param n Integer.
#' Sample size if \code{data} provided is a correlation matrix
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
#' @param uni.method Character.
#' What unidimensionality method should be used? 
#' Defaults to \code{"louvain"}.
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
#' {Applies the Leading Eigenvalue algorithm (\code{\link[igraph]{cluster_leading_eigen}})
#' on the empirical correlation matrix. If the number of dimensions is 1,
#' then the Leading Eigenvalue solution is used; otherwise, regular EGA
#' is used. This is the final method used in the Christensen, Garrido,
#' and Golino (2021) simulation.}
#' 
#' \item{\strong{\code{louvain}}}
#' {Applies the Louvain algorithm (\code{\link[igraph]{cluster_louvain}})
#' on the empirical correlation matrix using a resolution parameter = 0.95.
#' If the number of dimensions is 1, then the Louvain solution is used; otherwise,
#' regular EGA is used. This method was validated in the Christensen (2022) simulation.}
#' 
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
#' or \code{\link[EGAnet]{TMFG}}
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
#' \item{\strong{\code{leiden}}}
#' {Computes the Leiden algorithm using \code{\link[igraph]{cluster_leiden}}.
#' Defaults to \code{objective_function = "modularity"}}
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
#' @param consensus.iter Numeric.
#' Number of iterations to perform in consensus clustering for the Louvain algorithm
#' (see Lancichinetti & Fortunato, 2012).
#' Defaults to \code{100}
#' 
#' @param consensus.method Character.
#' What consensus clustering method should be used? 
#' Defaults to \code{"highest_modularity"}.
#' Current options are:
#' 
#' \itemize{
#' 
#' \item{\strong{\code{highest_modularity}}}
#' {Uses the community solution that achieves the highest modularity
#' across iterations}
#' 
#' \item{\strong{\code{most_common}}}
#' {Uses the community solution that is found the most
#' across iterations}
#' 
#' \item{\strong{\code{iterative}}}
#' {Identifies the most common community solutions across iterations
#' and determines how often nodes appear in the same community together.
#' A threshold of 0.30 is used to set low proportions to zero.
#' This process repeats iteratively until all nodes have a proportion of
#' 1 in the community solution.
#' }
#' 
#' \item{\code{lowest_tefi}}
#' {Uses the community solution that achieves the lowest \code{\link[EGAnet]{tefi}}
#' across iterations}
#' 
#' }
#'
#' @param plot.EGA Boolean.
#' If \code{TRUE}, returns a plot of the network and its estimated dimensions.
#' Defaults to \code{TRUE}
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
#' @param ... Additional arguments.
#' Used for deprecated arguments from previous versions of \code{\link{EGA}}
#'
#' @author Hudson Golino <hfg9s at virginia.edu>, Alexander P. Christensen <alexpaulchristensen at gmail.com>, Maria Dolores Nieto <acinodam at gmail.com> and Luis E. Garrido <garrido.luiseduardo at gmail.com>
#'
#' @return Returns a list containing:
#'
#' \item{network}{A symmetric network estimated using either the
#' \code{\link{EBICglasso.qgraph}} or \code{\link[EGAnet]{TMFG}}}
#'
#' \item{wc}{A vector representing the community (dimension) membership
#' of each node in the network. \code{NA} values mean that the node
#' was disconnected from the network}
#'
#' \item{n.dim}{A scalar of how many total dimensions were identified in the network}
#'
#' \item{cor.data}{The zero-order correlation matrix}
#'
#' @examples
#' # Obtain data
#' wmt <- wmt2[,7:24]
#' 
#' \dontrun{
#' # Estimate EGA
#' ega.wmt <- EGA(
#'   data = wmt,
#'   plot.EGA = FALSE # No plot for CRAN checks
#' )
#'
#' # Summary statistics
#' summary(ega.wmt)
#' 
#' # Produce Methods section
#' methods.section(ega.wmt)
#' 
#' # Estimate EGAtmfg
#' ega.wmt.tmfg <- EGA(
#'   data = wmt, model = "TMFG"
#' )
#'
#' # Estimate EGA with Louvain algorithm
#' ega.wmt.louvain <- EGA(
#'   data = wmt, algorithm = "louvain"
#' )
#' 
#' # Estimate EGA with Leiden algorithm
#' ega.wmt.leiden <- EGA(
#'   data = wmt, algorithm = "leiden"
#' )
#'
#' # Estimate EGA with Spinglass algorithm
#' ega.wmt.spinglass <- EGA(
#'   data = wmt,
#'   algorithm = igraph::cluster_spinglass
#' )}
#'
#' @seealso \code{\link{bootEGA}} to investigate the stability of EGA's estimation via bootstrap
#' and \code{\link{CFA}} to verify the fit of the structure suggested by EGA using confirmatory factor analysis.
#'
#' @references
#' # Louvain algorithm \cr
#' Blondel, V. D., Guillaume, J.-L., Lambiotte, R., & Lefebvre, E. (2008).
#' Fast unfolding of communities in large networks.
#' \emph{Journal of Statistical Mechanics: Theory and Experiment}, \emph{2008}, P10008.
#'
#' # Comprehensive unidimensionality simulation \cr
#' Christensen, A. P. (2022).
#' Unidimensional community detection: A Monte Carlo simulation, grid search, and comparison.
#' \emph{PsyArXiv}.
#'
#' # Compared all \emph{igraph} community detections algorithms, introduced Louvain algorithm, simulation with continuous and polytomous data \cr
#' # Also implements the Leading Eigenvalue unidimensional method \cr
#' Christensen, A. P., Garrido, L. E., & Golino, H. (2021).
#' Comparing community detection algorithms in psychological data: A Monte Carlo simulation.
#' \emph{PsyArXiv}.
#'
#' # Original simulation and implementation of EGA \cr
#' Golino, H. F., & Epskamp, S. (2017).
#' Exploratory graph analysis: A new approach for estimating the number of dimensions in psychological research.
#' \emph{PLoS ONE}, \emph{12}, e0174035.
#'
#' Golino, H. F., & Demetriou, A. (2017).
#' Estimating the dimensionality of intelligence like data using Exploratory Graph Analysis.
#' \emph{Intelligence}, \emph{62}, 54-70.
#'
#' # Current implementation of EGA, introduced unidimensional checks, continuous and dichotomous data \cr
#' Golino, H., Shi, D., Christensen, A. P., Garrido, L. E., Nieto, M. D., Sadana, R., & Thiyagarajan, J. A. (2020).
#' Investigating the performance of Exploratory Graph Analysis and traditional techniques to identify the number of latent factors: A simulation and tutorial.
#' \emph{Psychological Methods}, \emph{25}, 292-320.
#'
#' # Walktrap algorithm \cr
#' Pons, P., & Latapy, M. (2006).
#' Computing communities in large networks using random walks.
#' \emph{Journal of Graph Algorithms and Applications}, \emph{10}, 191-218.
#'
#' @importFrom stats cor rnorm runif na.omit
#'
#' @export
#'
# Updated 20.08.2022
# Louvain unidimensionality 27.07.2022
# Consensus clustering 13.05.2022
# LE adjustment 08.03.2021
EGA <- function (
    data, n = NULL,
    corr = c("cor_auto", "pearson", "spearman"),
    uni.method = c("expand", "LE", "louvain"),
    model = c("glasso", "TMFG"), model.args = list(),
    algorithm = c("walktrap", "leiden", "louvain"), algorithm.args = list(),
    consensus.method = c(
      "highest_modularity",
      "most_common",
      "iterative",
      "lowest_tefi"
    ), consensus.iter = 100, 
    plot.EGA = TRUE, plot.args = list(),
    ...
)
{

  # Get additional arguments
  add.args <- list(...)

  # Check if steps has been input as an argument
  if("steps" %in% names(add.args)){

    # Give deprecation warning
    warning(
      paste(
        "\nThe 'steps' argument has been deprecated in all EGA functions.\n\nInstead use: algorithm.args = list(steps = ", add.args$steps, ")\n",
        sep = ""
      )
    )

    # Handle the number of steps appropriately
    algorithm.args$steps <- add.args$steps
  }
  
  # Check if uni has been input as an argument
  if("uni" %in% names(add.args)){
    
    # Give deprecation warning
    warning(
      "\nThe 'uni' argument has been deprecated in all EGA functions.\n"
    )
  }

  #### ARGUMENTS HANDLING
  
  # if(missing(uni.method)){
  #   uni.method <- "LE"
  # }else{uni.method <- match.arg(uni.method)}
  # 
  # # Check if uni.method = "LE" has been used
  # if(uni.method == "LE"){
  #   # Give change warning
  #   warning(
  #     call. = FALSE,
  #     paste(
  #       "Previous versions of EGAnet (<= 0.9.8) checked unidimensionality using",
  #       styletext('uni.method = "expand"', defaults = "underline"),
  #       "as the default"
  #     )
  #   )
  # }else if(uni.method == "expand"){
  #   # Give change warning
  #   warning(
  #     call. = FALSE,
  #     paste(
  #       "Newer evidence suggests that",
  #       styletext('uni.method = "LE"', defaults = "underline"),
  #       'is more accurate than uni.method = "expand" (see Christensen, Garrido, & Golino, 2021 in references).',
  #       '\n\nIt\'s recommended to use uni.method = "LE"'
  #     )
  #   )
  # }
  
  if(missing(corr)){
    corr <- "cor_auto"
  }else{corr <- tolower(match.arg(corr))}
  
  if(missing(uni.method)){
    uni.method <- "louvain"
  }else{uni.method <- tolower(match.arg(uni.method))}

  if(missing(model)){
    model <- "glasso"
  }else{model <- match.arg(model)}

  if(missing(algorithm)){
    algorithm <- "walktrap"
  }else if(!is.function(algorithm)){
    algorithm <- tolower(match.arg(algorithm))
  }
  
  if(missing(consensus.method)){
    consensus.method <- "most_common"
  }else{consensus.method <- tolower(match.arg(consensus.method))}
  
  # Check for correlation matrix or data
  if(isSymmetric(unname(as.matrix(data)))){ ## Correlation matrix

    # Check for column names
    if(is.null(colnames(data))){
      colnames(data) <- paste("V", 1:ncol(data), sep = "")
    }
    
    # Force row names to be column names
    row.names(data) <- colnames(data)

    # Check for number of cases
    if(missing(n)){
      stop("There is no input for argument 'n'. Number of cases must be input when the matrix is square.")
    }
    
    # Set correlations as data
    correlation <- data
    
  }else{ # Data

    # Check for column names
    if(is.null(colnames(data))){
      colnames(data) <- paste("V", 1:ncol(data), sep = "")
    }

    # Get number of cases
    n <- nrow(data)
    
    # Compute correlation matrix
    correlation <- suppressMessages(
      switch(corr,
             "cor_auto" = qgraph::cor_auto(data, forcePD = TRUE),
             "pearson" = cor(data, use = "pairwise.complete.obs"),
             "spearman" = cor(data, method = "spearman", use = "pairwise.complete.obs")
      )
    )

  }
  
  # Unidimensional result
  uni.res <- try(
    unidimensionality.check(
      data = data, n = n, corr = corr,
      correlation = correlation, 
      uni.method = uni.method,
      model = model, model.args = model.args,
      algorithm = algorithm, algorithm.args = algorithm.args,
      consensus.method = consensus.method, consensus.iter = consensus.iter
    ),
    silent = TRUE
  )
  
  # Error check
  if(any(class(uni.res) == "try-error")){
    return(
      error.report(
        result = uni.res,
        SUB_FUN = "unidimensionality.check",
        FUN = "EGA"
      )
    )
  }

  # Multidimensional result
  multi.res <- try(
    EGA.estimate(
      data = correlation, n = n, corr = corr,
      model = model, model.args = model.args,
      algorithm = algorithm, algorithm.args = algorithm.args,
      consensus.method = consensus.method, consensus.iter = consensus.iter
    ),
    silent = TRUE
  )
  
  # Error check
  if(any(class(multi.res) == "try-error")){
    return(
      error.report(
        result = multi.res,
        SUB_FUN = "EGA.estimate",
        FUN = "EGA"
      )
    )
  }
  
  # Set up results
  if(uni.method == "expand"){
    
    if(uni.res$n.dim <= 2){ # Unidimensional 
      
      # Set results
      multi.res$wc <- uni.res$wc[!grepl("SIM", names(uni.res$wc))]
      multi.res$n.dim <- uni.res$n.dim
      
    }else if(multi.res$n.dim == 0){ # No dimensions
      
      # Set results
      multi.res$n.dim <- NA
      
    }
    
  }else{
    
    if(uni.res$n.dim == 1){ # Unidimensional 
      
      # Set results
      multi.res$wc <- uni.res$wc[!grepl("SIM", names(uni.res$wc))]
      multi.res$n.dim <- uni.res$n.dim
      
    }else if(multi.res$n.dim == 0){ # No dimensions
      
      # Set results
      multi.res$n.dim <- NA
      
    }
    
  }

  # Create dimension--variables output
  dim.variables <- data.frame(
    items = colnames(data),
    dimension = multi.res$wc
  )
  # Reorder by dimension
  multi.res$dim.variables <- dim.variables[order(dim.variables$dimension),]
  
  # Add unidimensional method
  multi.res$Methods$uni.method <- uni.method
  
  # Add type of EGA
  multi.res$EGA.type <- ifelse(multi.res$n.dim <= 2, "Unidimensional EGA", "Traditional EGA")

  # Replace cor.data with correlation
  multi.res$correlation <- multi.res$cor.data
  multi.res$cor.data <- NULL
  
  # Reorder output
  multi.res <- multi.res[
    c(
      "network", "wc", "n.dim", "correlation",
      "dim.variables", "EGA.type", "Methods"
    )
  ]

  class(multi.res) <- "EGA"

  if(isTRUE(plot.EGA)){
    
    multi.res$Plot.EGA <- suppressPackageStartupMessages(
      plot(multi.res, plot.args = plot.args)
    )
    
  }

  # Return EGA
  return(multi.res)
}
#----
