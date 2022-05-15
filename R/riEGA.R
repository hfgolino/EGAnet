#' Random-Intercept \code{\link[EGAnet]{EGA}}
#'
#' Estimates the number of substantive dimensions after controlling for wording effects. 
#' EGA is applied to a residual correlation matrix after subtracting and random intercept
#' factor with equal unstandardized loadings from all the regular and unrecoded
#' reversed items in the database
#'
#' @param data Matrix or data frame.
#' Variables (down columns) or correlation matrix.
#' If the input is a correlation matrix,
#' then argument \code{n} (number of cases) is \strong{required}.
#' Variables \strong{MUST} be unrecoded -- reversed items should
#' \strong{remain} reversed
#'
#' @param n Integer.
#' Sample size if \code{data} provided is a correlation matrix
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
#' or \code{TMFG}
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
#' @param estimator Character.
#' Estimator to use for random-intercept model (see \href{https://lavaan.ugent.be/tutorial/est.html}{Estimators}
#' for more details).
#' Defaults to \code{"auto"}, which selects \code{"MLR"} for continuous data and
#' \code{"WLSMV"} for mixed and categorical data.
#' Data are considered continuous data if they have 6 or
#' more categories (see Rhemtulla, Brosseau-Liard, & Savalei, 2012)
#' 
#' @param lavaan.args List.
#' If \code{reduce.method = "latent"}, then \code{\link{lavaan}}'s \code{\link[lavaan]{cfa}}
#' function will be used to create latent variables to reduce variables.
#' Arguments should be input as a list. Some example arguments 
#' (see \code{\link[lavaan]{lavOptions} for full details})
#' 
#' @return Returns a list containing:
#' 
#' \item{EGA}{Results from \code{\link[EGAnet]{EGA}}}
#' 
#' \item{RI}{A list containing information about the random-intercept
#' model (if the model converged):
#' 
#' \itemize{
#' 
#' \item{fit}
#' {The fit object for the random-intercept model using \code{\link[lavaan]{cfa}}}
#' 
#' \item{lavaan.args}
#' {The arguments used in \code{\link[lavaan]{cfa}}}
#' 
#' \item{loadings}
#' {Standardized loadings from the random-intercept model}
#' 
#' \item{correlation}
#' {Residual correlations after accounting for the random-intercept model}
#' 
#' }
#' 
#' 
#' 
#' }
#'
#' @author  
#' Alejandro Garcia-Pardina <alejandrogp97@gmail.com>,
#' Francisco J. Abad <fjose.abad@uam.es>,
#' Alexander P. Christensen <alexpaulchristensen@gmail.com>,
#' Hudson Golino <hfg9s at virginia.edu>,
#' Luis Eduardo Garrido <luisgarrido@pucmm.edu.do>, and
#' Robert Moulder <rgm4fd@virginia.edu>
#'
#' @examples
#' # Obtain example data
#' data <- optimism
#' 
#' \dontrun{
#' # riEGA example
#' opt.res <- riEGA(data = optimism)
#' }
#'
#' @references 
#' # Selection of CFA Estimator \cr
#' Rhemtulla, M., Brosseau-Liard, P. E., & Savalei, V. (2012).
#' When can categorical variables be treated as continuous? A comparison of robust continuous and categorical SEM estimation methods under suboptimal conditions.
#' \emph{Psychological Methods}, \emph{17}, 354-373.
#'
#' @export
#' 
# Random-Intercept EGA
# Changed from 'residualEGA.R' on 17.04.2022
# Updated 13.05.2022
riEGA <- function(
    data, n = NULL, uni.method = c("expand", "LE"),
    corr = c("cor_auto", "pearson", "spearman"),
    model = c("glasso", "TMFG"), model.args = list(),
    algorithm = c("walktrap", "louvain"), algorithm.args = list(),
    consensus.iter = 100,
    consensus.method = c(
      "highest_modularity",
      "most_common",
      "iterative",
      "lowest_tefi"
    ),
    plot.EGA = TRUE, plot.type = c("GGally", "qgraph"),
    plot.args = list(), estimator = c("auto", "WLSMV", "MLR"),
    lavaan.args = list()
  )
{
  
  #### ARGUMENTS HANDLING ####
  
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
  
  if(missing(consensus.method)){
    consensus.method <- "highest_modularity"
  }else{consensus.method <- match.arg(consensus.method)}
  
  if(missing(plot.type)){
    plot.type <- "GGally"
  }else{plot.type <- match.arg(plot.type)}
  
  #### ARGUMENTS HANDLING ####
  
  # Ensure data is a matrix
  data <- as.matrix(data)
  
  # Check for number of variables
  if(is.null(n)){
    
    # Check for sample
    if(nrow(data) != ncol(data)){
      
      # Obtain sample size
      n <- nrow(data)
      
    }else{
      stop("Argument 'n' must be specified for square matrices.")
    }
    
  }
  
  # Obtain variable names
  # (ensure these names are OK with lavaan)
  data <- lavaan.formula.names(data)
  variables <- colnames(data)
  
  # Configure random-intercept model
  ri_model <- paste(
    "RI =~",
    paste(
      "1*", variables, sep = "", collapse = " + "
    )
  )
  
  # Check for user changes to default arguments
  if(length(lavaan.args) == 0){
    # Set standard arguments
    lavaan.args <- formals(lavaan::cfa)
    lavaan.args[length(lavaan.args)] <- NULL
  }else{
    # Obtain default arguments
    lavaan.default <- formals(lavaan::cfa)
    lavaan.default[length(lavaan.default)] <- NULL
    
    # Switch out any defaults with user-specificed arguments
    if(any(names(lavaan.args) %in% names(lavaan.default))){
      lavaan.default[names(lavaan.args)] <- lavaan.args
    }
    
    # Change arguments
    lavaan.args <- lavaan.default
  }
  
  # Set random-intercept model
  lavaan.args$model <- ri_model
  
  # Set data
  lavaan.args$data <- data
  
  # Obtain estimator
  if(missing(estimator)){
    estimator <- "auto"
  }else{
    estimator <- match.arg(estimator)
  }
  
  # Check for "auto" estimator
  if(estimator == "auto"){
    
    # Obtain categories
    categories <- apply(data, 2, function(x){
      length(na.omit(unique(x)))
    })
    
    # Check categories
    if(sum(categories < 6) > 1){# Not all continuous
      lavaan.args$estimator <- "WLSMV"
      lavaan.args$missing <- "pairwise"
      lavaan.args$ordered <- TRUE
    }else{# All can be considered continuous
      lavaan.args$estimator <- "MLR"
      lavaan.args$missing <- "fiml"
    }
    
  }
  
  # Get CFA function from lavaan
  FUN <- lavaan::cfa
  
  # Fit model
  fit <- try(
    do.call(what = "FUN", args = as.list(lavaan.args)),
    silent = TRUE
  )
  
  # Check for "try-error"
  if(any(class(fit) == "try-error")){
  
    # Let user know that random-intercept model did not converge
    message(
      "Random-intercept model did not converge suggesting that there is not substantial evidence for wording effects."
    )
    
  }else{ # Proceed with random-intercept model
    
    # Obtain random-intercept loadings
    ri_loadings <- lavaan::inspect(fit, what = "std")$lambda
    
    # Obtain residual correlation matrix
    ri_correlations <- lavaan::residuals(fit, type = "cor")$cov
    
    # Change diagonal to 1
    diag(ri_correlations) <- 1
    
    # Ensure positive definite
    if(any(eigen(ri_correlations)$values < 0)){
      
      ri_correlations <- as.matrix(
        Matrix::nearPD(
          ri_correlations,
          corr = TRUE,
          keepDiag = TRUE
        )$mat
      )
      
    }
    
  }
  
  # Obtain EGA defaults
  ega_defaults <- formals(EGA)
  
  # Remove "..."
  ega_defaults <- ega_defaults[-length(ega_defaults)]
  
  # Insert values
  if(exists("ri_correlations")){
    
    # Set correlations
    ega_defaults$data <- ri_correlations
    
  }else{
    
    # Set data
    ega_defaults$data <- data
    
  }
  
  # Set the rest of the arguments
  ega_defaults$n <- n
  ega_defaults$uni.method <- uni.method
  ega_defaults$corr <- corr
  ega_defaults$model <- model
  ega_defaults$model.args <- model.args
  ega_defaults$algorithm <- algorithm
  ega_defaults$algorithm.args <- algorithm.args
  ega_defaults$consensus.method <- consensus.method
  ega_defaults$consensus.iter <- consensus.iter
  ega_defaults$plot.EGA <- plot.EGA
  ega_defaults$plot.type <- plot.type
  ega_defaults$plot.args <- plot.args
  
  # Get EGA
  suppressPackageStartupMessages(
    ega_result <- do.call(
      EGA, ega_defaults
    )
  )

  # Return results
  results <- list()
  results$EGA <- ega_result
  
  # Random-intercept model
  results$RI <- list()
  results$RI$fit <- fit
  results$RI$lavaan.args <- lavaan.args
  
  # Check for whether random-intercept model converged
  if(exists("ri_correlations")){
    
    # Report loadings and correlations
    results$RI$loadings <- ri_loadings
    results$RI$correlation <- ri_correlations
    
    # Message user about recoding
    recoding_message <- paste(
      "\nRandom-intercept model converged. Wording effects likely. Results are only valid if data are ",
      styletext(
        "unrecoded", defaults = "underline"
      ), ".\n",
      sep = ""
    )
    
    message(recoding_message)
  
  }
  
  # Obtain methods
  methods <- list()
  
  # Lavaan methods
  methods$estimator <- lavaan.args$estimator
  methods$corr <- corr
  methods$model <- model
  methods$model.args <- model.args
  methods$algorithm <- algorithm
  methods$algorithm.args <- algorithm.args
  methods$uni.method <- uni.method
  
  # Add methods to results
  results$Methods <- methods
  
  # Make class "riEGA"
  class(results) <- "riEGA"
  
  return(results)
}
