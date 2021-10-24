#' Item Stability Statistics from \code{\link[EGAnet]{bootEGA}}
#'
#' @description Based on the \code{\link[EGAnet]{bootEGA}} results, this function
#' computes and plots the number of times an item (variable) is estimated
#' in the same factor/dimension as originally estimated by \code{\link[EGAnet]{EGA}} (\code{item.replication}).
#' The output also contains each item's replication frequency (i.e., proportion of
#' bootstraps that an item appeared in each dimension; \code{item.dim.rep}) as well as the average
#' network loading for each item in each dimension (\code{item.loadings}).
#'
#' @param bootega.obj A \code{\link[EGAnet]{bootEGA}} object
#'
#' @param IS.plot Should the plot be produced for \code{item.replication}?
#' If \code{TRUE}, then a plot for the \code{item.replication} output will be produced.
#' Defaults to \code{TRUE}
#' 
#' @param ... Additional arguments.
#' Used for deprecated arguments from previous versions of \code{\link[EGAnet]{itemStability}}
#'
#' @return Returns a list containing:
#' 
#' \item{membership}{A list containing:
#'
#' \itemize{
#'
#' \item{\strong{\code{empirical}}}
#' {The empirical memberships from the empirical \code{\link[EGAnet]{EGA}} result}
#'
#' \item{\strong{\code{unique}}}
#' {The unique dimensions from the empirical \code{\link[EGAnet]{EGA}} result}
#'
#' \item{\strong{\code{bootstrap}}}
#' {The memberships from the replicate samples in the \code{\link[EGAnet]{bootEGA}} results}
#'    }
#' }
#' 
#' \item{item.stability}{A list containing:
#'
#' \itemize{
#'
#' \item{\strong{\code{empirical.dimensions}}}
#' {The proportion of times each item replicated
#' within the empirical \code{\link[EGAnet]{EGA}} defined dimension.
#' This EGA result is defined using the input from
#' \code{\link[EGAnet]{bootEGA}}}
#'
#' \item{\strong{\code{all.dimensions}}}
#' {The proportion of times each item replicated
#' in each of the empirical \code{\link[EGAnet]{EGA}} defined dimensions.
#' This EGA result is defined using the input from
#' \code{\link[EGAnet]{bootEGA}}}
#'
#'
#'    }
#' }
#'
#' \item{plot}{A plot of the number of times each item
#' replicated within the empirical \code{\link[EGAnet]{EGA}} defined dimension.}
#'
#' \item{mean.loadings}{Matrix of the average standardized network loading
#' (computed using \code{\link[EGAnet]{net.loads}}) for each item in each dimension}
#'
#' @examples
#'
#' # Load data
#' wmt <- wmt2[,7:24]
#'
#' \dontrun{# Estimate EGA network
#' ## plot.type = "qqraph" used for CRAN checks
#' ## plot.type = "GGally" is the default
#' ega.wmt <- EGA(data = wmt, model = "glasso", plot.type = "qgraph")
#'
#' # Estimate dimension stability
#' boot.wmt <- bootEGA(data = wmt, iter = 100, typicalStructure = TRUE,
#' plot.typicalStructure = TRUE, model = "glasso", plot.type = "qgraph",
#' type = "parametric", ncores = 2)
#' }
#' 
#' # Estimate item stability statistics
#' res <- itemStability(boot.wmt)
#' 
#' # Changing plot features (ggplot2)
#' ## Changing colors (ignore warnings)
#' ### qgraph Defaults
#' res$plot + 
#'     ggplot2::scale_color_manual(values = rainbow(max(res$membership$unique)))
#' 
#' ### Pastel
#' res$plot + 
#'     ggplot2::scale_color_brewer(palette = "Pastel1")
#'     
#' ## Changing Legend (ignore warnings)
#' res$plot + 
#'     ggplot2::scale_color_discrete(labels = "Intelligence")
#'
#' @references
#' Christensen, A. P., & Golino, H. (2021).
#' Estimating the stability of the number of factors via Bootstrap Exploratory Graph Analysis: A tutorial.
#' \emph{Psych}, \emph{3}(3), 479-500.
#'
#' Christensen, A. P., Golino, H., & Silvia, P. J. (2020).
#' A psychometric network perspective on the validity and validation of personality trait questionnaires.
#' \emph{European Journal of Personality}, \emph{34}(6), 1095-1108.
#'
#' @seealso \code{\link[EGAnet]{EGA}} to estimate the number of dimensions of an instrument using EGA and
#' \code{\link[EGAnet]{CFA}} to verify the fit of the structure suggested by EGA using confirmatory factor analysis.
#'
#' @author Hudson Golino <hfg9s at virginia.edu> and Alexander P. Christensen <alexpaulchristensen@gmail.com>
#'
#' @export
#Item Stability function
# Updated 27.02.2021
# Major revamp 27.02.2021 
itemStability <- function (bootega.obj, IS.plot = TRUE, ...){
  
  # Check for 'bootEGA' object
  if(class(bootega.obj) != "bootEGA")
  {stop("Input for 'bootega.obj' is not a 'bootEGA' object")}
  
  # Get additional arguments
  add.args <- list(...)
  
  # Check if 'orig.wc' has been input as an argument
  if("orig.wc" %in% names(add.args)){
    
    # Give deprecation warning
    warning(
        "The 'orig.wc' argument has been deprecated in itemStability.\n\nInstead, the empirical EGA estimated in bootEGA's results is used"
    )
  }
  
  # Check if 'item.freq' has been input as an argument
  if("item.freq" %in% names(add.args)){
    
    # Give deprecation warning
    warning(
      "The 'item.freq' argument has been deprecated in itemStability"
    )
  }
  
  # Check if 'plot.item.rep' has been input as an argument
  if("plot.item.rep" %in% names(add.args)){
    
    # Give deprecation warning
    warning(
      
      paste(
        "The 'plot.item.rep' argument has been deprecated in itemStability.\n\nInstead use: IS.plot =", add.args$plot.item.rep, sep = " "
      )
      
    )
    
    # Handle the plot appropriately
    IS.plot <- add.args$plot.item.rep
  }
  
  # Message function
  message(styletext(styletext("\nItem Stability Analysis", defaults = "underline"), defaults = "bold"))
  
  # Let user know results are being organized
  message("\nOrganizing data...", appendLF = FALSE)
  
  # Original EGA result
  ## Network
  empirical.EGA.network <- bootega.obj$EGA$network
  ## Community membership
  empirical.EGA.membership <- bootega.obj$EGA$wc
  
  # Get numeric memberships
  membership.numeric <- numeric.membership(empirical.EGA.membership)
  
  # Get unique memberships
  unique.membership <- unique(membership.numeric)
  
  # Obtain bootstrap membership matrix
  bootstrap.membership <- simplify2array(bootega.obj$boot.wc)
  
  # Homogenize memberships
  final.membership <- try(
    homogenize.membership(membership.numeric, bootstrap.membership),
    silent = TRUE
  )
  
  # Error check
  if(any(class(final.membership) == "try-error")){
    return(
      error.report(final.membership,
                   "homogenize.membership",
                   "itemStability")
    )
  }
  
  # Let user know results are done being organized
  message("done\n")
  
  # Let user know results are being computed
  message("Computing results...", appendLF = FALSE)
  
  # Get proportion table
  replication.proportion <- try(
    proportion.table(final.membership),
    silent = TRUE
  )
  
  # Error check
  if(any(class(replication.proportion) == "try-error")){
    return(
      error.report(replication.proportion,
                   "proportion.table",
                   "itemStability")
    )
  }
  
  # Adjust for missing and NA dimensions
  final.proportion <- try(
    missing.dimension.check(replication.proportion,
                            membership.numeric,
                            bootstrap.membership),
    silent = TRUE
  )
  
  # Error check
  if(any(class(final.proportion) == "try-error")){
    return(
      error.report(final.proportion,
                   "missing.dimension.check",
                   "itemStability")
    )
  }
  
  # Initialize results
  results <- list()
  
  # Add empirical results
  results$membership <- list()
  results$membership$empirical <- empirical.EGA.membership
  results$membership$unique <- unique.membership
  results$membership$bootstrap <- final.membership
  
  # Add item stability to the results
  ## Item stability for empirical dimension
  empirical.stability <- membership.numeric
  
  for(i in 1:nrow(final.proportion)){
    empirical.stability[i] <- final.proportion[i,paste(empirical.EGA.membership[i])]
  }
  
  results$item.stability <- list()
  results$item.stability$empirical.dimensions <- empirical.stability
  
  # Add full proportion matrix to the results
  final.proportion <- as.matrix(final.proportion)
  if(any(is.na(colnames(final.proportion))) | any(colnames(final.proportion) == "NA")){
    
    final.proportion <- final.proportion[,-which(is.na(colnames(final.proportion)) | colnames(final.proportion) == "NA")]
    
  }
  
  results$item.stability$all.dimensions <- final.proportion
  # as.matrix() resolves unidimensional structures
  
  # Plot
  results <- try(
    itemStability.plot(results, bootega.obj),
    silent = TRUE
  )
  
  # Error check
  if(any(class(results) == "try-error")){
    return(
      error.report(results,
                   "itemStability.plot",
                   "itemStability")
    )
  }
  
  # Reorder empirical and all dimensions stability
  results$item.stability$empirical.dimensions <- results$item.stability$empirical.dimensions[as.character(results$plot$data$Item)]
  results$item.stability$all.dimensions <- results$item.stability$all.dimensions[as.character(results$plot$data$Item),]
  
  # Compute network loadings
  loadings <- try(
    itemStability.loadings(results, bootega.obj),
    silent = TRUE
  )
  
  # Error check
  if(any(class(loadings) == "try-error")){
    return(
      error.report(loadings,
                   "itemStability.loadings",
                   "itemStability")
    )
  }
  
  # Insert loadings into results
  results$mean.loadings <- loadings
  
  # Reorder mean loadings
  results$mean.loadings <- results$mean.loadings[as.character(results$plot$data$Item),]
  
  # Let user know results are done being organized
  message("done\n")
  
  # Set class
  class(results) <- "itemStability"
  
  # Plot to user?
  if(IS.plot){
    plot(results$plot)
  }
  
  return(results)
  
}
#----
