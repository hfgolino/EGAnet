#' Detects Redundant Variables in a Multivariate Dataset
#' 
#' @description Identifies redundant variables in a multivariate dataset
#' using a number of different association methods and types of significance values
#' (see Christensen, Garrido, & Golino, 2020 for more details)
#'
#' @param data Matrix or data frame.
#' Input can either be data or a correlation matrix
#' 
#' @param n Numeric.
#' If input in \code{data} is a correlation matrix, 
#' then sample size is required.
#' Defaults to \code{NULL}
#' 
#' @param method Character.
#' Computes weighted topological overlap (\code{"wTO"} using \code{\link[qgraph]{EBICglasso}}),
#' partial correlations (\code{"pcor"}), and correlations (\code{"cor"}).
#' Defaults to \code{"wTO"}
#' 
#' @param type Character. Type of significance.
#' Computes significance using the standard \emph{p}-value (\code{"alpha"}),
#' adaptive alpha \emph{p}-value (\code{\link[NetworkToolbox]{adapt.a}}), 
#' or some threshold \code{"threshold"}.
#' Defaults to \code{"adapt"} 
#'
#' @param sig Numeric.
#' \emph{p}-value for significance of overlap (defaults to \code{.05}).
#' Defaults for \code{"threshold"} for each \code{method}:
#' 
#' \itemize{
#' 
#' \item{\code{"wTO"}}
#' {.20}
#' 
#' \item{\code{"pcor"}}
#' {.20}
#' 
#' \item{\code{"cor"}}
#' {.70}
#' 
#' } 
#' 
#' @param key Character vector.
#' A vector with variable descriptions that correspond
#' to the order of variables input into \code{data}
#' 
#' 
#' 
#' @param reduce Character.
#' How should data be reduced?
#' Defaults to \code{"latent"}
#' 
#' \itemize{
#' 
#' \item{\code{"latent"}}
#' {Redundant variables will be combined into a latent variable}
#' 
#' \item{\code{"remove"}}
#' {All but one redundant variable will be removed}
#' 
#' }
#' 
#' @param lavaan.args List.
#' If \code{reduce = "latent"}, then \code{\link{lavaan}}'s \code{\link[lavaan]{cfa}}
#' function will be used to create latent variables to reduce variables.
#' Arguments should be input as a list. Some example arguments 
#' (see \code{\link[lavaan]{lavOptions} for full details}:
#' 
#' \itemize{
#' 
#' \item{\code{estimator}}
#' {Estimator to use for latent variables (see \href{https://lavaan.ugent.be/tutorial/est.html}{Estimators})
#' for more details}
#' 
#' \item{\code{missing}}
#' {How missing data should be handled}
#' 
#' \item{\code{std.lv}}
#' {If \code{TRUE}, the metric of each latent variable is determined by fixing their (residual) variances to 1.0.
#' If \code{FALSE}, the metric of each latent variable is determined by fixing the factor loading of the first
#' indicator to 1.0. If there are multiple groups, \code{std.lv = TRUE} and \code{"loadings"} is included in the
#' \code{group.label} argument, then only the latent variances i of the first group will be fixed to 1.0, while
#' the latent variances of other groups are set free.}
#' 
#' }
#' 
#' @param adhoc Boolean.
#' Should adhoc check of redundancies be performed?
#' Defaults to \code{TRUE}.
#' If \code{TRUE}, adhoc check will run the redundancy analysis
#' on the reduced variable set to determine if there are any remaining
#' redundancies. This check is performed with the arguments:
#' \code{method = "wTO"}, \code{type = "threshold"}, and \code{sig = .20}.
#' This check is based on Christensen, Garrido, and Golino's (2020)
#' simulation where these parameters were found to be the most conservative,
#' demonstrating few false positives and false negatives
#'
#' @param plot.redundancy Boolean.
#' Should redundancies be plotted in a network plot?
#' Defaults to \code{FALSE}
#' 
#' @param plot.args List.
#' Arguments to be passed onto \code{\link[GGally]{ggnet2}}.
#' Defaults:
#' 
#' \itemize{
#' 
#' \item{\code{vsize = 6}}{}
#' 
#' \item{\code{alpha = 0.4}}{}
#' 
#' \item{\code{label.size = 5}}{}
#' 
#' \item{\code{edge.alpha = 0.7}}{}
#' 
#' }
#'
#' @return Returns a list:
#' 
#' \item{redundancy}{A list containing several objects:
#' 
#' \itemize{
#' 
#' \item{\code{redudant}}
#' {Vectors nested within the list corresponding
#' to redundant nodes with the name of object in the list}
#' 
#' \item{\code{data}}
#' {Original data}
#' 
#' \item{\code{correlation}}
#' {Correlation matrix of original data}
#' 
#' \item{\code{weights}}
#' {Weights determine by weighted topological overlap,
#' partial correlation, or zero-order correlation}
#' 
#' \item{\code{network}}
#' {If \code{method = "wTO"}, then
#' the network computed following \code{\link[EGAnet]{EGA}} with
#' \code{\link[qgraph]{EBICglasso}} network estimation}
#' 
#' \item{\code{plot}}
#' {If \code{redundancy.plot = TRUE}, then
#' a plot of all redundancies found}
#' 
#' \item{\code{descriptives}}{
#' 
#' \itemize{
#' 
#' \item{basic}
#' {A vector containing the mean, standard deviation,
#' median, median absolute deviation (MAD), 3 times the MAD, 6 times the MAD,
#' minimum, maximum, and critical value for the overlap measure
#' (i.e., weighted topological overlap, partial correlation, or threshold)}
#' 
#' \item{centralTendency}
#' {A matrix for all (absolute) non-zero values and their
#' respective standard deviation from the mean and median absolute deviation
#' from the median}
#' 
#' }
#' }
#' 
#' \item{\code{method}}
#' {Returns \code{method} argument}
#' 
#' \item{\code{type}}
#' {Returns \code{type} argument}
#' 
#' \item{\code{distribution}}
#' {If \code{type != "threshold"}, then 
#' distribution that was used to determine significance}
#' 
#' }
#' 
#' }
#' 
#' \item{reduced}{A list containing:
#' 
#' \itemize{
#' 
#' \item{\code{data}}
#' {New data with redundant variables merged or removed}
#'
#' \item{\code{merged}}{A matrix containing the variables that were
#' decided to be redundant with one another}
#' 
#' }
#' 
#' }
#' 
#' \item{adhoc}{If \code{adhoc = TRUE}, then
#' the adhoc check containing the same objects as in
#' the \code{redundancy} list object in the output
#' }
#'
#' @examples
#' # Select Five Factor Model personality items only
#' idx <- na.omit(match(gsub("-", "", unlist(psychTools::spi.keys[1:5])), colnames(psychTools::spi)))
#' items <- psychTools::spi[,idx]
#' 
#' # Change names in redundancy output to each item's description
#' key.ind <- match(colnames(items), as.character(psychTools::spi.dictionary$item_id))
#' key <- as.character(psychTools::spi.dictionary$item[key.ind])
#' 
#' if(interactive()){
#' redundancy.analysis(data = items, method = "wTO", type = "adapt",
#'                     key = key, reduce = "latent")
#' }
#'
#' @references
#' # Simulation using \code{redundancy.analysis} \cr
#' Christensen, A. P., Garrido, L. E., & Golino, H. (2020).
#' A novel approach for detecting redundant variables in multivariate data.
#' \emph{PsyArXiv}.
#' 
#' # Implementation of \code{redundancy.analysis} (formally \code{node.redundant}) \cr
#' Christensen, A. P., Golino, H., & Silvia, P. J. (in press).
#' A psychometric network perspective on the validity and validation of personality trait questionnaires.
#' \emph{European Journal of Personality}.
#' doi: \href{https://doi.org/10.1002/per.2265}{10.1002/per.2265}
#' 
#' # wTO measure \cr
#' Nowick, K., Gernat, T., Almaas, E., & Stubbs, L. (2009).
#' Differences in human and chimpanzee gene expression patterns define an evolving network of transcription factors in brain.
#' \emph{Proceedings of the National Academy of Sciences}, \emph{106}, 22358-22363.
#' doi: \href{https://doi.org/10.1073/pnas.0911376106}{10.1073/pnas.0911376106}
#'
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#'
#' @importFrom stats pgamma pnorm qgamma qnorm cov2cor mad
#'
#' @export
#
# Redundant Nodes Function
# Updated 12.12.2020
redundancy.analysis <- function(data, n = NULL,
                                method = c("cor", "pcor", "wTO"),
                                type = c("adapt", "alpha", "threshold"), sig,
                                key = NULL, reduce = c("latent", "remove"),
                                lavaan.args = list(), adhoc = TRUE,
                                plot.redundancy = FALSE, plot.args = list()
                                )
{
  # Missing and NULL arguments
  ## n
  if(nrow(data) == ncol(data)){
    if(is.null(n)){
      stop("Argument 'n' must be set for square matrices")
    }
  }else{### Get n
    n <- nrow(data)
    ### Compute correlation matrix
    cormat <- qgraph::cor_auto(data)
    ### Make sure it's positive definite
    if(any(eigen(cormat)$values < 0)){
      cormat <- as.matrix(Matrix::nearPD(cormat, keepDiag = TRUE)$mat)
    }
  }
  
  ## method
  if(missing(method)){
    method <- tolower("wTO")
  }else{method <- tolower(match.arg(method))}
  
  ## type
  if(missing(type)){
    type <- "adapt"
  }else{type <- tolower(match.arg(type))}
  
  ## sig
  if(missing(sig)){
    if(type == "threshold"){
      sig <- switch(method,
                    "cor" = .70,
                    "pcor" = .20,
                    "wto" = .20
      )
    }else{sig <- .05}
  }
  
  ## key
  if(is.null(key)){
    key <- colnames(data)
  }
  
  ## reduce
  if(missing(reduce)){
    reduce <- "latent"
  }else{reduce <- match.arg(reduce)}
  
  ## plot.args
  if(length(plot.args) == 0){
    plot.args <-list(vsize = 6, alpha = 0.4, label.size = 5, edge.alpha = 0.7)
  }else{
    plot.args <- plot.args
    plots.arg1 <- list(vsize = 6, label.size = 5, alpha = 0.4, edge.alpha = 0.7)
    plot.args.use <- plot.args
    
    if(any(names(plots.arg1) %in% names(plot.args.use))){
      
      plot.replace.args <- plots.arg1[na.omit(match(names(plot.args.use), names(plots.arg1)))]
      
      plot.args <- c(plot.args.use,plots.arg1[names(plots.arg1) %in% names(plot.args.use)==FALSE])}
  }
  
  # Perform redundancy analysis
  process <- redundancy.process(data = data, cormat = cormat,
                                n = n, method = method,
                                type = type, sig = sig,
                                plot.redundancy = plot.redundancy, plot.args = plot.args)
  
  
  # Get names
  if(!is.null(key)){
    process <- redund.names(node.redundant.obj = process, key = key)
  }
  
  # Run through redundancy reduction
  reduced <- redund.reduce(node.redundant.obj = process,
                           reduce = reduce,
                           plot.args = plot.args,
                           lavaan.args = lavaan.args)
  
  # Check for any remaining redundancies
  if(adhoc){
    ## Message user
    message("Running adhoc check for any potential redundancies remaining")
    
    ## Run check
    adhoc.check <- suppressMessages(
      redundancy.process(data = reduced$data, cormat = qgraph::cor_auto(reduced$data),
                         n = nrow(reduced$data), method = "wto",
                         type = "threshold", sig = .20,
                         plot.redundancy = FALSE, plot.args = plot.args)
    )
    
    if(all(is.na(adhoc.check$redundant))){
      
      message("Some redundancies may still exist. See `OUTPUT$adhoc`")
      
    }else{message("No redundancies reamin.")}
  }
  
  
  # Full results
  res <- list()
  res$redundancy <- process
  res$reduced <- reduced
  if(adhoc){res$adhoc <- adhoc.check}
  
  return(res)
  
}
  