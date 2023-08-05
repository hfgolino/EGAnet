#' Unique Variable Analysis
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
#' @param model Character.
#' A string indicating the method to use.
#' Current options are:
#'
#' \itemize{
#'
#' \item{\strong{\code{glasso}}}
#' {Estimates the Gaussian graphical model using graphical LASSO with
#' extended Bayesian information criterion to select optimal regularization parameter.
#' This is the default method}
#'
#' \item{\strong{\code{TMFG}}}
#' {Estimates a Triangulated Maximally Filtered Graph}
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
#' @param method Character.
#' Computes weighted topological overlap (\code{"wTO"} using \code{\link[qgraph]{EBICglasso}}),
#' partial correlations (\code{"pcor"}), or correlations (\code{"cor"})
#' Defaults to \code{"wTO"}
#' 
#' @param type Character. Type of significance.
#' Computes significance using the standard \emph{p}-value (\code{"alpha"}),
#' adaptive alpha \emph{p}-value (\code{adapt.a}), 
#' or some threshold \code{"threshold"}.
#' Defaults to \code{"threshold"} 
#'
#' @param sig Numeric.
#' \emph{p}-value for significance of overlap (defaults to \code{.05}).
#' Defaults for \code{"threshold"} for each \code{method}:
#' 
#' \itemize{
#' 
#' \item{\code{"wTO"}}
#' {.25}
#' 
#' \item{\code{"pcor"}}
#' {.35}
#' 
#' \item{\code{"cor"}}
#' {.50}
#' 
#' } 
#' 
#' @param key Character vector.
#' A vector with variable descriptions that correspond
#' to the order of variables input into \code{data}.
#' Defaults to \code{NULL} or the column names of \code{data}
#' 
#' @param reduce Boolean.
#' Should redundancy reduction be performed?
#' Defaults to \code{TRUE}.
#' Set to \code{FALSE} for redundancy analysis only
#' 
#' @param auto Boolean.
#' Should redundancy reduction be automated?
#' Defaults to \code{TRUE}.
#' Set to \code{FALSE} for manual selection
#' 
#' @param label.latent Boolean.
#' Should latent variables be labelled?
#' Defaults to \code{TRUE}.
#' Set to \code{FALSE} for arbitrary labelling (i.e., "LV_")
#' 
#' @param reduce.method Character.
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
#' \item{\code{"sum"}}
#' {Redundant variables are combined by summing across cases (rows)}
#' 
#' }
#' 
#' @param lavaan.args List.
#' If \code{reduce.method = "latent"}, then \code{\link{lavaan}}'s \code{\link[lavaan]{cfa}}
#' function will be used to create latent variables to reduce variables.
#' Arguments should be input as a list. Some example arguments 
#' (see \code{\link[lavaan]{lavOptions} for full details}):
#' 
#' \itemize{
#' 
#' \item{\code{estimator}}
#' {Estimator to use for latent variables (see \href{https://lavaan.ugent.be/tutorial/est.html}{Estimators})
#' for more details. Defaults to \code{"MLR"} for continuous data and \code{"WLSMV"} for mixed and categorical data.
#' Data are considered continuous data if they have 6 or more categories (see Rhemtulla, Brosseau-Liard, & Savalei, 2012)}
#' 
#' \item{\code{missing}}
#' {How missing data should be handled. Defaults to \code{"fiml"}}
#' 
#' \item{\code{std.lv}}
#' {If \code{TRUE}, the metric of each latent variable is determined by fixing their (residual) variances to 1.0.
#' If \code{FALSE}, the metric of each latent variable is determined by fixing the factor loading of the first
#' indicator to 1.0. If there are multiple groups, \code{std.lv = TRUE} and \code{"loadings"} is included in the
#' \code{group.label} argument, then only the latent variances i of the first group will be fixed to 1.0, while
#' the latent variances of other groups are set free.
#' Defaults to \code{TRUE}}
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
#' \item{\code{vsize = 6}}{Changes node size}
#' 
#' \item{\code{alpha = 0.4}}{Changes transparency}
#' 
#' \item{\code{label.size = 5}}{Changes label size}
#' 
#' \item{\code{edge.alpha = 0.7}}{Changes edge transparency}
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
#' \item{reduced}{If \code{reduce = TRUE}, then a list containing:
#' 
#' \itemize{
#' 
#' \item{\code{data}}
#' {New data with redundant variables merged or removed}
#'
#' \item{\code{merged}}{A matrix containing the variables that were
#' decided to be redundant with one another}
#' 
#' \item{\code{method}}{Method used to perform redundancy reduction}
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
#' \dontrun{
#' # Automated selection of local dependence (default)
#' uva.results <- UVA(data = items, key = key)
#' 
#' # Produce Methods section
#' methods.section(uva.results)}
#' 
#' # Manual selection of local dependence
#' if(interactive()){
#' uva.results <- UVA(data = items, key = key, auto = FALSE)}
#'
#' @references
#' # Simulation using \code{UVA} \cr
#' Christensen, A. P., Garrido, L. E., & Golino, H. (under review).
#' Unique Variable Analysis: A novel approach for detecting redundant variables in multivariate data.
#' \emph{PsyArXiv}.
#' 
#' # Implementation of \code{UVA} (formally \code{node.redundant}) \cr
#' Christensen, A. P., Golino, H., & Silvia, P. J. (2020).
#' A psychometric network perspective on the validity and validation of personality trait questionnaires.
#' \emph{European Journal of Personality}, \emph{34}, 1095-1108.
#' 
#' # wTO measure \cr
#' Nowick, K., Gernat, T., Almaas, E., & Stubbs, L. (2009).
#' Differences in human and chimpanzee gene expression patterns define an evolving network of transcription factors in brain.
#' \emph{Proceedings of the National Academy of Sciences}, \emph{106}, 22358-22363.
#' 
#' # Selection of CFA Estimator \cr
#' Rhemtulla, M., Brosseau-Liard, P. E., & Savalei, V. (2012).
#' When can categorical variables be treated as continuous? A comparison of robust continuous and categorical SEM estimation methods under suboptimal conditions.
#' \emph{Psychological Methods}, \emph{17}, 354-373.
#'
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#'
#' @noRd
#
# Unique Variable Analysis
# Superceded on 02.02.2023
# Updated 13.09.2022
oldUVA <- function(
  data, n = NULL,
  model = c("glasso", "TMFG"),
  corr = c("cor_auto", "pearson", "spearman"),
  method = c("cor", "pcor", "wTO"),
  type = c("adapt", "alpha", "threshold"), sig,
  key = NULL, reduce = TRUE, auto = TRUE, label.latent = FALSE,
  reduce.method = c("latent", "remove", "sum"),
  lavaan.args = list(), adhoc = TRUE,
  plot.redundancy = FALSE, plot.args = list()
)
{
  
  # Send not refactored message
  not_refactored("oldUVA")
  
  # Make sure data is a matrix
  data <- as.matrix(data)
  
  # Missing and NULL arguments
  ## corr
  if(missing(corr)){
    corr <- "cor_auto"
  }else{corr <- match.arg(corr)}
  
  ## model
  if(missing(model)){
    model <- tolower("glasso")
  }else{model <- tolower(match.arg(model))}
  
  ## method
  if(missing(method)){
    method <- tolower("wTO")
  }else{method <- tolower(match.arg(method))}
  
  ## type
  if(missing(type)){
    type <- "threshold"
  }else{type <- tolower(match.arg(type))}
  
  ## sig
  if(missing(sig)){
    if(type == "threshold"){
      sig <- switch(method,
                    "cor" = .50,
                    "pcor" = .35,
                    "wto" = .25
      )
    }else{sig <- .05}
  }
  
  ## reduce
  if(missing(reduce.method)){
    reduce.method <- "latent"
  }else{reduce.method <- match.arg(reduce.method)}
  
  ## n
  if(nrow(data) == ncol(data)){
    
    cormat <- data
    
    if(is.null(n)){
      stop("Argument 'n' must be set for square matrices")
    }
    
    ### Let user know that variables can't be combined
    if(isTRUE(reduce)){
      
      if(reduce.method != "remove"){
        reduce.method <- "remove"
        message("Input is a square matrix. Changing 'reduce.method' to \"remove\"")
      }
      
      
    }
    
  }else{### Get n
    n <- nrow(data)
    
    ## Compute correlation matrix
    cormat <- switch(corr,
                     "cor_auto" = auto.correlate(data),
                     "pearson" = cor(data, use = "pairwise.complete.obs"),
                     "spearman" = cor(data, method = "spearman", use = "pairwise.complete.obs")
    )
    ### Make sure it's positive definite
    if(any(eigen(cormat)$values < 0)){
      cormat <- as.matrix(Matrix::nearPD(cormat, keepDiag = TRUE)$mat)
    }
    
  }
  
  ## prepare arguments for lavaan
  if(isTRUE(reduce)){
    
    if(reduce.method == "latent"){
      
      ## lavaan.args
      if(length(lavaan.args) == 0){
        lavaan.args <- formals(lavaan::cfa)
        lavaan.args[length(lavaan.args)] <- NULL
        lavaan.args$std.lv <- TRUE
      }else{
        lavaan.default <- formals(lavaan::cfa)
        lavaan.default[length(lavaan.default)] <- NULL
        lavaan.default$std.lv <- TRUE
        
        if(any(names(lavaan.args) %in% names(lavaan.default))){
          lavaan.default[names(lavaan.args)] <- lavaan.args
        }
        
        lavaan.args <- lavaan.default
      }
      
      ## change key if NULL
      if(is.null(key)){
        data <- lavaan.formula.names(data)
      }
      
    }
    
  }
  
  ## check for automated procedure
  ### check for type first
  if(isTRUE(auto) & type != "threshold"){
    
    message("\nWarning: Automated UVA is not recommended for types other than 'reduce.method = \"threshold\"'\n'auto' set to 'FALSE'")
    auto <- FALSE
    
  }
  
  ### check for reduce method second
  if(isTRUE(auto) & isTRUE(reduce) & reduce.method == "remove"){
    
    message("\nWarning: Automated UVA is not available for 'reduce.method = \"remove\"'\n'auto' set to 'FALSE'")
    auto <- FALSE
    
  }

  ## plot.args
  plot.args <- suppressPackageStartupMessages(
    GGally_args(plot.args)
  )
  
  # Perform redundancy analysis
  process <- suppressWarnings(
    suppressMessages(
      redundancy.process(data = data, cormat = cormat,
                         n = n, model = model, method = method,
                         type = type, sig = sig,
                         plot.redundancy = plot.redundancy, plot.args = plot.args)
    )
  )
  
  ## key
  if(is.null(key)){
    
    if(is.null(colnames(data))){
      colnames(data) <- paste("V", 1:ncol(data), sep = "")
    }
    
    key <- colnames(data)
    
  }
  
  # Get names
  if(any(!is.na(process$redundant))){
    process <- redund.names(node.redundant.obj = process, key = key)
  }
  
  # Run through redundancy reduction
  if(all(is.na(process$redundant))){
    
    reduced <- process
    
    message("No redundant variables were identified.")
    
    adhoc.check <- NULL
    
  }else{
    
    ## Manual
    if(isTRUE(reduce) & !isTRUE(auto)){
      
      reduced <- redund.reduce(node.redundant.obj = process,
                               reduce.method = reduce.method,
                               plot.args = plot.args,
                               lavaan.args = lavaan.args,
                               corr = corr)
      
      # Check for any remaining redundancies
      if(adhoc){
        ## Message user
        message("Running adhoc check for any potential redundancies remaining...\n")
        
        ## Run check
        ## Compute correlation matrix
        if(isSymmetric(reduced$data)){
          cor.data <- reduced$data
        }else{
          
          cor.data <- switch(corr,
                             "cor_auto" = auto.correlate(reduced$data),
                             "pearson" = cor(reduced$data, use = "pairwise.complete.obs"),
                             "spearman" = cor(reduced$data, method = "spearman", use = "pairwise.complete.obs")
          )
          
        }
        
        adhoc.check <- suppressMessages(
          redundancy.process(data = reduced$data, cormat = cor.data,
                             n = n,
                             model = model,
                             method = method,
                             type = "threshold", sig = .20,
                             plot.redundancy = FALSE, plot.args = plot.args)
        )
        
        # Artificial pause for feel
        Sys.sleep(1)
        
      }
      
      # Artificial pause for feel
      Sys.sleep(1)
      
    }else if(isTRUE(reduce) & isTRUE(auto)){## Automated
      
      # Message user
      message("\nCombining variables...", appendLF = FALSE)
      
      # Initial reduction
      reduced <- redund.reduce.auto(
        node.redundant.obj = process,
        reduce.method = reduce.method,
        lavaan.args = lavaan.args,
        corr = corr
      )
      
      ## Run check
      ## Compute correlation matrix
      if(isSymmetric(reduced$data)){
        cor.data <- reduced$data
      }else{
        
        sink <- capture.output(
          cor.data <- suppressMessages(
            suppressWarnings(
              switch(corr,
                     "cor_auto" = auto.correlate(reduced$data),
                     "pearson" = cor(reduced$data, use = "pairwise.complete.obs"),
                     "spearman" = cor(reduced$data, method = "spearman", use = "pairwise.complete.obs")
              )
            )
          )
        )
        
      }
      
      ### Make sure it's positive definite
      if(any(eigen(cor.data)$values < 0)){
        cor.data <- as.matrix(Matrix::nearPD(cor.data, keepDiag = TRUE)$mat)
      }
      
      adhoc.check <- suppressMessages(
        redundancy.process(data = reduced$data, cormat = cor.data,
                           n = n,
                           model = model,
                           method = method,
                           type = "threshold", sig = sig,
                           plot.redundancy = FALSE, plot.args = plot.args)
      )
      
      # Check for names in key
      rename_check <- adhoc.check$redundant
      target_names <- names(rename_check) %in% names(key)
      if(any(target_names)){
        names(rename_check)[target_names] <- key[names(rename_check)[target_names]]
      }
      
      # Insert into adhoc.check
      adhoc.check$redundant <- lapply(rename_check, function(x){
        
        target_names <- x %in% names(key) 
        
        if(any(target_names)){
          x[target_names] <- key[x[target_names]]
        }
        
        return(x)
        
      })
      
      while(all(!is.na(adhoc.check$redundant))){
        
        # Adhoc reductions
        reduced <- redund.adhoc.auto(node.redundant.obj = adhoc.check,
                                     node.redundant.reduced = reduced,
                                     node.redundant.original = process,
                                     reduce.method = reduce.method,
                                     lavaan.args = lavaan.args,
                                     corr = corr)
        
        # Check for undimensionality
        if(ncol(reduced$data) == 1){
          break
        }else{
          
          ## Run check
          ## Compute correlation matrix
          if(isSymmetric(reduced$data)){
            cor.data <- reduced$data
          }else{
            
            sink <- capture.output(
              cor.data <- suppressMessages(
                suppressWarnings(
                  switch(corr,
                         "cor_auto" = auto.correlate(reduced$data),
                         "pearson" = cor(reduced$data, use = "pairwise.complete.obs"),
                         "spearman" = cor(reduced$data, method = "spearman", use = "pairwise.complete.obs")
                  )
                )
              )
            )
            
          }
          
          ### Make sure it's positive definite
          if(any(eigen(cor.data)$values < 0)){
            cor.data <- as.matrix(Matrix::nearPD(cor.data, keepDiag = TRUE)$mat)
          }
          
          adhoc.check <- suppressMessages(
            redundancy.process(data = reduced$data, cormat = cor.data,
                               n = n,
                               model = model,
                               method = "wto",
                               type = "threshold", sig = sig,
                               plot.redundancy = FALSE, plot.args = plot.args)
          )
          
          # Check for names in key
          rename_check <- adhoc.check$redundant
          target_names <- names(rename_check) %in% names(key)
          if(any(target_names)){
            names(rename_check)[target_names] <- key[names(rename_check)[target_names]]
          }
          
          # Insert into adhoc.check
          adhoc.check$redundant <- lapply(rename_check, function(x){
            
            target_names <- x %in% names(key) 
            
            if(any(target_names)){
              x[target_names] <- key[x[target_names]]
            }
            
            return(x)
            
          })
          
          
        }
        
      }
      
      # Message user
      message("done")
      
      # Message user
      if(ncol(reduced$data) == 1){
        message(
          "\nAfter combining local dependencies, data are determined to be unidimensional."
        )
      }
    
      # # Name latent variables
      # name_question <- readline(prompt = "Name latent variables? [Y/n]: ")
      # 
      # # Check for appropriate response
      # name_question <- tolower(name_question)
      # 
      # while(name_question != "y" & name_question != "n"){
      #   
      #   # Name latent variables
      #   name_question <- readline(prompt = "Inappropriate response. Try again. [Y/n]: ")
      #   
      #   # Check for appropriate response
      #   name_question <- tolower(name_question)
      #   
      # }
      
      # Name latent variables
      if(isTRUE(label.latent)){
        name_question <- "y"
      }else{
        name_question <- "n"
      }
      
      if(name_question == "y"){
        
        # Copy of reduced merged
        lats <- reduced$merged
        
        # Make it a list
        lats <- lapply(apply(lats, 1, as.list), unlist)
        lats <- lapply(lats, unname) # unname
        lats <- lapply(lats, function(x){ # make them homogeneous with spaces
          x <- na.omit(ifelse(x == "", NA, x))
          x <- c(x, "", "")
          return(x)
        })
        
        # Line break
        linebreak()
        
        # Loop through latent variables
        for(i in 1:length(lats)){
          
          # Variables in latent variable
          cat(
            paste(
              lats[[i]],
              collapse = "\n"
            )
          )
          
          # Ask for label
          lab <- readline(prompt = "New label for latent variable (no quotations): ")
          
          # Assign new label
          colnames(reduced$data)[
            which(colnames(reduced$data) == row.names(reduced$merged)[i])
          ] <- lab
          row.names(reduced$merged)[i] <- lab
          
          # Message user for progress
          message(
            paste(
              "\n", i, "of", nrow(reduced$merged), "latent variables named."
            )
          )
          
          # Line break
          linebreak()
          
        }
        
      }
      
    }else{reduced <- process}
    
  }
  
  # Full results
  res <- list()
  res$redundancy <- process
  if(reduce){
    res$reduced <- reduced
    if(adhoc){res$adhoc <- adhoc.check}
  }
  
  # Set up methods
  res$Methods <- list()
  res$Methods$method <- method
  res$Methods$type <- type
  res$Methods$sig <- sig
  if(reduce){
    
    res$Methods$reduce.method <- reduce.method
    
    if(reduce.method == "latent"){
      res$Methods$lavaan.args <- lavaan.args
    }
    
  }
  res$Methods$auto <- auto
  
  # Change reduced names for sum scores
  if(reduce.method == "sum"){
    
    colnames(res$reduced$data) <- gsub("LV_", "SUM_", colnames(res$reduced$data))
    row.names(res$reduced$merged) <- gsub("LV_", "SUM_", row.names(res$reduced$merged))
    
  }
  
    
  # Set class
  class(res) <- "UVA"
  
  return(res)
  
}
  
      