#' @title Unique Variable Analysis
#' 
#' @description Identifies locally dependent (redundant) variables in a 
#' multivariate dataset using the \code{\link[EGAnet]{EBICglasso.qgraph}} 
#' network estimation method and weighted topological overlap
#' (see Christensen, Garrido, & Golino, 2023 for more details)
#'
#' @param data Matrix or data frame.
#' Should consist only of variables to be used in the analysis.
#' Can be raw data or a correlation matrix.
#' Defaults to \code{NULL}
#' 
#' @param network Symmetric matrix or data frame.
#' A symmetric network.
#' Defaults to \code{NULL}
#' 
#' If both \code{data} and \code{network} are provided,
#' then \code{UVA} will use the \code{network}
#' with the \code{data} (rather than estimating a 
#' network from the \code{data})
#' 
#' @param n Numeric (length = 1).
#' Sample size if \code{data} provided is a correlation matrix.
#' Defaults to \code{NULL}
#' 
#' @param key Character vector (length = \code{ncol(data)}).
#' Item key for labeling variables in the results
#' 
#' @param uva.method Character (length = 1).
#' Whether the method described in Christensen, Garrido, and
#' Golino (2023) publication in \emph{Multivariate Behavioral Research}
#' (\code{"MBR"}) or Christensen, Golino, and Silvia (2020) publication
#' in \emph{European Journal of Personality} (\code{"EJP"}) should be used.
#' Defaults to \code{"MBR"}
#' 
#' Based on simulation and accumulating empirical evidence, the methods
#' described in Christensen, Golino, and Silvia (2020) such as 
#' adaptive alpha are \strong{outdated}. Evidence supports using a 
#' single cut-off value (regardless of continuous, polytomous, or
#' dichotomous data; Christensen, Garrido, & Golino, 2023)
#' 
#' @param cut.off Numeric (length = 1).
#' Cut-off used to determine when pairwise \code{\link[EGAnet]{wto}}
#' values are considered locally dependent (or redundant).
#' Must be values between \code{0} and \code{1}.
#' Defaults to \code{0.25}
#' 
#' This cut-off value is \strong{recommended} and based on extensive simulation
#' (Christensen, Garrido, & Golino, 2023). Printing the result will
#' provide a gradient of pairwise redundancies in increments of 0.20,
#' 0.25, and 0.30. Use \code{print} or \code{summary} on the output
#' rather than adjusting this cut-off value
#' 
#' @param reduce Logical (length = 1).
#' Whether redundancies should be reduced in data.
#' Defaults to \code{TRUE}
#' 
#' @param reduce.method Character (length = 1).
#' Method to reduce redundancies.
#' Available options:
#' 
#' \itemize{
#' 
#' \item{\code{"latent"} --- }
#' {Computes latent variables using \code{\link[lavaan]{cfa}} when 
#' there are three or more redundant variables. If variables are not 
#' all coded in the same direction, then they will be recoded as necessary.
#' A warning will be produced for all variables that are flipped}
#' 
#' \item{\code{"mean"} --- }
#' {Computes mean of redundant variables. If variables are not all coded in the
#' same direction, then they will be recoded as necessary.
#' A warning will be produced for all variables that are flipped}
#' 
#' \item{\code{"remove"} --- }
#' {Removes all but one variable from a set of redundant variables}
#' 
#' \item{\code{"sum"} --- }
#' {Computes sum of redundant variables. If variables are not all coded in the
#' same direction, then they will be recoded as necessary.
#' A warning will be produced for all variables that are flipped}
#' 
#' }
#' 
#' @param auto Logical (length = 1).
#' Whether \code{reduce} should occur automatically. For
#' \code{reduce.method = "remove"}, the automated decision
#' process is as follows:
#' 
#' \itemize{
#' 
#' \item{\code{Two variables} --- }
#' {The variable with the lowest maximum \code{\link[EGAnet]{wto}} 
#' to all other variables (other than the one it is redundant with)
#' is retained and the other is removed}
#' 
#' \item{\code{Three or more variables} --- }
#' {The variable with the highest mean \code{\link[EGAnet]{wto}}
#' to all other variables that are redundant with one another
#' is retained and all others are removed}
#' 
#' }
#' 
#' @param verbose Boolean (length = 1).
#' Whether messages and (insignificant) warnings should be output.
#' Defaults to \code{FALSE} (silent calls).
#' Set to \code{TRUE} to see all messages and warnings for every function call
#' 
#' @param ... Additional arguments that should be passed on to 
#' old versions of \code{UVA} or to
#' \code{\link[EGAnet]{EGA}} and
#' \code{\link[lavaan]{cfa}}
#' 
#' @examples 
#' # Perform UVA
#' uva.wmt <- UVA(wmt2[,7:24])
#' 
#' # Show summary
#' summary(uva.wmt)
#' 
#' @references 
#' \strong{Most recent simulation and implementation} \cr
#' Christensen, A. P., Garrido, L. E., & Golino, H. (2023).
#' Unique variable analysis: A network psychometrics method to detect local dependence.
#' \emph{Multivariate Behavioral Research}.
#' 
#' \strong{Conceptual foundation and outdated methods} \cr
#' Christensen, A. P., Golino, H., & Silvia, P. J. (2020).
#' A psychometric network perspective on the validity and validation of personality trait questionnaires.
#' \emph{European Journal of Personality}, \emph{34}(6), 1095-1108.
#' 
#' \strong{Weighted topological overlap} \cr
#' Nowick, K., Gernat, T., Almaas, E., & Stubbs, L. (2009).
#' Differences in human and chimpanzee gene expression patterns define an evolving network of transcription factors in brain.
#' \emph{Proceedings of the National Academy of Sciences}, \emph{106}, 22358-22363.
#' 
#' \strong{Selection of CFA Estimator} \cr
#' Rhemtulla, M., Brosseau-Liard, P. E., & Savalei, V. (2012).
#' When can categorical variables be treated as continuous? A comparison of robust continuous and categorical SEM estimation methods under suboptimal conditions.
#' \emph{Psychological Methods}, \emph{17}(3), 354-373.
#' 
#' @export
# Unique Variable Analysis ----
# Updated 08.08.2023
UVA <- function(
    data = NULL, network = NULL, n = NULL, key = NULL,
    uva.method = c("MBR", "EJP"),
    cut.off = 0.25, reduce = TRUE,
    reduce.method = c("latent", "mean", "remove", "sum"),
    auto = TRUE, verbose = FALSE, ... # `EGA` and {lavaan} arguments
)
{
  
  # Argument errors
  UVA_errors(data, network, n, cut.off, reduce, auto, verbose)
  
  # Set default method
  uva.method <- set_default(uva.method, "MBR", UVA)
  
  # Get ellipse
  ellipse <- list(...)
  
  # Check for method ("EJP" is old, not recommended)
  if(
    uva.method == "ejp" || !auto ||
    "type" %in% names(ellipse) && 
    ellipse$type %in% c("adapt", "alpha", "threshold")
  ){
    
    # Check for type = "adapt" or "alpha" 
    # (will not be supported after version 2.0.0)
    UVA_type_warning(ellipse)
    
    # Check for "manual" (will not be supported after version 2.0.0)
    if(!auto){UVA_manual_warning(auto)}
  
    # Return legacy UVA
    return(
      do.call( # Perform legacy UVA
        what = "oldUVA", # call old UVA function
        args = as.list( # force list into call
          legacy_UVA( # grab input from function calls
            data = data, n = n, key = key, reduce = reduce,
            reduce.method = reduce.method, auto = auto,
            FUN.args = list(...) # any other lingering arguments
          )
        )
      )
    )
    
  }
  
  # Set defaults
  reduce.method <- set_default(reduce.method, "remove", UVA)
  
  # Get EGA output (regardless)
  ega_output <- EGA(data, plot.EGA = FALSE, verbose = verbose, ...)
  
  # Get network
  if(!is.null(data) && is.null(network)){
    network <- ega_output$network
  }
  
  # Get key
  key <- swiftelse(is.null(key), dimnames(network)[[2]], key)
  
  # Compute weighted topological overlap
  wto_output <- abs(wto(network))
  
  # Compute descriptives
  descriptives <- wto_descriptives(wto_output)
  
  # Cut-off indices
  wto_indices <- descriptives$pairwise[
    descriptives$pairwise$wto >= cut.off,, drop = FALSE
  ]
  
  # Assign names from key (needs to be after cut-off)
  descriptives$pairwise$node_i <- key[descriptives$pairwise$node_i]
  descriptives$pairwise$node_j <- key[descriptives$pairwise$node_j]
  
  # Check for whether any redundancies exist
  if(dim(wto_indices)[1] == 0){
    
    # Return NULLs
    results <- list(
      redundant = NULL, network = network,
      wto = list(
        matrix = wto_output, # wTO matrix
        pairwise = descriptives$pairwise, # pairwise wTO
        descriptives = descriptives$basic, # basic statistics
        cut_off = cut.off # cut-off used
      )
    )
    
  }else{ # Branch for redundancies
    
    # Combine indices into a list
    redundant_variables <- get_redundancy_list(wto_output, wto_indices)
    
    # Check for reduction
    if(!reduce){
      
      # Assign names from key
      redundant_variables <- lapply(redundant_variables, function(x){key[x]})
      names(redundant_variables) <- cvapply(names(redundant_variables), function(x){key[as.numeric(x)]})
      
      # If no reduction, then return results
      results <- list(
        redundant = redundant_variables,
        network = network,
        wto = list(
          matrix = wto_output, # wTO matrix
          pairwise = descriptives$pairwise, # pairwise wTO
          descriptives = descriptives$basic, # basic statistics
          cut_off = cut.off # cut-off used
        )
      )
      
    }else{ # Proceed with reduction of data
      
      # All methods except "remove" require raw data (for now...)
      if(reduce.method != "remove" & is.null(data)){
        
        stop(
          paste0(
            "Raw data are necessary to perform reduction using reduce.method = \"",
            reduce.method, "\""
          ), call. = FALSE
        )
        
      }
      
      # Should be automatic reduction from here, so get proper function
      reduce_FUN <- switch(
        reduce.method,
        "latent" = reduce_latent,
        "mean" = reduce_mean,
        "remove" = reduce_remove,
        "sum" = reduce_sum
      )
      
      # Set up arguments
      reduce_ARGS <- list(
        data = data, wto_output = wto_output,
        redundant_variables = redundant_variables,
        correlation_matrix = ega_output$correlation,
        ellipse = ellipse
      )
      
      # Call reduction function
      reduced_results <- do.call(reduce_FUN, reduce_ARGS)
      
      # Check for "remove" method
      if(reduce.method == "remove"){
        
        # Get remove variables
        remove <- reduced_results$remove
        
        # Determine what data input were
        if(!is.null(data)){
          
          # Check for correlation matrix
          if(is_symmetric(data)){ # Correlation matrix
            data <- data[-remove, -remove]
          }else{ # Raw data
            data <- data[, -remove]
          }
          
        }else if(!is.null(network)){
          data <- network[-remove, -remove]
        }
        
      }else{ # Send reduced data
        data <- reduced_results
      }
      
      # Assign names from key
      redundant_variables <- lapply(redundant_variables, function(x){key[x]})
      names(redundant_variables) <- cvapply(names(redundant_variables), function(x){key[as.numeric(x)]})
      
      # Set up results to be returned
      results <- list(
        reduced_data = data,
        redundant = redundant_variables,
        network = network,
        wto = list(
          matrix = wto_output, # wTO matrix
          pairwise = descriptives$pairwise, # pairwise wTO
          descriptives = descriptives$basic, # basic statistics
          cut_off = cut.off # cut-off used
        )
      )
      
      # For "remove", send variables retained and removed
      if(reduce.method == "remove"){
        results$keep_remove <- list(
          keep = key[reduced_results$keep],
          remove = key[reduced_results$remove]
        )
      }
      
    }
    
  }
  
  # Add "methods" attribute
  attr(results, "methods") <- list(
    uva.method = uva.method, cut.off = cut.off, reduce = reduce
  )
  
  # Add "UVA" class
  class(results) <- "UVA"
  
  # Return results
  return(results)
  
}

# Bug checking ----
# # Select Five Factor Model personality items only
# idx <- na.omit(match(gsub("-", "", unlist(psychTools::spi.keys[1:5])), colnames(psychTools::spi)))
# items <- psychTools::spi[,idx]
# 
# # Change names in redundancy output to each item's description
# key.ind <- match(colnames(items), as.character(psychTools::spi.dictionary$item_id))
# key <- as.character(psychTools::spi.dictionary$item[key.ind])
# 
# data = items; network = NULL; n = NULL; key = key;
# cut.off = 0.25; reduce = TRUE; reduce.method = "remove";
# auto = TRUE; label.latent = FALSE; verbose = FALSE
# EGAnet.version = packageVersion("EGAnet"); uva.method = "MBR"
# ellipse = list()

#' @noRd
# Argument errors ----
# Updated 04.08.2023
UVA_errors <- function(data, network, n, cut.off, reduce, auto, verbose)
{
  
  # 'data' errors
  if(!is.null(data)){
    object_error(data, c("matrix", "data.frame"))
  }
  
  # 'network' errors
  if(!is.null(network)){
    object_error(network, c("matrix", "data.frame"))
  }
  
  # 'n' errors
  if(!is.null(n)){
    length_error(n, 1)
    typeof_error(n, "numeric")
  }
  
  # 'cut.off' errors
  length_error(cut.off, 1)
  typeof_error(cut.off, "numeric")
  range_error(cut.off, c(0, 1))
  
  # 'reduce' errors
  length_error(reduce, 1)
  typeof_error(reduce, "logical")
  
  # 'auto' errors
  length_error(auto, 1)
  typeof_error(auto, "logical")
  
  # 'verbose' errors
  length_error(verbose, 1)
  typeof_error(verbose, "logical")
  
}


#' @exportS3Method 
# S3Method Print Method ----
# Updated 04.08.2023
print.UVA <- function(x, ...)
{
  
  # Obtain wTO matrix
  wto_matrix <- x$wto$pairwise
  
  # Obtain wTO > 0.30
  large <- wto_matrix$wto > 0.30
  
  # Obtain wTO > 0.25
  moderate <- wto_matrix$wto < 0.30 & wto_matrix$wto > 0.25
  
  # Obtain wTO > 20
  small <- wto_matrix$wto < 0.25 & wto_matrix$wto > 0.20
  
  # Prepare for print
  ## Print 0.30
  cat("Variable pairs with wTO > 0.30 (large-to-very large redundancy)\n")
  if(any(large)){
    cat("\n")
    print(wto_matrix[large,], quote = FALSE, row.names = FALSE, digits = 3)
  }
  ## Print 0.25
  cat("\n----\n")
  cat("\nVariable pairs with wTO > 0.25 (moderate-to-large redundancy)\n")
  if(any(moderate)){
    cat("\n")
    print(wto_matrix[moderate,], quote = FALSE, row.names = FALSE, digits = 3)
  }
  ## Print 0.20
  cat("\n----\n")
  cat("\nVariable pairs with wTO > 0.20 (small-to-moderate redundancy)\n")
  if(any(small)){
    cat("\n")
    print(wto_matrix[small,], quote = FALSE, row.names = FALSE, digits = 3)
  }
  
}

#' @exportS3Method 
# S3Method Summary Method ----
# Updated 25.07.2023
summary.UVA <- function(object, ...)
{
  print(object, ...) # same as print
}

#' @noRd
# Obtain descriptives ----
# Updated 07.08.2023
wto_descriptives <- function(wto_output){
  
  # Get dimensions
  dimensions <- dim(wto_output)
  
  # Column sequence
  dimension_sequence <- seq_len(dimensions[2])
  
  # Initialize data frame
  wto_long <- fast.data.frame(
    c(
      rep(dimension_sequence, each = dimensions[2]),
      rep(dimension_sequence, times = dimensions[2]),
      as.vector(wto_output)
    ), nrow = length(wto_output), ncol = 3,
    colnames = c("node_i", "node_j", "wto")
  )
  
  # Subset to remove duplicates
  wto_long <- wto_long[wto_long$node_i < wto_long$node_j,]
  
  # Remove all values below zero
  wto_long <- wto_long[wto_long$wto > 0,]
  
  # Compute MAD, RANGE, QUANTILE
  MAD <- mad(wto_long$wto, constant = 1, na.rm = TRUE)
  RANGE <- range(wto_long$wto, na.rm = TRUE)
  QUANTILE <- quantile(wto_long$wto, probs = c(0.975, 0.995))
  names(QUANTILE) = c("95%", "99%")
  
  # Return list
  return(
    list(
      basic = round(
        c(
          "mean" = mean(wto_long$wto, na.rm = TRUE),
          "sd" = sd(wto_long$wto, na.rm = TRUE),
          "minimum" = RANGE[1],
          "maximum" = RANGE[2],
          "median" = median(wto_long$wto, na.rm = TRUE),
          "mad" = MAD,
          "mad3" = MAD * 3,
          "mad6" = MAD * 6,
          QUANTILE
        ), 3
      ),
      pairwise = wto_long[order(wto_long$wto, decreasing = TRUE),]
    )
  )
  
}

#' @noRd
# Get the redundancy list ----
# Updated 24.07.2023
get_redundancy_list <- function(wto_output, wto_indices)
{
  
  # Obtain the node columns of indices
  node_columns <- as.matrix(wto_indices[, c("node_i", "node_j")])
  
  # Obtain descending order of frequencies of each index
  index_frequencies <- sort(fast_table(node_columns), decreasing = TRUE)
  
  # Get index sums
  index_wto_sums <- nvapply(
    as.numeric(names(index_frequencies)), function(target_index){
      
      # Obtain vector of overlap indices
      overlapping_indices <- unlist(
        node_columns[
          node_columns[,"node_i"] == target_index |
          node_columns[,"node_j"] == target_index,
        ]
      )
      
      # Return sum
      return(
        sum(
          wto_output[
            target_index, 
            overlapping_indices[overlapping_indices != target_index]
          ], na.rm = TRUE
        )
      )
      
    }
  )
  
  # Create order based on:
  # 1. number of times node is redundant
  # 2. the total weight of redundancies
  ordered_redundancy <- names(
    index_frequencies[
      order(index_frequencies, index_wto_sums, decreasing = TRUE)
    ]
  )
  
  # Create list
  redundancy_list <- vector("list", length(ordered_redundancy))
  names(redundancy_list) <- ordered_redundancy
  
  # Extract redundancy list
  while(length(ordered_redundancy) != 0){
    
    # Check for node in node columns
    pairwise_exists <- ordered_redundancy[1] == node_columns
    
    # If any exist, then extract them
    if(any(pairwise_exists)){
      
      # Find where pairwise exists
      pairwise_row <- rowSums(pairwise_exists) != 0
      
      # Extract nodes from rows
      extracted_nodes <- as.vector(node_columns[pairwise_row,])
      
      # Add values to redundancy list
      redundancy_list[[ordered_redundancy[1]]] <-
        extracted_nodes[extracted_nodes != ordered_redundancy[1]]
      
      # Update node columns
      node_columns <- node_columns[!pairwise_row,, drop = FALSE]
      
    }
    
    # At the end, remove element from ordered redundancy
    ordered_redundancy <- ordered_redundancy[-1]
    
  }

  # Return list
  return(redundancy_list[!lvapply(redundancy_list, is.null)])
  
}

#' @noRd
# Obtain redundant variables ----
# Updated 25.07.2023
obtain_redundant_variables <- function(redundant_variables, index)
{
  
  # Check for whether all indices
  if(index == "all"){
    
    # Obtain named node
    named_node <- as.numeric(names(redundant_variables))
    
    # Obtain element node(s)
    element_node <- unname(unlist(redundant_variables))

  }else{
    
    # Obtain named node
    named_node <- as.numeric(names(redundant_variables)[index])
    
    # Obtain element node(s)
    element_node <- redundant_variables[[index]]
    
  }
  
  # Return all nodes (redundant but ensure unique)
  return(unique(c(named_node, element_node)))
  
}

#' @noRd
# Create composite variable matrix ----
# Updated 25.07.2023
create_composite_variable_matrix <- function(data, redundant_variables)
{
  
  # Get redundant variable length
  redundant_length <- length(redundant_variables)
  
  # Return matrix for composite variables
  return(
    matrix(
      NA, nrow = dim(data)[1], ncol = redundant_length,
      dimnames = list(
        NULL, paste0(
          "CV", 
          format_integer(
            seq_len(redundant_length),
            digits(redundant_length) - 1
          )
        )
      )
    )
  )
  
}

#' @noRd
# Recode ----
# Not ideal but must happen to ensure proper values
# Updated 25.07.2023
recode <- function(data, all_names, correlation_matrix, ellipse)
{
  
  # Get data dimensions
  dimensions <- dim(data)
  
  # Get variables
  variables <- data[, all_names, drop = FALSE]
  
  # Get variable categories
  variable_categories <- data_categories(variables)
  
  # Get correlations
  variable_correlations <- correlation_matrix[all_names, all_names]
  
  # Set diagonal to zero
  diag(variable_correlations) <- 0
  
  # Get signs (`obtain_signs` is in `net.loads` internals)
  variable_signs <- attr(obtain_signs(variable_correlations), "signs")
  
  # Get negative signs
  negative_signs <- variable_signs == -1
  
  # Check for any negative signs
  if(any(negative_signs)){
    
    # Check for dominant direction
    dominant_sums <- sum(variable_signs)
    
    # Get ordinal categories
    ordinal.categories <- swiftelse(
      "ordinal.categories" %in% names(ellipse),
      ellipse$ordinal.categories, 7
    )
    
    # If greater or equal to zero, then flip negative signs
    if(dominant_sums >= 0){
      
      # Get subtraction vector
      subtraction_vector <- swiftelse(
        variable_categories[negative_signs] <= ordinal.categories,
        ordinal.categories, 0
      )
      
      # Get variable names that will be flipped
      flipped_variables <- names(negative_signs)[negative_signs]
      
      # Loop over variables and flip them
      variables[, negative_signs] <- nvapply(
        seq_along(subtraction_vector), function(index){
        subtraction_vector[index] - 
        variables[, flipped_variables[index]]
      }, LENGTH = dimensions[1])
    
      
    }else{ # If less than zero, then flip positive signs
      
      # Get positive signs
      positive_signs <- !negative_signs
      
      # Get subtraction vector
      subtraction_vector <- swiftelse(
        variable_categories[positive_signs] <= ordinal.categories,
        ordinal.categories, 0
      )
      
      # Get variable names that will be flipped
      flipped_variables <- names(positive_signs)[positive_signs]

      # Loop over variables and flip them
      variables[, positive_signs] <- nvapply(
        seq_along(subtraction_vector), function(index){
          subtraction_vector[index] - 
            variables[, flipped_variables[index]]
        }, LENGTH = dimensions[1])
      
    }
    
    # Send warning about variables flipped
    warning(
      paste(
        "Some variables were recoded to ensure proper aggregation:",
        paste0(flipped_variables, collapse = ", ")
      ), call. = FALSE
    )
    
  }
  
  # Return variables
  return(variables)
  
}

#' @noRd
# "Latent" reduce method ----
# Updated 25.07.2023
reduce_latent <- function(
    data, wto_output, # not used
    redundant_variables, 
    correlation_matrix,
    ellipse
)
{
  
  # Get case sequence
  case_sequence <- nrow_sequence(data)
  
  # Get variable names
  variable_names <- dimnames(data)[[2]]
  
  # Get redundant length
  redundant_length <- length(redundant_variables)
  
  # Create matrix for composite variables
  composite_variables <- create_composite_variable_matrix(
    data, redundant_variables
  )
  
  # Get {lavaan}'s CFA function
  cfa_FUN <- silent_load(lavaan::cfa)
  
  # Get {lavaan} CFA arguments
  lavaan_ARGS_copy <- obtain_arguments(cfa_FUN, ellipse)
  
  # Send message about latent variable estimation
  message("Estimating latent variables...", appendLF = FALSE)
  
  # Loop over redundant variables
  for(index in seq_len(redundant_length)){
    
    # Refresh {lavaan} arguments
    lavaan_ARGS <- lavaan_ARGS_copy
    
    # Get variable names (from nodes)
    all_names <- variable_names[
      obtain_redundant_variables(
        redundant_variables, index = index
      )
    ]
    
    # Make CFA model (in `helpers.R`)
    lavaan_ARGS$model <- make_unidimensional_cfa(all_names)
    
    # Send data
    lavaan_ARGS$data <- recode(
      data, all_names, correlation_matrix, ellipse
    )
  
    # Get estimator arguments (in `helpers.R`)
    lavaan_ARGS <- estimator_arguments(lavaan_ARGS, ellipse)
    
    # Set `std.lv` to `TRUE`
    lavaan_ARGS$std.lv <- TRUE
    
    # Estimate latent variable model
    cfa_estimate <- do.call("cfa_FUN", lavaan_ARGS)
    
    # Identify available cases
    available_cases <- lavaan::inspect(cfa_estimate, what = "case.idx")
    
    # Obtain latent score
    latent_score <- lavaan::lavPredict(
      cfa_estimate, optim.method = "nlminb"
      # "nlminb" is faster than "bfgs" (with no difference in scores)
    )
    
    # Compute new composite
    composite_variables[
      case_sequence %in% available_cases, index
    ] <- latent_score
    
    # Update message
    message(
      paste0("\rEstimating latent variables... (", index, " of ", redundant_length, " complete)"),
      appendLF = FALSE
    )
    
  }
  
  # Update message
  message(
    paste0(
      "\rEstimating latent variables...done.",
      paste0(rep(" ", 13 + digits(redundant_length)), collapse = "")
    )
  )
  
  # Obtain all redundant nodes
  remove_variables <- obtain_redundant_variables(
    redundant_variables, index = "all"
  )

  # Return data with new composites
  return(cbind(data[,-remove_variables], composite_variables))
  
}

#' @noRd
# "Mean" reduce method ----
# Updated 25.07.2023
reduce_mean <- function(
    data, wto_output, # not used
    redundant_variables, 
    correlation_matrix,
    ellipse # not used
)
{
  
  # Create matrix for composite variables
  composite_variables <- create_composite_variable_matrix(
    data, redundant_variables
  )
  
  # Compute means into the composite variable matrix
  composite_variables[] <- nvapply(seq_along(redundant_variables), function(index){
    
    # Obtain all nodes
    all_nodes <- obtain_redundant_variables(
      redundant_variables, index
    )
    
    # Return new composite
    return(
      rowMeans(
        recode(data, all_nodes, correlation_matrix, ellipse),
        na.rm = TRUE
      )
    )
 
  }, LENGTH = dim(data)[1])

  # Obtain all redundant nodes
  remove_variables <- obtain_redundant_variables(
    redundant_variables, index = "all"
  )
  
  # Return data with new composites
  return(cbind(data[,-remove_variables], composite_variables))
  
}

#' @noRd
# "Remove" reduce method ----
# Updated 25.07.2023
reduce_remove <- function(
    data, # not used
    wto_output, redundant_variables, 
    correlation_matrix, ellipse
)
{
  
  # Loop over redundant variables and return keep and remove
  selection_list <- lapply(seq_along(redundant_variables), function(index){
    
    # Obtain all nodes
    all_nodes <- obtain_redundant_variables(
      redundant_variables, index
    )
    
    # Determine whether to use wTO or standard deviation
    if(length(all_nodes) > 2){
      
      # Selection index based on maximum average 
      # wTO value to other redundant variables
      selection_index <- which.max(
        colMeans(wto_output[all_nodes, all_nodes], na.rm = TRUE)
      )
      
    }else{ # Only two nodes
      
      # Selection index based on lowest maximum 
      # wTO value to all other variables
      selection_index <- which.min(
        apply(wto_output[all_nodes, -all_nodes], 1, max, na.rm = TRUE)
      )
      
    }
    
    # Return list
    return(
      list(
        keep = all_nodes[selection_index], 
        remove = all_nodes[-selection_index]
      )
    )
  
  })
  
  # Return results
  return(
    list(
      keep = sort(ulapply(selection_list, function(x){x$keep})),
      remove = sort(ulapply(selection_list, function(x){x$remove}))
    )
  )
  
}

#' @noRd
# "Sum" reduce method ----
# Updated 25.07.2023
reduce_sum <- function(
    data, wto_output, # not used
    redundant_variables, 
    correlation_matrix,
    ellipse # not used
)
{
  
  # Create matrix for composite variables
  composite_variables <- create_composite_variable_matrix(
    data, redundant_variables
  )
  
  # Compute means into the composite variable matrix
  composite_variables[] <- nvapply(seq_along(redundant_variables), function(index){
    
    # Obtain all nodes
    all_nodes <- obtain_redundant_variables(
      redundant_variables, index
    )
    
    # Return new composite
    return(
      rowSums(
        recode(data, all_nodes, correlation_matrix, ellipse),
        na.rm = TRUE
      )
    )
    
  }, LENGTH = dim(data)[1])
  
  # Obtain all redundant nodes
  remove_variables <- obtain_redundant_variables(
    redundant_variables, index = "all"
  )
  
  # Return data with new composites
  return(cbind(data[,-remove_variables], composite_variables))
  
}

#' @noRd
# Legacy arguments for "oldUVA.R" ----
# Updated 08.08.2023
legacy_UVA <- function(
    data, n, key, reduce,
    reduce.method, auto,
    FUN.args
)
{
  
  # Old UVA arguments
  oldUVA.args <- obtain_arguments(
    FUN = oldUVA, FUN.args = FUN.args # Need to set up defaults
  )
  
  # Replace necessary arguments
  oldUVA.args$data <- data; oldUVA.args$n <- n; oldUVA.args$key <- key
  oldUVA.args$reduce <- reduce; oldUVA.args$reduce.method <- reduce.method
  oldUVA.args$auto <- auto;
  
  # Return arguments
  return(oldUVA.args)
  
}
  
#' @noRd
# "type" warning ---
# Updated 24.07.2023
UVA_type_warning <- function(ellipse)
{
  # Check for "type"
  if("type" %in% names(ellipse) && ellipse$type != "threshold"){
  
    warning(
      paste0(
        "Argument `type = \"", ellipse$type, "\"` will not be supported ",
        "in future versions of {EGAnet}. Recent evidence suggests that ",
        "`cut.off = 0.25` is best practice:",
        "\n\nChristensen, A. P., Garrido, L. E., & Golino, H. (2023). ",
        "Unique variable analysis: A network psychometrics method to ",
        "detect local dependence. ",
        styletext("Multivariate Behavioral Research", defaults = "italics"),
        ", 1-18. https://doi.org/10.1080/00273171.2023.2194606",
        "\n\nDo not submit error reports. Bugs will not be fixed"
      ), call. = FALSE
    )
    
  }
}

#' @noRd
# "auto" is `FALSE` warning ---
# Updated 07.08.2023
UVA_manual_warning <- function(auto)
{
  
  # Check for manual
  if(!auto){
    warning(
      paste0(
        "Manual decisions (`auto = FALSE`) will not be supported ",
        "in future versions of {EGAnet}.",
        "\n\nUse `reduce = FALSE` to perform manual inspection instead",
        "\n\nDo not submit error reports. Bugs will not be fixed"
      ), call. = FALSE
    )
  }
  
}
