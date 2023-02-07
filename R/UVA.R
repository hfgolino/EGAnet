#' Unique Variable Analysis
#' 
#' @description Identifies redundant variables in a multivariate dataset
#' using the \code{\link[qgraph]{EBICglasso}} network estimation method
#' and weighted topological overlap
#' (see Christensen, Garrido, & Golino, 2020 for more details)
#'
#' @param data Matrix, data frame, or symmetric matrix.
#' Input can either be data or a correlation matrix
#' 
#' @param network Symmetric matrix.
#' Input for a symmetric network matrix.
#' 
#' If both \code{data} and \code{network} are provided,
#' then \code{UVA} will proceed to use the input network
#' with the data (rather than estimating a network)
#' 
#' @param n Numeric vector (length = 1).
#' If input in \code{data} is a correlation matrix, 
#' then sample size is required.
#' Defaults to \code{NULL}
#' 
#' @param cut.off Numeric vector (length = 1).
#' Must be values between \code{0} and \code{1}.
#' Defaults to \code{0.25}
#' 
#' @param key Character vector (length = \code{ncol(data)}).
#' Item key for labeling items
#' 
#' @param reduce Logical (length = 1).
#' Whether redundancies should be reduced in data.
#' Defaults to \code{TRUE}
#' 
#' @param reduce.method Character.
#' Method to reduce redundancies:
#' 
#' \itemize{
#' 
#' \item{\code{"latent"}}
#' {Computes latent variables when there are three or more
#' redundant variables. Computes sum otherwise.
#' }
#' 
#' \item{\code{"remove"}}
#' {Removes all but one variable from a set of redundant variables
#' }
#' 
#' \item{\code{"sum"}}
#' {Computes sum of redundant variables
#' }
#' 
#' }
#' 
#' @param auto Logical (length = 1).
#' Whether \code{reduce} should occur automatically using different rules
#' depending on reduction method:
#' 
#' \itemize{
#' 
#' \item{\code{"latent"}}
#' {Computes latent variables when there are three or more
#' redundant variables. Computes sum otherwise.
#' }
#' 
#' \item{\code{"remove"}}
#' {Removes all but one variable from a set of redundant variables.
#' Keeps the variable with the highest variable-total correlation using
#' the correlation between each variable and the sum of all other
#' variables with the target variable removed. For ties and two variables,
#' the variable with the largest standard deviation is kept.
#' }
#' 
#' \item{\code{"sum"}}
#' {Computes sum of redundant variables
#' }
#' 
#' }
#' 
#' @param label.latent Boolean (length = 1).
#' When \code{reduce.method = "latent"}, should
#' latent variables be labelled?
#' Defaults to \code{FALSE}.
#' Set to \code{TRUE} to type your own labels 
#' 
#' @param EGAnet.version Character (length = 1).
#' \code{\link{EGAnet}} version used to perform previous
#' \code{UVA}.
#' Defaults to current version.
#' Set version to previous version to reproduce results
#' from an older version
#' 
#' @param ... Additional arguments.
#' Arguments that should be passed onto old versions of \code{UVA}
#' (use the \code{EGAnet.version} argument to set the version),
#' \code{\link[EGAnet]{EGA}}, and \code{\link[lavaan]{cfa}}
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
#' # Results with no reduction
#' no_reduce_uva <- UVA(
#'   data = items, key = key, reduce = FALSE
#' )
#' 
#' # Show summary
#' summary(no_reduce_uva)
#' 
#' \dontrun{
#' # Results with automatic reduction
#' reduced_uva <- UVA(
#'   data = items, key = key, reduce = TRUE,
#'   reduce.method = "latent" 
#' )}
#' 
#' @references 
#' # Simulation using UVA
#' Christensen, A. P., Garrido, L. E., & Golino, H. (under review).
#' Unique variable analysis: A network psychometrics method to detect local dependence.
#' \emph{PsyArXiv}.
#' 
#' # Implementation of UVA (formally node.redundant)
#' Christensen, A. P., Golino, H., & Silvia, P. J. (2020).
#' A psychometric network perspective on the validity and validation of personality trait questionnaires.
#' \emph{European Journal of Personality}, \emph{34}(6), 1095-1108.
#' 
#' # wTO measure
#' Nowick, K., Gernat, T., Almaas, E., & Stubbs, L. (2009).
#' Differences in human and chimpanzee gene expression patterns define an evolving network of transcription factors in brain.
#' \emph{Proceedings of the National Academy of Sciences}, \emph{106}, 22358-22363.
#' 
#' # Selection of CFA Estimator
#' Rhemtulla, M., Brosseau-Liard, P. E., & Savalei, V. (2012).
#' When can categorical variables be treated as continuous? A comparison of robust continuous and categorical SEM estimation methods under suboptimal conditions.
#' \emph{Psychological Methods}, \emph{17}(3), 354-373.
#' 
#' @export
# Unique Variable Analysis ----
# Updated 07.02.2023
UVA <- function(
    data = NULL, network = NULL, n = NULL, key = NULL, cut.off = 0.25,
    reduce = TRUE, reduce.method = c("latent", "mean", "remove", "sum"),
    auto = TRUE, label.latent = FALSE,
    EGAnet.version = packageVersion("EGAnet"),
    ... # `EGA`, {lavaan}, and "EGAnet.version" arguments
)
{
  
  # If version <= 1.2.4, push to old UVA
  ## Keeps legacy code alive
  if(version_conversion(EGAnet.version) <= 124){
    
    return(
      do.call( # Perform legacy UVA
        what = "oldUVA", # call old UVA function
        args = as.list( # force list into call
          legacy_UVA( # grab input from function calls
            data, n, key, cut.off, reduce,
            reduce.method, auto, label.latent,
            FUN.args = list(...) # any other lingering arguments
          )
        )
      )
    )
    
  }
  
  # Set default for "reduce.method"
  if(missing(reduce.method)){
    reduce.method <- "latent"
  }else{reduce.method <- match.arg(reduce.method)}
  
  # Check for valid inputs
  ## Objects (see "helpers-errors.R" for functions)
  if(!is.null(data)){
    object_error(input = data, expected_type = c("matrix", "data.frame"))
  }else if(!is.null(network)){
    object_error(input = network, expected_type = c("matrix", "data.frame"))
  }else{
    stop("Both 'data' and 'network' are `NULL`. Expected at least one input.")
  }
  
  ## Values
  if(!is.null(n)){value_error(input = n, expected_value = "numeric")};
  if(!is.null(key)){value_error(key, "character")}; value_error(cut.off, "numeric"); 
  value_error(reduce, "logical"); value_error(reduce.method, "character"); 
  value_error(auto, "logical"); value_error(label.latent, "logical")
  
  ## Lengths
  if(!is.null(n)){length_error(input = n, expected_lengths = 1)};
  if(!is.null(key)){length_error(key, ncol(data))}; length_error(cut.off, 1);
  length_error(reduce, 1); length_error(reduce.method, 1);
  length_error(auto, 1); length_error(label.latent, 1);
  
  ## Ranges
  if(!is.null(n)){range_error(input = n, expected_ranges = c(2, Inf))};
  range_error(cut.off, c(0, 1));
  
  # Get data to network
  if(!is.null(data)){
    
    ## Check for "n"
    if(is.null(n) & !is_symmetric(data)){
      n <- nrow(data)
    }
    
    ## Symmetric matrix
    symmetric_matrix_error(data = data, n = n)
    
    ## Ensure dimension names (see "helpers-general.R" for functions)
    data <- ensure_dimension_names(data)
    
    ## Check for network
    if(is.null(network)){
      
      ## EGA arguments (see "helpers-functions.R" for functions)
      ega_ARGS <- ega_arguments(arguments = list(...))
      
      ## Update "data", "n", and "plot.EGA" arguments
      ega_ARGS$data <- data; ega_ARGS$n <- n; ega_ARGS$plot.EGA <- FALSE;
      
      ## Perform EGA
      ega_output <- do.call(
        what = "EGA",
        args = ega_ARGS
      )
      
      ## Extract output for weighted topological overlap
      network <- ega_output$network
      
    }
    
  }
  
  # At this point, there should be a network somewhere
  ## Symmetric matrix
  symmetric_matrix_error(data = network, n = n)
  
  ## Ensure dimension names (see "helpers-general.R" for functions)
  network <- ensure_dimension_names(network)
  
  # Compute weighted topological overlap
  wto_output <- abs(wto(network, signed = FALSE))
  
  # Compute descriptives
  descriptives <- wto_descriptives(
    wto_output = wto_output, key = key
  )
  
  # Cut-off indices
  wto_indices <- wto_cut_off(
    wto_output = wto_output,
    cut_off = cut.off
  )
  
  # Check for whether any redundancies exist
  if(nrow(wto_indices) == 0){
    
    # Return NULLs
    results <- list(
      redundant = NULL,
      network = network,
      wto = list(
        matrix = wto_output, # wTO matrix
        cut_off = cut.off, # cut-off used
        descriptives = descriptives # descriptives of wTO matrix
      )
    )
    
    # Ensure "UVA" class
    class(results) <- "UVA"
    
    # Return results
    return(results)
    
  }
  
  # Combine indices into a list
  overlapping_list <- wto_index_list(
    wto_output = wto_output, wto_indices = wto_indices
  )
  
  # Condense list to remove duplicate entries
  # And remove elements that are included elsewhere
  redundant_variables <- condense_overlap(
    overlapping_list = overlapping_list
  )
  
  # Determine whether data should be reduced
  if(isTRUE(reduce)){
    
    # Force lowercase for character input
    reduce.method <- tolower(reduce.method)
    
    # All methods except "remove" require raw data (for now...)
    if(reduce.method != "remove" & is.null(data)){
      
      stop(
        paste0(
          "Raw data are necessary to perform reduction using reduce.method = \"",
          reduce.method, "\""
        )
      )
      
    }
    
    # Check for whether {lavaan} arguments are necessary
    if(reduce.method == "latent"){
      
      ## {lavaan} arguments (see "helpers-functions.R" for functions)
      lavaan_ARGS <- lavaan_arguments(arguments = list(...))
    
    }
    
    # Determine whether to use automated procedure
    if(isTRUE(auto)){
      
      # Use switch to grab function
      reduce_FUN <- switch(
        reduce.method,
        "latent" = reduce_latent,
        "mean" = reduce_mean,
        "remove" = reduce_remove,
        "sum" = reduce_sum
      )
      
      # Set up arguments
      reduce_ARGS <- list(
        data = data,
        wto_output = wto_output,
        redundant_variables = redundant_variables
      )
      
      # Check for whether {lavaan} arguments are necessary
      if(reduce.method == "latent"){
        reduce_ARGS$lavaan_ARGS <- lavaan_ARGS
      }
      
      # Call reduction function
      reduced_results <- do.call(
        what = reduce_FUN, args = reduce_ARGS 
      )
      
      # Organize output for each method
      if(reduce.method != "remove"){
        
        # Obtain reduced 
        reduced_data <- reduced_results
        
      }else{
        
        # Obtain `keep` and `remove` variables
        keep <- reduced_results$keep; remove <- reduced_results$remove;
        
        # Check for whether raw data was used
        if(!is.null(data) & !is_symmetric(data)){
          
          # Remove variables from data
          reduced_data <- data[,-remove]
          
        }else if(is_symmetric(data)){
          
          # Remove variables from correlation matrix
          reduced_data <- data[-remove, -remove]
          
        }else if(!is.null(network)){
          
          # Remove variables from network
          reduced_network <- network[-remove, -remove]
          
          # Original data can be returned
          results <- list(
            reduced_network = reduced_network,
            original_network = network,
            redundant = redundant_variables,
            wto = list(
              matrix = wto_output, # wTO matrix
              cut_off = cut.off, # cut-off used
              descriptives = descriptives # descriptives of wTO matrix
            )
          )
          
        }
        
      }
      
      # Generalized reduction results return
      if(!is.null(data)){
        
        # Original data can be returned
        results <- list(
          reduced_data = reduced_data,
          original_data = data,
          redundant = redundant_variables,
          network = network,
          wto = list(
            matrix = wto_output, # wTO matrix
            cut_off = cut.off, # cut-off used
            descriptives = descriptives # descriptives of wTO matrix
          )
        )
        
      }else{
        
        # No data -- network should not be reduced
        results <- list(
          redundant = redundant_variables,
          network = network,
          wto = list(
            matrix = wto_output, # wTO matrix
            cut_off = cut.off, # cut-off used
            descriptives = descriptives # descriptives of wTO matrix
          )
        )
        
      }
      
      # Check for "remove" method
      if(reduce.method == "remove"){
        
        # For "remove", add `keep` and `remove`
        results$keep_remove <- list(
          keep = keep, remove = remove
        )
        
      }
      
      # Ensure "UVA" class
      class(results) <- "UVA"
      
      # Return results
      return(results)
  
    }else{
      
      # Enter into user procedure...
      
      
      
      
      
      
      
      
      # User guided procedure...
      ## Re-direct to older version (for now...)
      return(
        do.call( # Perform legacy UVA
          what = "oldUVA", # call old UVA function
          args = as.list( # force list into call
            legacy_UVA( # grab input from function calls
              data, n, key, cut.off, reduce,
              reduce.method, auto, label.latent,
              FUN.args = list(...) # any other lingering arguments
            )
          )
        )
      )
    
    }
    
  }else{
    
    # Check for key
    if(!is.null(key)){
      
      ## Assign names from key
      redundant_variables <- assign_redundancy_names(
        redundant_variables = redundant_variables,
        name_vector = key
      )
      
    }else{
      
      ## Assign names from data
      redundant_variables <- assign_redundancy_names(
        redundant_variables = redundant_variables,
        name_vector = colnames(data)
      )
      
    }
    
    # If no reduction, then return results
    results <- list(
      redundant = redundant_variables,
      network = network,
      wto = list(
        matrix = wto_output, # wTO matrix
        cut_off = cut.off, # cut-off used
        descriptives = descriptives # descriptives of wTO matrix
      )
    )
    
    # Ensure "UVA" class
    class(results) <- "UVA"
    
    # Return results
    return(results)
    
  }
  
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
# auto = TRUE; label.latent = FALSE;
# EGAnet.version = packageVersion("EGAnet");
# ega_ARGS <- ega_arguments(arguments = list())
# ... # `EGA`, {lavaan}, and "EGAnet.version" arguments


#' @noRd
# Obtain UVA cut-off values ----
# Updated 07.02.2023
wto_cut_off <- function(wto_output, cut_off)
{
  
  # Obtain indices that meet cut-off
  wto_indices <- which(wto_output >= cut_off, arr.ind = TRUE)
  
  # Check whether there are values
  if(nrow(wto_indices) > 0){
    
    # Obtain unique pairs only
    wto_indices <- wto_indices[wto_indices[,"row"] < wto_indices[,"col"],]
   
    # Ensure matrix
    wto_indices <- matrix(
      wto_indices,
      ncol = 2
    )
    
    # Add column names
    colnames(wto_indices) <- c("row", "col")
     
  }
  
  # Return indices
  return(wto_indices)
  
}

#' @noRd
# Obtain descriptives ----
# Updated 02.02.2023
wto_descriptives <- function(wto_output, key = NULL){
  
  # Obtain node names
  if(!is.null(key)){
    node_names <- key
  }else{
    node_names <- colnames(wto_output)
  }
  
  # Initialize data frame
  wto_long <- data.frame(
    node_i = rep(1:ncol(wto_output), each = ncol(wto_output)),
    node_j = rep(1:ncol(wto_output), times = ncol(wto_output)),
    wto = as.vector(as.matrix(wto_output))
  )
  
  # Subset to remove duplicates
  wto_long <- wto_long[wto_long$node_i < wto_long$node_j,]
  
  # Remove all values below zero
  wto_long <- wto_long[
    wto_long$wto != 0,
  ]
  
  # Replace node names
  wto_long$node_i <- node_names[wto_long$node_i]
  wto_long$node_j <- node_names[wto_long$node_j]
  
  # Compute MAD, RANGE, QUANTILE
  MAD <- mad(wto_long$wto, constant = 1, na.rm = TRUE)
  RANGE <- range(wto_long$wto, na.rm = TRUE)
  QUANTILE <- quantile(wto_long$wto, probs = c(0.975, 0.995))
  names(QUANTILE) = c("95%", "99%")
  
  # Compute summary statistics (rounded)
  summary_statistics <- round(
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
  )
  
  # Order long data frame
  wto_long <- wto_long[order(wto_long$wto, decreasing = TRUE),]
  
  # Compute standard deviation from mean (rounded)
  wto_long$sd_from_mean <- round(
    (wto_long$wto - summary_statistics["mean"]) / summary_statistics["sd"], 3
  )
  
  # Compute MAD from median (rounded)
  wto_long$mad_from_median <- round(
    (wto_long$wto - summary_statistics["median"]) / summary_statistics["mad"], 3
  )
  
  # Round wTO values
  wto_long$wto <- round(wto_long$wto, 3)
  
  # Return list
  return(
    list(
      basic = summary_statistics,
      centralTendency = wto_long
    )
  )
  
}

#' @noRd
# Condense into index list ----
# Updated 02.02.2023
wto_index_list <- function(wto_output, wto_indices)
{
  
  # Obtain descending order of frequencies of each index
  index_frequencies <- sort(table(wto_indices), decreasing = TRUE)
  
  # Initialize index wTO sums
  index_wto_sums <- numeric(length(index_frequencies))
  
  # Obtain sums
  for(i in seq_along(index_frequencies)){
    
    # Target index
    target_index <- as.numeric(names(index_frequencies)[i])
    
    # Obtain overlapping indices
    overlapping_indices <- overlap_indices(
      wto_indices = wto_indices, target_index = target_index
    )
    
    # Obtain sum
    index_wto_sums[i] <- sum(wto_output[target_index, overlapping_indices], na.rm = TRUE)
    
  }
  
  # Add names to sums
  names(index_wto_sums) <- names(index_frequencies)
  
  # Obtain unique frequenices (should already be descending)
  unique_frequencies <- unique(index_frequencies)
  
  # Loop order unique frequencies
  ordered_frequencies <- unlist(
    lapply(unique_frequencies, function(freq){
      
      # Obtain target indices
      target_indices <- index_frequencies[index_frequencies == freq]
      
      # Obtain target sums
      target_sums <- index_wto_sums[names(target_indices)]
      
      # Reorder based on sums
      ordered_nodes <- names(target_sums)[
        order(target_sums, decreasing = TRUE)
      ]
      
      # Return as numeric
      return(as.numeric(ordered_nodes))
      
    })
  )
  
  
  # Based on ordered frequencies, extract 
  # all indices that are tied to each other
  overlap_list <- lapply(
    ordered_frequencies,
    FUN = overlap_indices,
    wto_indices = wto_indices
  )
  
  # Assign names
  names(overlap_list) <- ordered_frequencies
  
  # Return list
  return(overlap_list)
  
}

#' @noRd
# Obtain indices that overlap ----
# Updated 02.02.2023
overlap_indices <- function(wto_indices, target_index)
{
  
  # Obtain overlap indices
  overlapping_indices <- wto_indices[
    wto_indices[,"row"] == target_index |
      wto_indices[,"col"] == target_index,
  ]
  
  # Obtain vector of overlap indices
  overlapping_indices <- as.vector(as.matrix(overlapping_indices))
  
  # Remove target index
  overlapping_indices <- overlapping_indices[overlapping_indices != target_index]
  
  # Return vector
  return(overlapping_indices)
  
}

#' @noRd
# Condense the overlapping list ----
# Updated 02.02.2023
condense_overlap <- function(overlapping_list)
{
  
  # Convert all values to character
  character_list <- lapply(overlapping_list, as.character)
  
  # Initialize current index
  current_index <- 1
  
  # Remove name of list (target node) from overlap indices
  while(current_index < length(character_list)){
    
    # Obtain current target node
    target_node <- names(character_list)[current_index]
    
    # Update character list
    character_list <- condense_target_node(
      character_list, target_node
    )
    
    # Remove empty lists
    character_list <- remove_empty_elements(character_list)
    
    # Obtain current element node(s)
    element_node <- character_list[[current_index]]
    
    # Update character list
    character_list <- condense_element_node(
      character_list, target_node, element_node
    )
    
    # Remove empty lists
    character_list <- remove_empty_elements(character_list)
    
    # Update current index
    current_index <- current_index + 1
    
  }
  
  # Convert back to numeric values
  numeric_list <- lapply(character_list, as.numeric)
  
  # Return list
  return(numeric_list)
  
}

#' @noRd
# Condense target node ----
# Updated 02.02.2023
condense_target_node <- function(character_list, target_node)
{
  
  # Determine whether target node is indexed
  ## Update character list
  updated_list <- lapply(character_list, function(x){
    
    # Check whether target node is in the list
    if(target_node %in% x){
      x <- x[x != target_node]
    }
    
    # Return value
    return(x)
    
  })
  
  # Return updated list
  return(updated_list)

  
}

#' @noRd
# Condense element node ----
# Updated 02.02.2023
condense_element_node <- function(character_list, target_node, element_node)
{
  
  # Obtain target lists
  target_lists <- !names(character_list) %in% target_node
  
  # Determine whether element node is indexed
  for(node in element_node){
    
    ## Update character list (element-wise)
    character_list[target_lists] <- lapply(character_list[target_lists], function(x){
      
      # Check whether element node is in the list
      if(node %in% x){
        x <- x[x != node]
      }
      
      # Return value
      return(x)
      
    })
    
    ## Update character list (name-wise)
    if(node %in% names(character_list)){
      
      # Set list name to NULL
      character_list[
        names(character_list) %in% node
      ] <- NULL
      
    }
    
  }
  
  # Return kept elements
  return(character_list)
  
}


#' @noRd
# Remove empty elements ----
# Updated 02.02.2023
remove_empty_elements <- function(object)
{
  
  # Loop over object to find elements to keep
  keep_elements <- unlist(
    lapply(object, function(x){
      
      # Check whether values exist
      if(is.null(x)){
        return(FALSE)
      }else if(length(x) == 0){
        return(FALSE)
      }else{
        return(TRUE)
      }
      
    })
  )
  
  # Return kept elements
  return(object[keep_elements])
  
}

#' @noRd
# Assign names to redundancy list ----
# Updated 02.02.2023
assign_redundancy_names <- function(redundant_variables, name_vector)
{
  
  # Re-assign names
  ## List names
  names(redundant_variables) <- name_vector[
    as.numeric(names(redundant_variables))
  ]
  
  ## Element names
  redundant_variables <- lapply(redundant_variables, function(x){
    name_vector[x]
  })
  
  # Return list
  return(redundant_variables)
  
  
}

#' @noRd
# Obtain redundant variables ----
# Updated 02.02.2023
obtain_redundant_variables <- function(
    redundant_variables, index
)
{
  
  # Check for whether all indices
  if(index == "all"){
    
    # Obtain named node
    named_node <- as.numeric(names(redundant_variables))
    
    # Obtain element node(s)
    element_node <- unname(unlist(redundant_variables))
    
    # Obtain all nodes
    all_nodes <- c(named_node, element_node)
    
  }else{
    
    # Obtain named node
    named_node <- as.numeric(names(redundant_variables)[index])
    
    # Obtain element node(s)
    element_node <- redundant_variables[[index]]
    
    # Obtain all nodes
    all_nodes <- c(named_node, element_node)
    
  }
  
  # Return all nodes
  return(all_nodes)
  
}

#' @noRd
# Create composite variable matrix ----
# Updated 02.02.2023
create_composite_variable_matrix <- function(data, redundant_variables)
{
  
  # Create matrix for composite variables
  composite_variables <- matrix(
    NA, nrow = nrow(data),
    ncol = length(redundant_variables)
  )
  
  # Name composite variables
  colnames(composite_variables) <- paste0(
    "CV", formatC(
      1:ncol(composite_variables),
      digits = digits(ncol(composite_variables)) - 1,
      format = "d", flag = "0"
    )
  )
  
  # Return matrix
  return(composite_variables)
  
}

#' @noRd
# "Latent" reduce method ----
# Updated 02.02.2023
reduce_latent <- function(
    data,
    wto_output, 
    redundant_variables,
    lavaan_ARGS
)
{
  
  # Create matrix for composite variables
  composite_variables <- create_composite_variable_matrix(
    data = data, redundant_variables = redundant_variables
  )
  
  # Obtain CFA function
  lavaan_cfa <- lavaan::cfa
  
  # Original {lavaan} arguments
  original_lavaan_ARGS <- lavaan_ARGS
  
  # Loop over redundant variables
  for(i in seq_along(redundant_variables)){
    
    # Re-set {lavaan} arguments (unnecessary but to be sure)
    lavaan_ARGS <- original_lavaan_ARGS
    
    # Obtain all nodes
    all_nodes <- obtain_redundant_variables(
      redundant_variables = redundant_variables,
      index = i
    )
    
    # Obtain variable names for nodes
    variable_names <- colnames(data)[all_nodes]
    
    # Make CFA model
    ## (see "helper-functions.R" for function)
    model <- make_unidimensional_cfa(variable_names)
    
    # Replace "model" and "data" arguments
    lavaan_ARGS$model <- model
    lavaan_ARGS$data <- data[,all_nodes]
    
    # Obtain estimator arguments
    ## (see "helper-functions.R" for function)
    lavaan_ARGS <- estimator_arguments(lavaan_ARGS)
    
    # Estimate latent variable model
    cfa_estimate <- do.call(
      what = "lavaan_cfa",
      args = lavaan_ARGS
    )
    
    # Identify available cases
    available_cases <- lavaan::inspect(
      cfa_estimate, what = "case.idx"
    )
    
    # Obtain latent score
    latent_score <- lavaan::lavPredict(cfa_estimate)
    
    # Obtain corresponding cases for scores
    corresponding_cases <- intersect(1:nrow(data), available_cases)
    
    # Compute new composite
    composite_variables[corresponding_cases, i] <- latent_score 
    
  }
  
  # Obtain all redundant nodes
  all_nodes <- obtain_redundant_variables(
    redundant_variables = redundant_variables,
    index = "all"
  )
  
  # Remove redundant nodes from data
  new_data <- data[,-all_nodes]
  
  # Add composite variables to new data
  new_data <- cbind(
    new_data, composite_variables
  )
  
  # Return new data
  return(new_data)
  
}

#' @noRd
# "Mean" reduce method ----
# Updated 02.02.2023
reduce_mean <- function(
    data,
    wto_output, # wto_output is not actually used
    # used as placeholder for other `reduce_remove` function
    redundant_variables
)
{
  
  # Create matrix for composite variables
  composite_variables <- create_composite_variable_matrix(
    data = data, redundant_variables = redundant_variables
  )
  
  # Loop over redundant variables
  for(i in seq_along(redundant_variables)){
    
    # Obtain all nodes
    all_nodes <- obtain_redundant_variables(
      redundant_variables = redundant_variables,
      index = i
    )
    
    # Compute new composite
    composite_variables[,i] <- rowMeans(
      data[,all_nodes], na.rm = TRUE
    )
    
  }
  
  # Obtain all redundant nodes
  all_nodes <- obtain_redundant_variables(
    redundant_variables = redundant_variables,
    index = "all"
  )
  
  # Remove redundant nodes from data
  new_data <- data[,-all_nodes]
  
  # Add composite variables to new data
  new_data <- cbind(
    new_data, composite_variables
  )
  
  # Return new data
  return(new_data)
  
}

#' @noRd
# "Remove" reduce method ----
# Updated 02.02.2023
reduce_remove <- function(
    data, # data is not actually used
    # used as placeholder for other `reduce_*` functions
    wto_output, 
    redundant_variables
)
{
  
  # Initialize lists: keep and remove
  remove <- keep <- list()
  
  # Loop over redundant variables to determine
  # which variables to keep and remove
  for(i in seq_along(redundant_variables)){
    
    # Obtain all nodes
    all_nodes <- obtain_redundant_variables(
      redundant_variables = redundant_variables,
      index = i
    )
    
    # Determine whether to use wTO or standard deviation
    if(length(all_nodes) > 2){
      
      # Selection index based on maximum average wTO value
      # to other redundant variables
      selection_index <- which.max(
        colMeans(wto_output[all_nodes, all_nodes], na.rm = TRUE)
      )
      
    }else{ # Only two nodes
      
      # Selection index based on lowest maximum wTO value
      # to all other variables
      selection_index <- which.min(
        apply(wto_output[all_nodes,-all_nodes], 1, max, na.rm = TRUE)
      )
    
    }
    
    # Determine which to keep
    keep[[i]] <- all_nodes[selection_index]
    
    # Determine which to remove
    remove[[i]] <- all_nodes[setdiff(seq_along(all_nodes), selection_index)]
    
  }
  
  # Unlist lists
  keep <- unlist(keep)
  remove <- unlist(remove)
  
  # Set up results
  results <- list(
    keep = sort(keep),
    remove = sort(remove)
  )
  
  # Return results
  return(results)
  
}

#' @noRd
# "Sum" reduce method ----
# Updated 02.02.2023
reduce_sum <- function(
    data,
    wto_output, # wto_output is not actually used
    # used as placeholder for other `reduce_remove` function
    redundant_variables
)
{
  
  # Create matrix for composite variables
  composite_variables <- create_composite_variable_matrix(
    data = data, redundant_variables = redundant_variables
  )
  
  # Loop over redundant variables
  for(i in seq_along(redundant_variables)){
    
    # Obtain all nodes
    all_nodes <- obtain_redundant_variables(
      redundant_variables = redundant_variables,
      index = i
    )
    
    # Compute new composite
    composite_variables[,i] <- rowSums(
      data[,all_nodes], na.rm = TRUE
    )
    
  }
  
  # Obtain all redundant nodes
  all_nodes <- obtain_redundant_variables(
    redundant_variables = redundant_variables,
    index = "all"
  )
  
  # Remove redundant nodes from data
  new_data <- data[,-all_nodes]
  
  # Add composite variables to new data
  new_data <- cbind(
    new_data, composite_variables
  )
  
  # Return new data
  return(new_data)
  
}

#' @noRd
# Legacy arguments for "oldUVA.R" ----
# Updated 02.02.2023
legacy_UVA <- function(
    data, n, key, cut.off, reduce,
    reduce.method, auto, label.latent,
    FUN.args
)
{
  
  # Old UVA arguments
  oldUVA.args <- obtain.arguments(
    FUN = oldUVA, FUN.args = FUN.args # Need to set up defaults
  )
  
  # Replace necessary arguments
  oldUVA.args$data <- data; oldUVA.args$n <- n; oldUVA.args$key <- key;
  oldUVA.args$sig <- cut.off; oldUVA.args$reduce <- reduce;
  oldUVA.args$reduce.method <- reduce.method; oldUVA.args$auto <- auto;
  oldUVA.args$label.latent <- label.latent
  
  # Return arguments
  return(oldUVA.args)
  
}
  
#' @noRd
# Format redundant list into matrix ----
# Updated 03.02.2023
list_to_matrix <- function(redundant_list)
{
  
  # Find maximum length
  lengths <- sapply(redundant_list, length)
  
  # Initialize matrix
  redundant_matrix <- matrix(
    "", nrow = length(redundant_list),
    ncol = max(lengths, na.rm = TRUE) + 1
  ) 
  
  # Populate matrix
  for(i in seq_along(redundant_list)){
    
    # Obtain values
    values <- c(
      as.numeric(names(redundant_list)[i]),
      unname(unlist(redundant_list[[i]]))
    )
    
    redundant_matrix[i,1:length(values)] <- values
    
  }
  
  # Return redundant matrix
  return(redundant_matrix)
  
}

#' @noRd
# Convert indices to names ----
# Updated 03.02.2023
indices_to_names <- function(wto_names, redundant_matrix)
{
  
  # Loop over matrix
  for(i in 1:nrow(redundant_matrix)){
    
    # Replace indices with names
    numeric_indices <- as.numeric(redundant_matrix[i,])
    
    # Convert to names
    named_indices <- wto_names[numeric_indices]
    
    # Convert `NA` to `""`
    named_indices <- ifelse(
      is.na(named_indices), "", named_indices
    )
    
    # Replace in matrix
    redundant_matrix[i,] <- named_indices
    
  }
  
  # Return redundant matrix
  return(redundant_matrix)
  
}
  
#' @export
# S3Method for `summary()` ----
# Updated 03.02.2023
summary.UVA <- function(object, ...)
{
  
  # Obtain wTO matrix
  wto_matrix <- object$wto$descriptives$centralTendency
  
  # Obtain wTO > 30
  wto_30_named <- wto_matrix[
    wto_matrix$wto > 0.30,
    c("node_i", "node_j", "wto")
  ]
  
  # Obtain wTO > 25
  wto_25_named <- wto_matrix[
    wto_matrix$wto < 0.30 &
      wto_matrix$wto > 0.25,
    c("node_i", "node_j", "wto")
  ]
  
  # Obtain wTO > 20
  wto_20_named <- wto_matrix[
    wto_matrix$wto < 0.25 &
      wto_matrix$wto > 0.20,
    c("node_i", "node_j", "wto")
  ]
  
  # # Obtain indices with wTO greater than 0.20
  # wto_indices <- which(wto_matrix > 0.20, arr.ind = TRUE)
  # 
  # # Limit indices to unique elements
  # wto_elements <- wto_indices[
  #   wto_indices[,"row"] < wto_indices[,"col"],
  # ]
  # 
  # # Obtain values
  # values <- numeric(nrow(wto_elements))
  # 
  # # Loop over values
  # for(i in 1:nrow(wto_elements)){
  #   values[i] <- wto_matrix[
  #     wto_elements[i,"row"],
  #     wto_elements[i,"col"]
  #   ]
  # }
  # 
  # # Add values to elements
  # wto_sparse <- cbind(wto_elements, values)
  # 
  # # Order by values
  # wto_ordered <- wto_sparse[
  #   order(wto_sparse[,"values"], decreasing = TRUE),
  # ]
  # 
  # # Cut into chunks
  # ## 0.30
  # wto_30 <- wto_ordered[
  #   wto_ordered[,"values"] > 0.30,
  # ]
  # ## 0.25
  # wto_25 <- wto_ordered[
  #   wto_ordered[,"values"] > 0.25 &
  #     wto_ordered[,"values"] < 0.30,
  # ]
  # ## 0.20
  # wto_20 <- wto_ordered[
  #   wto_ordered[,"values"] > 0.20 &
  #     wto_ordered[,"values"] < 0.25,
  # ]
  # 
  # # Remove values and condense into lists
  # ## 0.30
  # wto_30_list <- condense_overlap(
  #   overlapping_list = wto_index_list(
  #     wto_output = wto_matrix,
  #     wto_indices = wto_30[
  #       , -which(colnames(wto_30) == "values")
  #     ]
  #   )
  # )
  # ## 0.25
  # wto_25_list <- condense_overlap(
  #   overlapping_list = wto_index_list(
  #     wto_output = wto_matrix,
  #     wto_indices = wto_25[
  #       , -which(colnames(wto_25) == "values")
  #     ]
  #   )
  # )
  # ## 0.20
  # wto_20_list <- condense_overlap(
  #   overlapping_list = wto_index_list(
  #     wto_output = wto_matrix,
  #     wto_indices = wto_20[
  #       , -which(colnames(wto_20) == "values")
  #     ]
  #   )
  # )
  # 
  # # Format into matrices
  # wto_30_matrix <- list_to_matrix(wto_30_list)
  # wto_25_matrix <- list_to_matrix(wto_25_list)
  # wto_20_matrix <- list_to_matrix(wto_20_list)
  # 
  # # Replace indices with variable names
  # wto_30_named <- indices_to_names(
  #   colnames(wto_matrix), wto_30_matrix
  # )
  # wto_25_named <- indices_to_names(
  #   colnames(wto_matrix), wto_25_matrix
  # )
  # wto_20_named <- indices_to_names(
  #   colnames(wto_matrix), wto_20_matrix
  # )
  
  # Prepare for print
  ## Print 0.30
  cat("Variable pairs with wTO > 0.30 (large-to-very large redundancy)\n\n")
  print(wto_30_named, quote = FALSE, row.names = FALSE)
  # no_name_print(wto_30_named)
  ## Print 0.25
  cat("\n----\n")
  cat("\nVariable pairs with wTO > 0.25 (moderate-to-large redundancy)\n\n")
  print(wto_25_named, quote = FALSE, row.names = FALSE)
  # no_name_print(wto_25_named)
  ## Print 0.20
  cat("\n----\n")
  cat("\nVariable pairs with wTO > 0.20 (small-to-moderate redundancy)\n\n")
  print(wto_20_named, quote = FALSE, row.names = FALSE)
  # no_name_print(wto_20_named)
  
}  
  
#' @export
# S3Method for `print()` ----
# Updated 03.02.2023
print.UVA <- function(x, ...)
{
  
  # Obtain wTO matrix
  wto_matrix <- x$wto$descriptives$centralTendency
  
  # Obtain wTO > 30
  wto_30_named <- wto_matrix[
    wto_matrix$wto > 0.30,
    c("node_i", "node_j", "wto")
  ]
  
  # Obtain wTO > 25
  wto_25_named <- wto_matrix[
    wto_matrix$wto < 0.30 &
      wto_matrix$wto > 0.25,
    c("node_i", "node_j", "wto")
  ]
  
  # Obtain wTO > 20
  wto_20_named <- wto_matrix[
    wto_matrix$wto < 0.25 &
      wto_matrix$wto > 0.20,
    c("node_i", "node_j", "wto")
  ]
  
  # Prepare for print
  ## Print 0.30
  cat("Variable pairs with wTO > 0.30 (large-to-very large redundancy)\n\n")
  print(wto_30_named, quote = FALSE, row.names = FALSE)
  # no_name_print(wto_30_named)
  ## Print 0.25
  cat("\n----\n")
  cat("\nVariable pairs with wTO > 0.25 (moderate-to-large redundancy)\n\n")
  print(wto_25_named, quote = FALSE, row.names = FALSE)
  # no_name_print(wto_25_named)
  ## Print 0.20
  cat("\n----\n")
  cat("\nVariable pairs with wTO > 0.20 (small-to-moderate redundancy)\n\n")
  print(wto_20_named, quote = FALSE, row.names = FALSE)
  # no_name_print(wto_20_named)
  
}  
  
  
  











