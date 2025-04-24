#' @title Diagnostics Analysis for Low Stability Items
#'
#' @description Computes item diagnostics to determine whether there exist several
#' potential psychometric issues that affect \code{\link[EGAnet]{itemStability}}.
#' These issues include local dependence (using \code{\link[EGAnet]{UVA}}),
#' minor dimensions (using \code{\link[EGAnet]{cosine}} on the stability patterns),
#' multidimensional items (using \code{\link[EGAnet]{net.loads}}), and
#' items with low loadings (using \code{\link[EGAnet]{net.loads}})
#'
#' @param data Matrix or data frame.
#' Should consist only of variables to be used in the analysis
#'
#' @param ... Additional arguments to pass on to
#' \code{\link[EGAnet]{bootEGA}},
#' \code{\link[EGAnet]{net.loads}}, and
#' \code{\link[EGAnet]{UVA}}
#'
#' @return Returns a list containing:
#'
#' \item{diagnostics}{A data frame containing the diagnostics of low item stabilities
#' (see \code{\link[EGAnet]{itemStability}})}
#'
#' \item{boot}{Output from \code{\link[EGAnet]{bootEGA}}}
#'
#' \item{uva}{Output from \code{\link[EGAnet]{UVA}}}
#'
#' \item{minor}{A list containing suggested items to \code{keep}, \code{remove}, and
#' a matrix for probable minor dimensions (\code{minor.matrix})}
#'
#' \item{loadings}{Output from \code{\link[EGAnet]{net.loads}}}
#'
#' \item{suggested}{Variables that are suggested to be retained to increase item stability}
#'
#' @examples
#' # Load data
#' wmt <- wmt2[,7:24]
#'
#' \dontrun{
#' # Obtain diagnostics
#' diagnostics <- itemDiagnostics(wmt, ncores = 2)}
#'
#' @author Alexander P. Christensen <alexpaulchristensen@gmail.com>, Hudson Golino <hfg9s at virginia.edu>, and Luis Eduardo Garrido <garrido.luiseduardo@gmail.com>
#'
#' @export
#'
# Item diagnostics ----
# Updated 24.04.2025
itemDiagnostics <- function(data, ...)
{

  # Check for errors
  data <- itemDiagnostics_errors(data, ...)

  # Send message
  message("Performing bootstrap...")

  # Perform bootEGA
  boot <- bootEGA(data, ...)

  # Get node names
  node_names <- dimnames(boot$EGA$network)[[2]]

  # Identify two node communities
  ## Get node counts
  node_counts <- fast_table(boot$EGA$wc)

  ## Check for two node communities
  communities <- names(node_counts[node_counts == 2])

  ## Set up data frame
  two_node <- as.data.frame(
    do.call(rbind, lapply(
      communities, function(community){
        names(boot$EGA$wc)[boot$EGA$wc == community]
      }
    ))
  )

  ## Set row names to communities and remove column names
  row.names(two_node) <- communities; colnames(two_node) <- NULL

  # Obtain item stabilities
  stabilities <- boot$stability$item.stability$item.stability$empirical.dimensions

  # Select for low stabilities
  low_stabilities <- stabilities[stabilities < 0.75]
  low_names <- names(low_stabilities)
  n_low_stabilities <- length(low_stabilities)

  # Check for good stability
  if(n_low_stabilities < 2){

    # Catch ellipse
    ellipse <- list(...)

    # Set up message
    message <- "All items have good stability (>= 0.75) "

    # Add except message
    if(n_low_stabilities == 1){
      message <- paste0(
        message, "except for '",
        low_names, "'.\nRemoval is suggested\n")
    }

    # Check for verbose
    if("verbose" %in% names(ellipse) && ellipse$verbose){
      cat(message)
    }else{
      cat(message)
    }

    # Return shell of results
    return(
      list(
        boot = boot, uva = UVA(data, reduce = FALSE, ...),
        loadings = silent_call(net.loads(boot$EGA, ...)),
        keep = swiftelse(
          n_low_stabilities == 1,
          node_names[node_names != low_names],
          node_names
        )
      )
    )

  }

  # Perform UVA (for wTO)
  uva <- UVA(data, reduce = FALSE, ...)

  # Check for local dependence
  local_dependence <- low_names[
    colSums(uva$wto$matrix[low_names, low_names, drop = FALSE] >= 0.25) > 0
  ]

  # Check for minor dimensions
  minor <- minor_dimensions(
    ega = boot$EGA, wto_output = uva$wto$matrix,
    stabilities = boot$stability$item.stability$item.stability$all.dimensions[low_names,, drop = FALSE],
    cut_off = 0.95
  )

  # Obtain loadings
  loadings <- silent_call(net.loads(boot$EGA, ...))

  # Obtain unstable loadings and make them absolute for the following checks
  loadings_unstable <- abs(loadings$std[low_names,, drop = FALSE])

  # Check for multidimensional
  multidimensional <- low_names[rowSums(loadings_unstable >= 0.20) > 1]

  # Check for low loadings
  low_loadings <- low_names[
    lvapply(low_names, function(name){loadings_unstable[name, boot$EGA$wc[name]] < 0.20})
  ]

  # Set up poor stability diagnostics
  diagnostic_df <- data.frame(node = low_names)
  diagnostic_df$diagnostic <- swiftelse(low_names %in% local_dependence, "LD ", "")
  diagnostic_df$diagnostic <- paste0(
    diagnostic_df$diagnostic, swiftelse(low_names %in% as.matrix(minor), "MiD ", "")
  )
  diagnostic_df$diagnostic <- paste0(
    diagnostic_df$diagnostic, swiftelse(low_names %in% multidimensional, "MuD ", "")
  )
  diagnostic_df$diagnostic <- paste0(
    diagnostic_df$diagnostic, swiftelse(low_names %in% low_loadings, "Low ", "")
  )
  diagnostic_df$diagnostic <- paste0(
    diagnostic_df$diagnostic, swiftelse(low_names %in% two_node, "Two ", "")
  )

  # Set names for diagnostics
  row.names(diagnostic_df) <- low_names
  diagnostic_df$node <- NULL

  # Set results
  results <- list(
    diagnostics = diagnostic_df,
    boot = boot, uva = uva, minor = minor,
    loadings = loadings, two_node = two_node
  )

  # Add suggested results
  results$suggested <- automated_selection(results)

  # Set class
  class(results) <- "itemDiagnostics"

  # Return results
  return(results)

}

#' @noRd
# Errors ----
# Updated 05.04.2025
itemDiagnostics_errors <- function(data, ...)
{

  # 'data' errors
  object_error(data, c("matrix", "data.frame", "tibble"), "auto.correlate")

  # Check for tibble
  if(get_object_type(data) == "tibble"){
    data <- as.data.frame(data)
  }

  # Check for usable data
  if(needs_usable(list(...))){
    data <- usable_data(data, FALSE)
  }

  # Return usable data (in case of tibble)
  return(data)

}

#' @exportS3Method
# S3 Print Method ----
# Updated 08.04.2025
print.itemDiagnostics <- function(x, ...)
{

  # Print title
  cat(styletext("Low Item Stability Diagnostics", defaults = "bold"), "\n")

  # Remove column names
  colnames(x$diagnostics) <- NULL

  # Print diagnostics
  print(x$diagnostics)

  # Collect diagnostics
  diagnostics <- strsplit(trimws(as.matrix(x$diagnostics)), split = " ")

  # Add breakspace
  cat("----\n")

  # Build diagnostic codes
  start <- "Diagnostic codes: "
  text <- c(
    swiftelse(
      any(lvapply(diagnostics, function(x){"LD" %in% x})),
      "'LD' = local dependence (see `$uva`)", ""
    ),
    swiftelse(
      any(lvapply(diagnostics, function(x){"MiD" %in% x})),
      "'MiD' = minor dimension (see `$minor`)", ""
    ),
    swiftelse(
      any(lvapply(diagnostics, function(x){"MuD" %in% x})),
      "'MuD' = multidimensional (see `$loadings`)", ""
    ),
    swiftelse(
      any(lvapply(diagnostics, function(x){"Low" %in% x})),
      "'Low' = low loading (see `$loadings`)", ""
    ),
    swiftelse(
      any(lvapply(diagnostics, function(x){"Two" %in% x})),
      "'Two' = two-node community (see `$two_node`)", ""
    )
  )

  # Print codes
  cat(paste0(start, paste0(text[text != ""], collapse = ", ")))

}

#' @exportS3Method
# S3 Summary Method ----
# Updated 05.04.2025
summary.itemDiagnostics <- function(object, ...)
{
  print(object, ...) # same as print
}

#' @noRd
# Cosine for minor dimensions stabilities ----
# Updated 18.04.2025
minor_dimensions <- function(ega, wto_output, stabilities, cut_off = 0.95)
{

  # Transpose stabilities
  stabilities <- t(stabilities)

  # Obtain cosines
  cosines <- cosine(stabilities)
  key <- colnames(cosines)

  # Identify cosines greater than cut-off
  index <- which(cosines >= cut_off, arr.ind = TRUE)
  index <- index[index[,1] < index[,2],, drop = FALSE]
  dimnames(index)[[2]] <- c("node_i", "node_j")

  # Get redundancy list
  redundant_variables <- get_redundancy_list(cosines, index)

  # Get lengths
  lengths <- nvapply(redundant_variables, length)

  # Length of lengths
  n_lengths <- length(lengths)

  # Check for lengths
  if(n_lengths == 0){
    return(matrix(0, nrow = n_lengths, ncol = 0))
  }

  # Get max length
  max_length <- max(lengths)

  # Order from least to most
  redundant_variables <- redundant_variables[order(lengths)]

  # Initial columns
  minor_columns <- max_length + 1

  # Create matrix
  minor_matrix <- matrix(0, nrow = n_lengths, ncol = minor_columns)

  # Populate matrix
  for(i in seq_len(n_lengths)){

    # Get redundant variables
    variables <- obtain_redundant_variables(redundant_variables, i)

    # Populate
    minor_matrix[i, seq_along(variables)] <- variables

  }

  # # Only select rows that have 3 or less nodes
  # less_than <- rowSums(minor_matrix != 0) < 4
  # minor_matrix <- minor_matrix[
  #   less_than, seq_len(swiftelse(minor_columns < 3, minor_columns, 3)), drop = FALSE
  # ]
  # n_lengths <- sum(less_than)

  # Check for more than one
  if(n_lengths > 1){

    # Remove if variables are available in higher sets
    for(i in 1:(n_lengths - 1)){

      # Identify whether all numbers exist at a higher set
      current_set <- minor_matrix[i,]
      current_set <- current_set[current_set != 0]

      # Map current set to other sets
      map_set <- lvapply((i + 1):n_lengths, function(j){
        current_set[current_set != 0] %in% minor_matrix[j,]
      }, LENGTH = length(current_set))

      # Check sums
      minor_matrix[i, rowSums(map_set) > 0] <- 0

    }

    # Remove sets
    minor_matrix <- minor_matrix[rowSums(minor_matrix != 0) > 1,, drop = FALSE]

  }

  # Set up keys
  zeros <- minor_matrix == 0
  minor_matrix[!zeros] <- key[minor_matrix]
  minor_matrix[zeros] <- ""

  # Double-check true minor dimensions
  total <- dim(minor_matrix)[1]
  total_sequence <- seq_len(total)
  extra_dimensions <- seq.int(ega$n.dim + 1, ega$n.dim + total)

  # Create new structure with minor dimensions
  for(i in total_sequence){

    # Identify variables
    variables <- minor_matrix[i,]

    # Update EGA dimensions
    ega$wc[variables[variables != ""]] <- extra_dimensions[i]

  }

  # Compute loadings (make absolute for check below)
  loadings <- abs(silent_call(net.loads(ega)$std))

  # Numeric communities
  numeric_communities <- as.numeric(dimnames(loadings)[[2]])

  # Check over minor dimensions
  for(i in total_sequence){

    # Identify variables
    variables <- minor_matrix[i,]
    variables <- variables[variables != ""]

    # Check loadings
    greater_than <- loadings[
      variables, numeric_communities == extra_dimensions[i]
    ] > 0.20

    # Check whether to keep
    keep <- !(is.na(greater_than) | !greater_than)

    # Update minor matrix
    minor_matrix[i, seq_along(keep)] <- swiftelse(keep, variables, "")

  }

  # Ensure at least two variables
  minor_matrix <- minor_matrix[
    apply(minor_matrix, 1, function(x){sum(x != "")}) > 1,, drop = FALSE
  ]

  # Re-organize
  minor_matrix <- t(apply(minor_matrix, 1, sort, decreasing = TRUE))
  minor_matrix <- as.data.frame(minor_matrix)
  colnames(minor_matrix) <- NULL

  # Return results
  return(minor_matrix)

}

#' @noRd
# Loadings for minor dimensions stabilities ----
# Updated 07.04.2025
loadings_remove <- function(boot, stabilities, cut_off = 0.35, ...)
{

  # Compute loadings
  loadings <- silent_call(
    net.loads(boot$EGA, ...)$std[dimnames(boot$EGA$network)[[2]],, drop = FALSE]
  )

  # Get residuals from implied - empirical correlations
  residuals <- abs(nload2cor(loadings) - boot$EGA$correlation)
  diag(residuals) <- 1

  # Get node names
  node_names <- dimnames(stabilities)[[1]]

  # Get low stability residuals
  low_residuals <- residuals[node_names, node_names]

  # Estimate EGA
  ega <- EGA(
    low_residuals, n = boot$EGA$n,
    algorithm = "louvain", order = "lower",
    plot.EGA = FALSE, verbose = FALSE, ...
  )

  # Get loadings
  loadings <- silent_call(net.loads(ega, ...)$std[node_names,, drop = FALSE])

  # Get maximum loadings
  max_loadings <- nvapply(as.data.frame(abs(t(loadings))), function(x){max(x)})

  # Send labels
  return(names(max_loadings[max_loadings < cut_off]))

}

#' @noRd
# Automated node selection ----
# Updated 08.04.2025
automated_selection <- function(result)
{

  # Collect diagnostics
  diagnostics <- strsplit(result$diagnostics$diagnostic, split = " ")

  # Set tracker and diagnostics names
  tracker <- names(diagnostics) <- row.names(result$diagnostics)

  # Identify unique tags
  unique_tags <- unique(unlist(diagnostics))

  # Track diagnostics
  tracker <- names(diagnostics)

  # Set column indices
  index <- seq_along(result$boot$EGA$wc)
  names(index) <- names(result$boot$EGA$wc)

  # Set good nodes
  good_nodes <- setdiff(names(index), tracker)

  # Initialize keep list
  keep_list <- vector("list", length = length(unique_tags))
  names(keep_list) <- unique_tags

  # Community sequence
  community_sequence <- seq_len(result$boot$EGA$n.dim)

  # 1. Check for two node communities
  if("Two" %in% unique_tags){

    # Get number of two node communities
    n_two <- dim(result$two_node)[1]

    # Set keep vector
    keep <- character(length = n_two)

    # Loop over to select node to keep
    for(i in seq_len(n_two)){

      # Get index
      current_index <- index[unlist(result$two_node[i,])]

      # Selection index based on lowest maximum wTO value to all other variables
      keep[i] <- names(
        which.min(
          apply(result$uva$wto$matrix[current_index, -current_index], 1, max, na.rm = TRUE)
        )
      )

      # Remove from tracker
      tracker <- tracker[!tracker %in% names(current_index)]

    }

    # Update keep list
    keep_list$Two <- keep

  }

  # 2. Check for minor dimensions
  if("MiD" %in% unique_tags){

    # Get number of minor dimensions
    n_minor <- dim(result$minor)[1]

    # Set up keep list
    keep <- vector("list", length = n_minor)

    # Loop over to select node to keep
    for(i in seq_len(n_minor)){

      # Get index
      current_index <- as.character(result$minor[i,])

      # Check for tracker
      check_tracker <- current_index %in% tracker
      check_tracker <- current_index[check_tracker]

      # Continue with remaining variables
      if(length(check_tracker) != 0){

        # Keep the highest stability
        keep[[i]] <- check_tracker[
          which.max(
            apply(
              result$boot$stability$item.stability$item.stability$all.dimensions[
                check_tracker, community_sequence
              ], 1, max
            )
          )
        ]

        # Remove from tracker
        tracker <- tracker[!tracker %in% check_tracker]

      }

    }

    # Update keep list
    keep_list$MiD <- keep

  }

  # 3. Check for multidimensional
  if("MuD" %in% unique_tags){

    # Get multidimensional
    multidimesnional <- diagnostics[lvapply(diagnostics, function(x){"MuD" %in% x})]

    # Remove indices from tracker
    tracker <- tracker[!tracker %in% names(multidimesnional)]

  }

  # 4. Check for low loadings
  if("Low" %in% unique_tags){

    # Get low loadings
    low_loadings <- diagnostics[lvapply(diagnostics, function(x){"Low" %in% x})]

    # Remove indices from tracker
    tracker <- tracker[!tracker %in% names(low_loadings)]

  }

  # Return remaining good, tracker, and keep nodes
  return(c(good_nodes, tracker, unlist(keep_list, use.names = FALSE)))

}
