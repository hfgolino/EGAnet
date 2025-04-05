#' @title Diagnostics Analysis for Low Stability Items
#'
#' @description Computes the between- and within-community
#' \code{strength} of each variable for each community
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
#' \item{uva}{Output from \code{\link[EGAnet]{UVA}}}
#'
#' \item{minor}{A list containing suggested items to \code{keep}, \code{remove}, and
#' a matrix for probable minor dimensions (\code{minor.matrix})}
#'
#' \item{loadings}{Output from \code{\link[EGAnet]{net.loads}}}
#'
#' @examples
#' # Load data
#' wmt <- wmt2[,7:24]
#'
#' \dontrun{
#' # Obtain diagnostics
#' diagnostics <- itemDiagnostics(wmt, ncores = 2)}
#'
#' @author Alexander P. Christensen <alexpaulchristensen@gmail.com> and Hudson Golino <hfg9s at virginia.edu>
#'
#' @export
#'
# Item diagnostics ----
# Updated 05.04.2025
itemDiagnostics <- function(data, ...)
{

  # Check for errors
  data <- itemDiagnostics_errors(data, ...)

  # Get node names
  node_names <- dimnames(data)[[2]]

  # Send message
  message("Performing bootstrap...")

  # Perform bootEGA
  boot <- bootEGA(data, ...)

  # Obtain item stabilities
  stabilities <- boot$stability$item.stability$item.stability$empirical.dimensions

  # Select for low stabilities
  low_stabilities <- stabilities[stabilities < 0.75]
  low_names <- names(low_stabilities)

  # Check for good stability
  if(length(low_stabilities) == 0){
    message <- "All items have good stability (>= 0.75)"
    cat(message)
    return(message)
  }

  # Perform UVA (for wTO)
  uva <- UVA(data, reduce = FALSE, ...)

  # Check for low dependence
  local_dependence <- low_names[
    colSums(uva$wto$matrix[low_names, low_names] >= 0.25) > 0
  ]

  # Check for minor dimensions
  minor <- minor_dimensions(
    wto_output = uva$wto$matrix,
    stabilities = boot$stability$item.stability$item.stability$all.dimensions[low_names,],
    cut_off = 0.95
  )

  # Obtain loadings
  loadings <- silent_call(net.loads(boot$EGA, ...))
  loadings_unstable <- loadings$std[low_names,, drop = FALSE]

  # Check for multidimensional
  multidimensional <- low_names[rowSums(abs(loadings_unstable) >= 0.20) > 1]

  # Check for low loadings
  low_loadings <- low_names[
    lvapply(low_names, function(name){loadings_unstable[name, boot$EGA$wc[name]] < 0.20})
  ]

  # Set up poor stability diagnostics
  diagnostic_df <- data.frame(node = low_names)
  diagnostic_df$diagnostic <- swiftelse(low_names %in% local_dependence, "LD ", "")
  diagnostic_df$diagnostic <- paste0(
    diagnostic_df$diagnostic, swiftelse(
      low_names %in% as.matrix(minor$minor.matrix), "MiD ", ""
    )
  )
  diagnostic_df$diagnostic <- paste0(
    diagnostic_df$diagnostic, swiftelse(low_names %in% multidimensional, "MuD ", "")
  )
  diagnostic_df$diagnostic <- paste0(
    diagnostic_df$diagnostic, swiftelse(low_names %in% low_loadings, "Low ", "")
  )

  # Set names for diagnostics
  row.names(diagnostic_df) <- low_names
  diagnostic_df$node <- NULL

  # Set results
  results <- list(
    diagnostics = diagnostic_df,
    uva = uva, minor = minor,
    loadings = loadings$std
  )

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
# Updated 05.04.2025
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
# Updated 05.04.2025
minor_dimensions <- function(wto_output, stabilities, cut_off = 0.95)
{

  # Transpose stabilities
  stabilities <- t(stabilities)

  # Obtain cosines
  cosines <- EGAnet::cosine(stabilities)
  key <- colnames(cosines)

  # Identify cosines greater than cut-off
  index <- which(cosines >= cut_off, arr.ind = TRUE)
  index <- index[index[,1] < index[,2],, drop = FALSE]
  dimnames(index)[[2]] <- c("node_i", "node_j")

  # Get redundancy list
  redundant_variables <- get_redundancy_list(cosines, index)

  # Get lengths
  lengths <- nvapply(redundant_variables, length)

  # Get max length
  max_length <- max(lengths)

  # Order from least to most
  redundant_variables <- redundant_variables[order(lengths)]

  # Length of lengths
  n_lengths <- length(lengths)

  # Create matrix
  minor_matrix <- matrix(0, nrow = n_lengths, ncol = max_length + 1)

  # Populate matrix
  for(i in seq_len(n_lengths)){

    # Get redundant variables
    variables <- obtain_redundant_variables(redundant_variables, i)

    # Populate
    minor_matrix[i, seq_along(variables)] <- variables

  }

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
    keep_index <- rowSums(minor_matrix != 0) > 1
    redundant_variables <- redundant_variables[keep_index]
    minor_matrix <- minor_matrix[keep_index,, drop = FALSE]

  }

  # Loop over redundant variables and return keep and remove
  selection_list <- lapply(seq_along(redundant_variables), function(index){

    # Obtain all nodes
    all_nodes <- obtain_redundant_variables(redundant_variables, index)

    # Determine whether to use wTO or standard deviation
    if(length(all_nodes) > 2){

      # Selection index based on maximum average
      # wTO value to other redundant variables
      selection_index <- which.min(
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

  # Set up keys
  minor_matrix[] <- apply(minor_matrix, 2, function(x){key[x]})
  minor_matrix <- as.data.frame(minor_matrix)
  colnames(minor_matrix) <- NULL

  # Return results
  return(
    list(
      keep = key[ulapply(selection_list, function(x){x$keep})],
      remove = key[ulapply(selection_list, function(x){x$remove})],
      minor.matrix = minor_matrix
    )
  )

}
