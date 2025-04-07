#' @title Diagnostics Analysis for Low Stability Items
#'
#' @description Computes the between- and within-community
#' \code{strength} of each variable for each community
#'
#' @param data Matrix or data frame.
#' Should consist only of variables to be used in the analysis
#'
#' @param minor.method Character (length = 1).
#' Method to detect the likelihood of minor dimensions.
#' Available options include:
#'
#' \itemize{
#'
#' \item \code{"cosine"} (default) --- uses \code{\link[EGAnet]{cosine}} similarity to determine
#' the similarity of the item stability patterns from \code{\link[EGAnet]{bootEGA}}.
#' More consistent patterns across multiple dimensions reflect an increased likelihood
#' of a minor dimension (particularly with respect to extra dimensions forming).
#' By default, values greater than or equal to \code{0.95} are selected
#'
#' \item \code{"residuals"} --- subtracts the network-implied zero-order correlation matrix
#' from the empirical zero-order correlation matrix to determine residuals. These residuals
#' are submitted to an \code{\link[EGAnet]{EGA}} and network loadings (\code{\link[EGAnet]{net.loads}})
#' are re-computed. By default, loadings greater than or equal to \code{0.35} are selected
#'
#' }
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
# Updated 07.04.2025
itemDiagnostics <- function(data, minor.method = c("cosine", "residuals"), ...)
{

  # Check for missing arguments (argument, default, function)
  minor.method <- set_default(minor.method, "cosine", itemDiagnostics)

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

  # Check for local dependence
  local_dependence <- low_names[
    colSums(uva$wto$matrix[low_names, low_names, drop = FALSE] >= 0.25) > 0
  ]

  # Check for minor dimensions
  minor <- switch(
    minor.method,
    "cosine" = minor_dimensions(
      wto_output = uva$wto$matrix,
      stabilities = boot$stability$item.stability$item.stability$all.dimensions[low_names,],
      cut_off = 0.95
    ),
    "residuals" = loadings_remove(
      boot = boot,
      stabilities = boot$stability$item.stability$item.stability$all.dimensions[low_names,],
      cut_off = 0.35, ...
    )
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
    diagnostic_df$diagnostic, swiftelse(low_names %in% as.matrix(minor), "MiD ", "")
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
    boot = boot, uva = uva,
    minor = minor, loadings = loadings
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
# Updated 07.04.2025
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

  # Length of lengths
  n_lengths <- length(lengths)

  # Check for lengths
  if(length(lengths) == 0){
    return(matrix(0, nrow = n_lengths, ncol = 0))
  }

  # Get max length
  max_length <- max(lengths)

  # Order from least to most
  redundant_variables <- redundant_variables[order(lengths)]

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
    minor_matrix <- minor_matrix[rowSums(minor_matrix != 0) > 1,, drop = FALSE]

  }

  # Set up keys
  minor_matrix <- as.data.frame(ifelse(minor_matrix == 0, "", key[minor_matrix]))
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
