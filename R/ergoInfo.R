#' Ergodicity Information Index
#'
#' @description Computes the Ergodicity Information Index
#'
#' @param dynEGA.object A \code{\link[EGAnet]{dynEGA.ind.pop}} object
#'
#' @param use Character.
#' A string indicating what network element will be used
#' to compute the algorithm complexity, the list of edges or the weights of the network.
#' Defaults to \code{use = "edge.list"}.
#' Current options are:
#'
#' \itemize{
#'
#' \item{\strong{\code{"edge.list"}}}
#' {Calculates the algorithm complexity using the list of edges.}
#'
#' \item{\strong{\code{"unweighted"}}}
#' {Calculates the algorithm complexity using the binary weights of the network.
#' 0 = edge absent and 1 = edge present}
#'
#' \item{\strong{\code{"weighted"}}}
#' {Calculates the algorithm complexity using the weights of the network.}
#' }
#'
#' @return Returns a list containing:
#'
#' \item{PrimeWeight}{The prime-weight encoding of the individual networks}
#'
#' \item{PrimeWeight.pop}{The prime-weight encoding of the population network}
#'
#' \item{Kcomp}{The Kolmogorov complexity of the prime-weight encoded individual networks}
#'
#' \item{Kcomp.pop}{The Kolmogorov complexity of the prime-weight encoded population network}
#'
#' \item{complexity}{The complexity metric proposed by Santora and Nicosia (2020)}
#'
#' \item{EII}{The Ergodicity Information Index}
#'
#' @author Hudson Golino <hfg9s at virginia.edu> and Alexander Christensen <alexpaulchristensen@gmail.com>
#'
#' @examples
#' # Obtain data
#' sim.dynEGA <- sim.dynEGA # bypasses CRAN checks
#'
#' \dontrun{
#' # Dynamic EGA individual and population structure
#' dyn.ega1 <- dynEGA.ind.pop(
#'   data = sim.dynEGA[,-26], n.embed = 5, tau = 1,
#'   delta = 1, id = 25, use.derivatives = 1,
#'   ncores = 2, corr = "pearson"
#' )
#'
#' # Compute empirical ergodicity information index
#' eii <- ergoInfo(
#'   dynEGA.object = dyn.ega1,
#'   use = "weighted"
#' )}
#'
#' @export
#'
# Ergodicity Information Index
# Updated 10.07.2023
ergoInfo <- function(
    dynEGA.object,
    use = c("edge.list", "unweighted")
    # , seed = 1234
)
{
  
  # Check for missing arguments (argument, default, function)
  use <- set_default(use, "edge.list", ergoInfo)
  
  # Check for appropriate class ("dynEGA.ind.pop" defunct to legacy)
  if(!is(dynEGA.object, "dynEGA") & !is(dynEGA.object, "dynEGA.ind.pop")){
    class_error(dynEGA.object, "dynEGA")
  }
  
  # Get proper objects (if not, send an error)
  dynega_objects <- get_dynEGA_object(dynEGA.object)
  
  # Get individual networks
  individual_networks <- lapply(dynega_objects$individual, function(x){x$network})
  
  # Get sequence for number of individuals
  individual_sequence <- seq_along(individual_networks)
  
  # Get adjacency networks
  adjacency_networks <- lapply(
    individual_networks, # NAs occur when there is zero variance 
    function(x){
      
      # Get adjacency
      adjacency <- x != 0 & !is.na(x)
      
      # Attach number of edges as attribute
      attr(adjacency, "edges") <- sum(adjacency) / 2
      
      # Return adjacency
      return(adjacency)
      
    }
  )
  
  # Get ANY edges across the individual networks
  edges <- symmetric_matrix_lapply(adjacency_networks, any)
  
  # "unweighted" needs canonical prime association
  if(use == "unweighted"){
  
    # Get order of each counts
    edge_ordering <- order(
      nvapply(adjacency_networks, attr, "edges"), 
      decreasing = FALSE
    )
    
    # Reorder networks and adjacencys
    individual_networks <- individual_networks[edge_ordering]
    adjacency_networks <- adjacency_networks[edge_ordering]
    
    # Get prime numbers
    prime_numbers <- get(
      data("prime.num", package = "EGAnet", envir = environment())
    )[individual_sequence]
    
    # Get prime weights
    prime_weights <- lapply(individual_sequence, function(case){
      
      # Assign primes
      prime_network <- adjacency_networks[[case]] * prime_numbers[case]
      
      # Assign 1s to 0s
      prime_network[prime_network == 0] <- 1
      
      # Return prime network
      return(prime_network)
      
    })
    
    # Get encoding matrix
    encoding_matrix <- Reduce("*", prime_weights)
    
    # Revert 1s to 0s
    encoding_matrix[encoding_matrix == 1] <- 0
    
  }
  
  # Store upper triangle indices (use for population as well)
  upper_triangle <- upper.tri(edges)
  
  # Set upper triangle to FALSE
  edges[upper_triangle] <- FALSE
  
  # Get weights
  # Use `keep_weights` for quick indexing
  if(use == "edge.list"){
    weights <- 1 # not used!
    keep_weights <- c(1L, 2L)
  }else{
    weights <- encoding_matrix[edges]
    keep_weights <- c(1L, 2L, 3L)
    # original: keep_weights <- 3L 
  }
  
  # Get edge list ("col" then "row" matches {igraph})
  edge_list <- cbind(which(edges, arr.ind = TRUE)[,c("col", "row")], weights)
  # Order matters!! (see pasting in `k_complexity`)
  
  # Get edge list rows
  edge_rows <- dim(edge_list)[1]
  
  # Get edge list sequence
  edge_sequence <- seq_len(edge_rows)
  
  # Get seeds for reproducible results
  # Includes default number of iterations (1000)
  # (defined in Santoro & Nicosia, 2020)
  # seed_values <- reproducible_seed(n = 1000, seed = seed)
  iter_sequence <- seq_len(1000)

  # Get k-complexity
  individual_kcomplexity <- nvapply( # seed_values,
    iter_sequence, function(iteration){
      
      # Return k-complexity
      return(
        k_complexity(
          edge_list[ # rows
            # reproducible_sample( # reproducible `sample`
            sample(
              x = edge_sequence, size = edge_rows,
              replace = TRUE# , seed = single_seed
            ),
            keep_weights # either pairwise edges or weights
          ]
        )
      )
      
    }
  )
  
  # Set up population
  
  # Get population adjacency
  population_edges <- dynega_objects$population$network != 0
  
  # Branch based on "use"
  if(use == "unweighted"){
    
    # Prime will always be equal 2
    population_encoding <- population_edges * 2
 
    # Revert 1s to 0s
    population_encoding[population_encoding == 1] <- 0
    
  }
  
  # Set upper triangle to FALSE
  population_edges[upper_triangle] <- FALSE
  
  # Get weights
  # `keep_weights` was defined above
  if(use == "edge.list"){
    population_weights <- 1 # not used!
  }else{
    population_weights <- population_encoding[population_edges]
  }
  
  # Get edge list ("col" then "row" matches {igraph})
  population_edge_list <- cbind(
    which(population_edges, arr.ind = TRUE)[,c("col", "row")], 
    population_weights
  )
  # Order matters!! (see pasting in `k_complexity`)
  
  # Get edge list rows
  population_edge_rows <- dim(population_edge_list)[1]
  
  # Get edge list sequence
  population_edge_sequence <- seq_len(population_edge_rows)
  
  # Get k-complexity
  population_kcomplexity <- nvapply( # seed_values,
    iter_sequence, function(single_seed){
      
      # Return k-complexity
      return(
        k_complexity(
          edge_list[ # rows
            # reproducible_sample( # reproducible `sample`
            sample(
              x = population_edge_sequence, size = population_edge_rows,
              replace = TRUE# , seed = single_seed
            ),
            keep_weights # either pairwise edges or weights
          ]
        )
      )
      
    }
  )
  
  # Pre-compute values
  mean_individual_complexity <- mean(individual_kcomplexity, na.rm = TRUE)
  mean_population_complexity <- mean(population_kcomplexity, na.rm = TRUE)
  
  # Set up results
  results <- list(
    KComp = mean_individual_complexity,
    KComp.pop = mean_population_complexity,
    EII = sqrt(dynega_objects$population$n.dim)^(
      (mean_individual_complexity / mean_population_complexity) / log(edge_rows)
    )
  )
  
  # Check for prime weights
  if(use == "unweighted"){
    results$PrimeWeight <- remove_attributes(encoding_matrix)
    results$PrimeWeight.pop <- remove_attributes(population_encoding)
  }
  
  # Add "methods" attribute
  attr(results, "methods") <- list(use = use)
  
  # Add class
  class(results) <- "EII"
  
  # Return results
  return(results)
  
  
}

# Bug checking ----
# DATA
# Population, group, and individual structure
# dynEGA.object <- dynEGA(
#   data = sim.dynEGA,
#   level = c("individual", "group", "population"),
#   ncores = 8, verbose = TRUE
# )
# use = "edge.list"; seed = 1234
# r_sample_seeds <- EGAnet:::r_sample_seeds
# r_sample_with_replacement <- EGAnet:::r_sample_with_replacement
# r_sample_without_replacement <- EGAnet:::r_sample_without_replacement
# Need above functions for testing!

#' @noRd
# k-complexity ----
# Updated 10.07.2023
k_complexity <- function(values)
{
  
  # Streamlined form
  return(
    length( # length of compression
      memCompress( # bit compression
        paste0( # bits (matches `toString`)
          t(values), collapse = ", "
        ), type = "gzip" # type of compression
      )
    )
  )
  
  # Original (first) method
  # return(
  #   length( # length of compression
  #     memCompress( # bit compression
  #       paste0( # bits (matches `toString`)
  #         values, collapse = ", "
  #       ), type = "gzip" # type of compression
  #     )
  #   )
  # )
  
}

#' @noRd
# Structural edge overlap ----
# Updated 10.07.2023
structural_overlap <- function(adjacency_networks)
{
  
  # Convert to numeric
  numeric_adjacency <- lapply(adjacency_networks, function(x){x * 1})
  
  # Get edge overlap (defined as "o")
  edge_overlap <- symmetric_matrix_lapply(numeric_adjacency, sum)
  
  # Number of layers (defined as "M")
  layers <- length(adjacency_networks)
  
  # Get copy of edge overlap
  M_edge_overlap <- edge_overlap
  
  # Set non-zero values to 1 (defined as "Heaviside o")
  M_edge_overlap[M_edge_overlap != 0] <- 1
  
  # Return structural overlap
  (layers / (layers - 1)) *
  ((edge_overlap / (layers * M_edge_overlap)) - (1 / layers))

}

# A key question on order of edge list:
# Should edges be in order of their pairwise correspondence; for example:
#
# 1 2
# 3 5
# 4 5
# 5 6
#
# Reads: "1, 2, 3, 5, 4, 5, 5, 6"
#
# OR
#
# Reads (current implementation): "1, 3, 4, 5, 2, 5, 5, 6"
