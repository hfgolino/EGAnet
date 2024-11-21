#' @title Ergodicity Information Index
#'
#' @description Computes the Ergodicity Information Index
#'
#' @param dynEGA.object A \code{\link[EGAnet]{dynEGA.ind.pop}} object
#'
#' @param use Character (length = 1).
#' A string indicating what network element will be used
#' to compute the algorithm complexity, the list of edges or the weights of the network.
#' Defaults to \code{use = "unweighted"}.
#' Current options are:
#'
#' \itemize{
#'
#' \item \code{"edge.list"} --- Calculates the algorithm complexity using the list of edges
#'
#' \item \code{"unweighted"} --- Calculates the algorithm complexity using the binary weights of the encoded prime 
#' transformed network. 0 = edge absent and 1 = edge present
#' 
#' \item \code{"weighted"} --- Calculates the algorithm complexity using the weights of encoded prime-weight transformed network
#' 
#' }
#' 
#' @param shuffles Numeric.
#' Number of shuffles used to compute the Kolmogorov complexity.
#' Defaults to \code{5000}
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
#' eii <- ergoInfo(dyn.ega1)}
#' 
#' @references 
#' \strong{Original Implementation} \cr
#' Golino, H., Nesselroade, J. R., & Christensen, A. P. (2022).
#' Toward a psychology of individuals: The ergodicity information index and a bottom-up approach for finding generalizations.
#' \emph{PsyArXiv}.
#'
#' @export
#'
# Ergodicity Information Index ----
# Updated 12.11.2023
ergoInfo <- function(
    dynEGA.object,
    use = c("edge.list", "unweighted", "weighted"),
    shuffles = 5000
)
{
  
  # Check for missing arguments (argument, default, function)
  use <- set_default(use, "unweighted", ergoInfo)
  
  # Check for appropriate class ("dynEGA.ind.pop" defunct to legacy)
  if(!is(dynEGA.object, "dynEGA") & !is(dynEGA.object, "dynEGA.ind.pop")){
    class_error(dynEGA.object, "dynEGA", "ergoInfo")
  }
  
  # Control exponential notation of numbers
  user_scipen <- options("scipen")$scipen
  
  # Set option
  options(scipen = 0)
  
  # Get proper objects (if not, send an error)
  dynEGA.object <- get_dynEGA_object(dynEGA.object)
  
  # Get individual networks
  individual_networks <- lapply(dynEGA.object$individual, function(x){x$network})
  
  # Get sequence for number of individuals
  individual_sequence <- seq_along(individual_networks)
  
  # Initialize encoding matrix
  dimensions <- dim(individual_networks[[1]])
  encoding_matrix <- matrix(1, nrow = dimensions[1], ncol = dimensions[2])
  
  # "unweighted" and "weighted" needs canonical prime association
  if(use != "edge.list"){
    
    # Set individual networks based on use
    if(use == "unweighted"){
      
      # Binarize networks
      individual_networks <- lapply(
        individual_networks, # NAs occur when there is zero variance 
        function(x){return(x != 0 & !is.na(x))}
      )
      
    }
    
    # Order networks
    individual_networks <- individual_networks[
      order(nvapply(individual_networks, edge_count), decreasing = FALSE)
    ]
    
    # Get prime numbers
    prime_numbers <- get(
      data("prime.num", package = "EGAnet", envir = environment())
    )[individual_sequence]
    
    # Get prime weights
    prime_weights <- lapply(individual_sequence, function(case){
      return(prime_numbers[[case]] ^ individual_networks[[case]])
    })
    
    # Get encoding matrix
    encoding_matrix <- Reduce("*", prime_weights)

    # Revert 1s to 0s
    encoding_matrix[encoding_matrix == 1] <- 0
    
  }
  
  # Get edge list (matches {igraph})
  edge_list <- sparse_network(encoding_matrix)
  
  # Get edge list rows (used at the end)
  edge_rows <- dim(edge_list)[1]
  
  # Get edge list sequence
  edge_sequence <- seq_len(edge_rows)
  
  # Set up population
  
  # Get population adjacency
  population_edges <- dynEGA.object$population$network != 0
  
  # Initialize population encoding matrix
  population_encoding <- matrix(1, nrow = dimensions[1], ncol = dimensions[2])
  
  # Branch based on "use"
  if(use != "edge.list"){
    
    # Check for weighted vs. unweighted
    population_encoding <- 2 ^ swiftelse(
      use == "unweighted", population_edges,
      dynEGA.object$population$network
    )
    
    # Revert 1s to 0s
    population_encoding[population_encoding == 1] <- 0
    
  }
  
  # Set upper triangle to FALSE
  population_edges[upper.tri(population_edges)] <- FALSE
  
  # Get edge list ("col" then "row" matches {igraph})
  population_edge_list <- cbind(
    which(population_edges, arr.ind = TRUE)[,c("col", "row")], 
    population_encoding[population_edges]
  )
  
  # Get edge list sequence
  population_edge_sequence <- nrow_sequence(population_edge_list)
  
  # K-complexity
  
  # Use `keep_weights` for quick indexing
  keep_weights <- swiftelse(
    use == "edge.list",
    c(1L, 2L), c(1L, 2L, 3L)
  )
  
  # Transpose edge lists (more efficient)
  edge_list <- t(edge_list[,keep_weights])
  population_edge_list <- t(population_edge_list[,keep_weights])
  
  # Get k-complexity for individuals and population
  kcomplexities <- lapply(
    seq_len(shuffles), function(iteration){
      
      # Set up return
      return(
        c(
          individual = k_complexity(
            edge_list[,shuffle_replace(edge_sequence)]
          ),
          population = k_complexity(
            population_edge_list[,shuffle_replace(population_edge_sequence)]
          )
        )
      )
      
    }
  )
  
  # Pre-compute values
  mean_individual_complexity <- mean(
    nvapply(kcomplexities, function(x){x["individual"]}), 
    na.rm = TRUE
  )
  mean_population_complexity <- mean(
    nvapply(kcomplexities, function(x){x["population"]}), 
    na.rm = TRUE
  )
  
  # Set up results
  results <- list(
    KComp = mean_individual_complexity,
    KComp.pop = mean_population_complexity,
    EII = sqrt(dynEGA.object$population$n.dim + 1)^(
      (mean_individual_complexity / mean_population_complexity) / log(edge_rows)
    )
  )
  
  # Check for prime weights
  if(use != "edge.list"){
    results$PrimeWeight <- remove_attributes(encoding_matrix)
    results$PrimeWeight.pop <- remove_attributes(population_encoding)
  }
  
  # Add "methods" attribute
  attr(results, "methods") <- list(use = use, shuffles = shuffles)
  
  # Add class
  class(results) <- "EII"
  
  # Restore user's "scipen" option
  options(scipen = user_scipen)
  
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
# r_sample_seeds <- r_sample_seeds
# r_sample_with_replacement <- r_sample_with_replacement
# r_sample_without_replacement <- r_sample_without_replacement
# Need above functions for testing!

#' @exportS3Method 
# S3 Print Method
# Updated 12.11.2023
print.EII <- function(x, ...)
{
  
  # Print EII method
  cat(
    paste0(
      "EII Method: ",
      switch(
        attr(x, "methods")$use,
        "edge.list" = "Edge List",
        "unweighted" = "Unweighted",
        "weighted" = "Weighted"
      ), "\n",
      "Shuffles: ", attr(x, "methods")$shuffles, "\n"
    )
  )
  
  # Print EII value
  cat("EII: ", x$EII)
  
}

#' @exportS3Method 
# S3 Summary Method
# Updated 14.07.2023
summary.EII <- function(object, ...)
{
  print(object, ...) # same as print
}

#' @noRd
# k-complexity ----
# Updated 12.11.2023
k_complexity <- function(values)
{
  
  # Streamlined form
  return(
    length( # length of compression
      memCompress( # bit compression
        paste0( # bits (matches `toString`)
          values, collapse = ", "
        ), type = "gzip" # type of compression
      )
    )
  )
  
}