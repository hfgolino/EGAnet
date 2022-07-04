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
#' \item{\strong{\code{edge.list}}}
#' {Calculates the algorithm complexity using the list of edges.}
#'
#' \item{\strong{\code{weights}}}
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
#' @export
#' 
# Ergodicity Information Index
# Updated 04.07.2022
ergoInfo <- function(
    dynEGA.object,
    use = c(
      "edge.list",
      "unweighted",
      "weighted"
    )
)
{
  
  #### MISSING ARGUMENTS HANDLING ####
  if(missing(use)){use <- "edge.list"}
  
  # Check for class
  if(!is(dynEGA.object, "dynEGA.ind.pop")){
    stop(
      paste(
        "Input into the `dynEGA.object` argument's class is not `dynEGA.ind.pop`.\n\n",
        "Class of dynEGA.object = ", paste(
          class(dynEGA.object), sep = "", collapse = ", "
        ),
        sep = ""
      )
    )
  }
  
  # Sort population- and individual-level outputs
  dynEGA.pop <- dynEGA.object$dynEGA.pop
  dynEGA.ind <- dynEGA.object$dynEGA.ind
  
  # Remove Methods
  if("Methods" %in% names(dynEGA.ind$dynEGA)){
    dynEGA.ind$dynEGA <- dynEGA.ind$dynEGA[-which(names(dynEGA.ind$dynEGA) == "Methods")]
  }
  
  # Edge list
  if(use == "edge.list"){
    
    # Obtain individual networks
    individual_networks <- lapply(dynEGA.ind$dynEGA, function(x){x$network})
    
    # Stack networks in an array
    arrayed_networks <- simplify2array(individual_networks)
    
    # Obtain edges from networks
    edges <- apply(arrayed_networks, 1:2, sum) != 0
    
    # Make upper triangle equal FALSE
    edges[upper.tri(edges)] <- FALSE
    
    # Obtain edge list
    edge_list <- which(edges, arr.ind = TRUE)
    
    # Set number of edges (10^3 based on Santora & Nicosia, 2020)
    iter <- 1:1000
    
    # Obtain bits
    individual_bits <- lapply(iter, function(i){
      toString(edge_list[
        sample(
          1:nrow(edge_list),
          size = nrow(edge_list), 
          replace = TRUE
        ),
      ])
    })
    
    # Compress bits
    individual_compression <- lapply(iter, function(i){
      memCompress(individual_bits[[i]], "gzip")
    })
    
    # Obtain complexity
    individual_kcomp <- lapply(iter, function(i){
      length(individual_compression[[i]])
    })
    
    # Obtain population network
    population_network <- dynEGA.pop$dynEGA$network
    
    # Obtain edges from networks
    population_edges <- population_network != 0
    
    # Make upper triangle equal FALSE
    population_edges[upper.tri(population_edges)] <- FALSE
    
    # Obtain edge list
    population_edge_list <- which(population_edges, arr.ind = TRUE)
    
    # Obtain bits
    population_bits <- lapply(iter, function(i){
      toString(population_edge_list[
        sample(
          1:nrow(population_edge_list),
          size = nrow(population_edge_list), 
          replace = TRUE
        ),
      ])
    })
    
    # Compress bits
    population_compression <- lapply(iter, function(i){
      memCompress(population_bits[[i]], "gzip")
    })
    
    # Obtain complexity
    population_kcomp <- lapply(iter, function(i){
      length(population_compression[[i]])
    })
    
    # Kolmogorov Complexity:
    results <- list()
    results$KComp <- mean(unlist(individual_kcomp))
    results$KComp.pop <- mean(unlist(population_kcomp))
    results$complexity <- results$KComp / results$KComp.pop
    results$EII  <- sqrt(dynEGA.pop$dynEGA$n.dim)^((results$KComp/results$KComp.pop)/log(nrow(population_edge_list)))
    results$use <- use
    class(results) <- "EII"
    
  }else if(use == "unweighted"){
    
    # Obtain individual networks
    individual_networks <- lapply(dynEGA.ind$dynEGA, function(x){x$network})
    
    # Obtain adjacency matrices
    individual_adjacency <- lapply(individual_networks, function(x){
      ifelse(x != 0, 1, 0)
    })
    
    # Obtain number of edges
    edge_count <- unlist(lapply(individual_adjacency, function(x){sum(x) / 2}))
    
    # Order networks and adjacency matrices by size
    ordering <- order(edge_count, decreasing = FALSE)
    individual_networks <- individual_networks[ordering]
    individual_adjacency <- individual_adjacency[ordering]
    
    # Obtain prime numbers
    prime.num <- get(data(
      "prime.num",
      package = "EGAnet",
      envir = environment()
    ))
    
    # Get prime numbers equal to networks
    prime_numbers <- prime.num[1:length(individual_networks)]
    
    # Obtain prime weights
    prime_weights <- lapply(seq_along(individual_networks), function(i){
      
      # Assign prime
      prime_network <- individual_adjacency[[i]] * prime_numbers[i]
      
      # Assign 1s to 0s
      prime_network <- ifelse(prime_network == 0, 1, prime_network)
      
      # Raise to weight power
      power_prime_network <- prime_network
      
    })
    
    # Reduce to encoding matrix
    encoding_matrix <- Reduce("*", prime_weights)
    
    # Make 1 encodings 0
    encoding_matrix <- ifelse(encoding_matrix == 1, 0, encoding_matrix)
    
    # Obtain edges from encoding matrix
    edges <- encoding_matrix != 0
    
    # Make upper triangle equal FALSE
    edges[upper.tri(edges)] <- FALSE
    
    # Obtain edge list
    edge_list <- as.data.frame(which(edges, arr.ind = TRUE))
    edge_list$weight <- encoding_matrix[edges]
    
    # Set number of edges (10^3 based on Santora & Nicosia, 2020)
    iter <- 1:1000
    
    # Obtain bits
    individual_bits <- lapply(iter, function(i){
      toString(edge_list$weight[
        sample(
          1:nrow(edge_list),
          size = nrow(edge_list), 
          replace = TRUE
        )
      ])
    })
    
    # Compress bits
    individual_compression <- lapply(iter, function(i){
      memCompress(individual_bits[[i]], "gzip")
    })
    
    # Obtain complexity
    individual_kcomp <- lapply(iter, function(i){
      length(individual_compression[[i]])
    })
    
    # Obtain population network
    population_network <- dynEGA.pop$dynEGA$network
    
    # Obtain population adjacency
    population_adjacency <- ifelse(population_network != 0, 1, 0)
    
    # Obtain prime matrix
    population_prime <- population_adjacency * 2 # smallest prime
    
    # Obtain edges from networks
    population_edges <- population_network != 0
    
    # Make upper triangle equal FALSE
    population_edges[upper.tri(population_edges)] <- FALSE
    
    # Obtain edge list
    population_edge_list <- as.data.frame(which(population_edges, arr.ind = TRUE))
    population_edge_list$weight <- population_prime[population_edges]
    
    # Set number of edges (10^3 based on Santora & Nicosia, 2020)
    iter <- 1:1000
    
    # Obtain bits
    population_bits <- lapply(iter, function(i){
      toString(population_edge_list$weight[
        sample(
          1:nrow(population_edge_list),
          size = nrow(population_edge_list), 
          replace = TRUE
        )
      ])
    })
    
    # Compress bits
    population_compression <- lapply(iter, function(i){
      memCompress(population_bits[[i]], "gzip")
    })
    
    # Obtain complexity
    population_kcomp <- lapply(iter, function(i){
      length(population_compression[[i]])
    })
    
    # Kolmogorov Complexity:
    results <- list()
    results$KComp <- mean(unlist(individual_kcomp))
    results$KComp.pop <- mean(unlist(population_kcomp))
    results$complexity <- results$KComp / results$KComp.pop
    results$EII  <- sqrt(dynEGA.pop$dynEGA$n.dim)^((results$KComp/results$KComp.pop)/log(nrow(population_edge_list)))
    results$use <- use
    class(results) <- "EII"
    
    
  }else if(use == "weighted"){
    
    # Obtain individual networks
    individual_networks <- lapply(dynEGA.ind$dynEGA, function(x){x$network})
    
    # Obtain adjacency matrices
    individual_adjacency <- lapply(individual_networks, function(x){
      ifelse(x != 0, 1, 0)
    })
    
    # Obtain number of edges
    edge_count <- unlist(lapply(individual_adjacency, function(x){sum(x) / 2}))
    
    # Order networks and adjacency matrices by size
    ordering <- order(edge_count, decreasing = FALSE)
    individual_networks <- individual_networks[ordering]
    individual_adjacency <- individual_adjacency[ordering]
    
    # Obtain prime numbers
    prime.num <- get(data(
      "prime.num",
      package = "EGAnet",
      envir = environment()
    ))
    
    # Get prime numbers equal to networks
    prime_numbers <- prime.num[1:length(individual_networks)]
    
    # Obtain prime weights
    prime_weights <- lapply(seq_along(individual_networks), function(i){
      
      # Assign prime
      prime_network <- individual_adjacency[[i]] * prime_numbers[i]
      
      # Assign 1s to 0s
      prime_network <- ifelse(prime_network == 0, 1, prime_network)
      
      # Raise to weight power
      power_prime_network <- prime_network^individual_networks[[i]]
      
    })
    
    # Reduce to encoding matrix
    encoding_matrix <- Reduce("*", prime_weights)
    
    # Make 1 encodings 0
    encoding_matrix <- ifelse(encoding_matrix == 1, 0, encoding_matrix)
    
    # Obtain edges from encoding matrix
    edges <- encoding_matrix != 0
    
    # Make upper triangle equal FALSE
    edges[upper.tri(edges)] <- FALSE
    
    # Obtain edge list
    edge_list <- as.data.frame(which(edges, arr.ind = TRUE))
    edge_list$weight <- encoding_matrix[edges]
    
    # Set number of edges (10^3 based on Santora & Nicosia, 2020)
    iter <- 1:1000
    
    # Obtain bits
    individual_bits <- lapply(iter, function(i){
      toString(edge_list$weight[
        sample(
          1:nrow(edge_list),
          size = nrow(edge_list), 
          replace = TRUE
        )
      ])
    })
    
    # Compress bits
    individual_compression <- lapply(iter, function(i){
      memCompress(individual_bits[[i]], "gzip")
    })
    
    # Obtain complexity
    individual_kcomp <- lapply(iter, function(i){
      length(individual_compression[[i]])
    })
    
    # Obtain population network
    population_network <- dynEGA.pop$dynEGA$network
    
    # Obtain population adjacency
    population_adjacency <- ifelse(population_network != 0, 1, 0)
    
    # Obtain prime matrix
    population_prime <- (population_adjacency * 2)^population_network
    
    # Obtain edges from networks
    population_edges <- population_network != 0
    
    # Make upper triangle equal FALSE
    population_edges[upper.tri(population_edges)] <- FALSE
    
    # Obtain edge list
    population_edge_list <- as.data.frame(which(population_edges, arr.ind = TRUE))
    population_edge_list$weight <- population_prime[population_edges]
    
    # Set number of edges (10^3 based on Santora & Nicosia, 2020)
    iter <- 1:1000
    
    # Obtain bits
    population_bits <- lapply(iter, function(i){
      toString(population_edge_list$weight[
        sample(
          1:nrow(population_edge_list),
          size = nrow(population_edge_list), 
          replace = TRUE
        )
      ])
    })
    
    # Compress bits
    population_compression <- lapply(iter, function(i){
      memCompress(population_bits[[i]], "gzip")
    })
    
    # Obtain complexity
    population_kcomp <- lapply(iter, function(i){
      length(population_compression[[i]])
    })
    
    # Kolmogorov Complexity:
    results <- list()
    results$KComp <- mean(unlist(individual_kcomp))
    results$KComp.pop <- mean(unlist(population_kcomp))
    results$complexity <- results$KComp / results$KComp.pop
    results$EII  <- sqrt(dynEGA.pop$dynEGA$n.dim)^((results$KComp/results$KComp.pop)/log(nrow(population_edge_list)))
    results$use <- use
    class(results) <- "EII"
    
  }
  
  return(results)
}
#----
