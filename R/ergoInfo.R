#' Ergodicity Information Index
#' @description Computes the Ergodicity Information Index
#'
#' @param data A dynEGA.ind.pop object
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
#' \item{EII}{The Ergodicity Information Index}
#'
#' @author Hudson Golino <hfg9s at virginia.edu>
#'
#'
#' @export
# Ergodicity Information Index
# Updated 21.10.2020
ergoInfo <- function(data){
  # dynEGA (Individual)
  ids <- unique(dplyr::last(data$Derivatives$EstimatesDF))

  data.ind <- data$dynEGA.ind

  # List of Individual Networks:

  net.list <- vector("list")
  lists <- vector("list")

  for(i in 1:length(ids)){
    net.list[[i]] <- lists
  }

  for(i in 1:length(ids)){
    net.list[[i]]$IDs <- ids[i]
    net.list[[i]]$Network <- data$dynEGA.ind[[ids[i]]]$network
  }

  data.pop <- data$dynEGA.pop
  ids.pop <- 1
  net.list.pop <- vector("list")
  lists.pop <- vector("list")

  for(i in 1:length(ids.pop)){
    net.list.pop[[i]] <- lists.pop
  }


  for(i in 1:1){
    net.list.pop[[i]]$IDs <- ids.pop
    net.list.pop[[i]]$Network <- data$dynEGA.pop$network
  }

  ## Complexity - Individual Networks:
  # Number of edges:

  for(i in 1:length(ids)){
    net.list[[i]]$igraph.Network <- igraph::as.igraph(qgraph::qgraph(net.list[[i]]$Network, DoNotPlot=TRUE))
  }


  for(i in 1:length(ids)){
    net.list[[i]]$gsize.net <- igraph::gsize(net.list[[i]]$igraph.Network)
  }

  gsize.net.vec <- vector()
  for(i in 1:length(ids)){
    gsize.net.vec[i] <- net.list[[i]]$gsize.net
  }

  # order by size:
  mat.size <- data.frame(Ind = 1:length(ids), Size = gsize.net.vec)
  mat.size <- mat.size[order(mat.size$Size, decreasing = FALSE),]

  # Associate prime number:
  prime.num <- primes::generate_primes(min = 0, max = 10^5)
  mat.size$Prime <- prime.num[1:nrow(mat.size)]

  # Get the adjacency matrix:
  for(i in 1:length(ids)){
    net.list[[i]]$adj.net <- as.matrix(igraph::get.adjacency(net.list[[i]]$igraph.Network,type="both"))
  }

  # Organize by size:

  net.list.ordered <- net.list[mat.size$Ind]

  # Transforming the matrices into a grid (ncol and nrow = ID 1)

  grid <- expand.grid(x = 1:ncol(net.list.ordered[[1]]$Network), y = 1:ncol(net.list.ordered[[1]]$Network))


  # Computing the Prime-Weight Transformation

  mat.encode.list <- vector("list")
  for(i in 1:length(ids)){
    mat.encode.list[[i]] <- apply(grid, 1,
                                  function(x){return(
                                    ifelse(net.list.ordered[[i]]$Network[x[1], x[2]]==0,1,(net.list.ordered[[i]]$adj.net[x[1], x[2]]*mat.size[i,"Prime"])^(net.list.ordered[[i]]$Network[x[1], x[2]]))
                                  )})
  }

  encode <- Reduce("*", mat.encode.list)
  mat.encode <- matrix(encode, ncol = ncol(net.list.ordered[[1]]$Network), nrow = nrow(net.list.ordered[[1]]$Network))
  mat.encode <- ifelse(mat.encode==1, 0, mat.encode)
  mat.encode.igraph <- igraph::as.igraph(qgraph::qgraph(mat.encode, layout = "spring", DoNotPlot = TRUE))
  edge.list.mat.encode <- igraph::get.edgelist(mat.encode.igraph)


  bits <- vector("list")
  compression <- vector("list")
  kcomp <- vector("list")


  for(i in 1:1000){
    bits[[i]] <- toString(edge.list.mat.encode[sample(1:nrow(edge.list.mat.encode), size = nrow(edge.list.mat.encode), replace = TRUE)])
    compression[[i]] <- memCompress(bits[[i]], "gzip")
    kcomp[[i]] <- length(compression[[i]])
  }


  ## Complexity - Population Network:
  # Number of edges:

  for(i in 1:length(ids.pop)){
    net.list.pop[[i]]$igraph.Network <- igraph::as.igraph(qgraph::qgraph(net.list.pop[[i]]$Network, DoNotPlot=TRUE))
  }


  for(i in 1:length(ids.pop)){
    net.list.pop[[i]]$gsize.net <- igraph::gsize(net.list.pop[[i]]$igraph.Network)
  }

  gsize.net.vec.pop <- vector()
  for(i in 1:length(ids.pop)){
    gsize.net.vec.pop[i] <- net.list.pop[[i]]$gsize.net
  }

  # order by size:
  mat.size.pop <- data.frame(Ind = 1:length(ids.pop), Size = gsize.net.vec.pop)
  mat.size.pop <- mat.size.pop[order(mat.size.pop$Size, decreasing = FALSE),]

  # Associate prime number:
  prime.num.pop <- primes::generate_primes(min = 0, max = 10^5)
  mat.size.pop$Prime <- prime.num.pop[1:nrow(mat.size.pop)]

  # Get the adjacency matrix:
  for(i in 1:length(ids.pop)){
    net.list.pop[[i]]$adj.net <- as.matrix(igraph::get.adjacency(net.list.pop[[i]]$igraph.Network,type="both"))
  }

  # Organize by size:

  net.list.ordered.pop <- net.list.pop[mat.size.pop$Ind]

  # Transforming the matrices into a grid (ncol and nrow = ID 1)

  grid.pop <- expand.grid(x = 1:ncol(net.list.ordered.pop[[1]]$Network), y = 1:ncol(net.list.ordered.pop[[1]]$Network))


  # Computing the Prime-Weight Transformation

  mat.encode.list.pop <- vector("list")
  for(i in 1:length(ids.pop)){
    mat.encode.list.pop[[i]] <- apply(grid.pop, 1,
                                      function(x){return(
                                        ifelse(net.list.ordered.pop[[i]]$Network[x[1], x[2]]==0,1,(net.list.ordered.pop[[i]]$adj.net[x[1], x[2]]*mat.size.pop[i,"Prime"])^(net.list.ordered.pop[[i]]$Network[x[1], x[2]]))
                                      )})
  }

  encode.pop <- Reduce("*", mat.encode.list.pop)
  mat.encode.pop <- matrix(encode.pop, ncol = ncol(net.list.ordered.pop[[1]]$Network), nrow = nrow(net.list.ordered.pop[[1]]$Network))
  mat.encode.pop <- ifelse(mat.encode.pop==1, 0, mat.encode.pop)
  mat.encode.igraph.pop <- igraph::as.igraph(qgraph::qgraph(mat.encode.pop, layout = "spring", DoNotPlot = TRUE))
  edge.list.mat.encode.pop <- igraph::get.edgelist(mat.encode.igraph.pop)


  bits.pop <- vector("list")
  compression.pop <- vector("list")
  kcomp.pop <- vector("list")


  for(i in 1:1000){
    bits.pop[[i]] <- toString(edge.list.mat.encode.pop[sample(1:nrow(edge.list.mat.encode.pop), size = nrow(edge.list.mat.encode.pop), replace = TRUE)])
    compression.pop[[i]] <- memCompress(bits.pop[[i]], "gzip")
    kcomp.pop[[i]] <- length(compression.pop[[i]])
  }


  # Kolmogorov Complexity:
  results <- list()
  results$PrimeWeight <- mat.encode
  results$PrimeWeight.pop <- mat.encode.pop
  results$KComp <- mean(unlist(kcomp))
  results$KComp.pop <- mean(unlist(kcomp.pop))

  ergo.info.index<- sqrt(data$dynEGA.pop$n.dim)^((results$KComp.pop/results$KComp)/log(sum(!results$PrimeWeight.pop==0)))
  results$EII <- ergo.info.index
  return(results)
}
#----
