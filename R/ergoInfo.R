#' Ergodicity Information Index
#' @description Computes the Ergodicity Information Index
#'
#' @param data A dynEGA.ind.pop object
#'
#' @param use Character.
#' A string indicating what network element will be used to compute the algorithm complexity, the list of edges or the weights of the network.
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
#' \item{EII}{The Ergodicity Information Index}
#'
#' @author Hudson Golino <hfg9s at virginia.edu>
#'
#'
#' @export
# Ergodicity Information Index
# Updated 04.27.2022
ergoInfo <- function(data, use = c("edge.list", "weights")){

  #### MISSING ARGUMENTS HANDLING ####
  if(missing(use))
  {use <- "edge.list"
  }else{use}

  # dynEGA (Individual)
  ids <- unique(data$Derivatives$EstimatesDF[,ncol(data$Derivatives$EstimatesDF)])
  ids <- paste0("ID", ids)
  # List of Individual Networks ----
  net.list <- lapply(seq_along(ids), function(i){
    res_list <- list()
    res_list$IDs <- ids[i] # IDs
    res_list$Network <- data$dynEGA.ind[[ids[i]]]$network # Networks
    res_list$igraph.Network <- convert2igraph(res_list$Network) # igraph Network
    res_list$gsize.net <- igraph::gsize(res_list$igraph.Network) # igraph Size
    res_list$adj.net <- as.matrix(igraph::get.adjacency(res_list$igraph.Network, type="both"))
    return(res_list)
  })

  gsize.net.vec <- unlist(
    lapply(net.list, function(x){
      x$gsize.net
    })
  )

  # order by size:
  mat.size <- data.frame(Ind = seq_along(ids), Size = gsize.net.vec)
  mat.size <- mat.size[order(mat.size$Size, decreasing = FALSE),]

  # Load deep learning neural network weights
  prime.num <- get(data("prime.num", envir = environment()))

  # Associate prime number:
  mat.size$Prime <- prime.num[1:nrow(mat.size)]

  # Organize by size:
  net.list.ordered <- net.list[mat.size$Ind]

  # Transforming the matrices into a grid (ncol and nrow = ID 1)
  grid <- expand.grid(x = 1:ncol(net.list.ordered[[1]]$Network), y = 1:ncol(net.list.ordered[[1]]$Network))

  # Computing the Prime-Weight Transformation
  mat.encode.list <- lapply(seq_along(ids), function(i){

    apply(grid, 1,function(x){
      return(
        ifelse(
          net.list.ordered[[i]]$Network[x[1],x[2]]==0,
          1,
          (net.list.ordered[[i]]$adj.net[x[1], x[2]]*mat.size[i,"Prime"])^(net.list.ordered[[i]]$Network[x[1], x[2]])
        )
      )
    })

  })

  encode <- Reduce("*", mat.encode.list)
  mat.encode <- matrix(encode, ncol = ncol(net.list.ordered[[1]]$Network), nrow = nrow(net.list.ordered[[1]]$Network))
  mat.encode <- ifelse(mat.encode==1, 0, mat.encode)
  mat.encode.igraph <- convert2igraph(mat.encode)
  edge.list.mat.encode <- igraph::get.edgelist(mat.encode.igraph)

  if(use == "edge.list"){

    iter <- 1:1000

    bits <- lapply(iter, function(i){
      toString(edge.list.mat.encode[sample(1:nrow(edge.list.mat.encode), size = nrow(edge.list.mat.encode), replace = TRUE)])
    })

    compression <- lapply(iter, function(i){
      memCompress(bits[[i]], "gzip")
    })

    kcomp <- lapply(iter, function(i){
      length(compression[[i]])
    })

  }else{

    edge.list.mat.encode2 <- edge.list.mat.encode

    edge.list.mat.encode2 <- cbind(edge.list.mat.encode2, NA)

    for(i in 1:nrow(edge.list.mat.encode2)){
      edge.list.mat.encode2[i,3] <- mat.encode[edge.list.mat.encode2[i,1],edge.list.mat.encode2[i,2]]
    }

    iter <- 1:1000

    bits <- lapply(iter, function(i){
      toString(edge.list.mat.encode2[sample(1:nrow(edge.list.mat.encode2), size = nrow(edge.list.mat.encode2), replace = TRUE),3])
    })

    compression <- lapply(iter, function(i){
      memCompress(bits[[i]], "gzip")
    })

    kcomp <- lapply(iter, function(i){
      length(compression[[i]])
    })

  }

  # List of Population Networks ----
  net.list.pop <- lapply(1, function(i){
    res_list <- list()
    res_list$IDs <- i # ID
    res_list$Network <- data$dynEGA.pop$network # Network
    res_list$igraph.Network <- convert2igraph(res_list$Network) # igraph Network
    res_list$gsize.net <- igraph::gsize(res_list$igraph.Network) # igraph Size
    res_list$adj.net <- as.matrix(igraph::get.adjacency(res_list$igraph.Network,type="both"))
    return(res_list)
  })

  gsize.net.vec.pop <- net.list.pop[[1]]$gsize.net

  # order by size:
  mat.size.pop <- data.frame(Ind = 1, Size = gsize.net.vec.pop)
  mat.size.pop <- mat.size.pop[order(mat.size.pop$Size, decreasing = FALSE),]

  # Associate prime number:
  mat.size.pop$Prime <- prime.num[1:nrow(mat.size.pop)]

  # Organize by size:
  net.list.ordered.pop <- net.list.pop[mat.size.pop$Ind]

  # Transforming the matrices into a grid (ncol and nrow = ID 1)
  grid.pop <- expand.grid(x = 1:ncol(net.list.ordered.pop[[1]]$Network), y = 1:ncol(net.list.ordered.pop[[1]]$Network))

  # Computing the Prime-Weight Transformation
  mat.encode.list.pop <- lapply(1, function(i){
    apply(grid.pop, 1 , function(x){
      return(
        ifelse(
          net.list.ordered.pop[[i]]$Network[x[1], x[2]]==0,
          1,
          (net.list.ordered.pop[[i]]$adj.net[x[1], x[2]]*mat.size.pop[i,"Prime"])^(net.list.ordered.pop[[i]]$Network[x[1], x[2]])
        )
      )
    })
  })

  encode.pop <- Reduce("*", mat.encode.list.pop)
  mat.encode.pop <- matrix(encode.pop, ncol = ncol(net.list.ordered.pop[[1]]$Network), nrow = nrow(net.list.ordered.pop[[1]]$Network))
  mat.encode.pop <- ifelse(mat.encode.pop==1, 0, mat.encode.pop)
  mat.encode.igraph.pop <- convert2igraph(mat.encode.pop)
  edge.list.mat.encode.pop <- igraph::get.edgelist(mat.encode.igraph.pop)

  if(use == "edge.list"){

    iter <- 1:1000

    bits.pop <- lapply(iter, function(i){
      toString(edge.list.mat.encode.pop[sample(1:nrow(edge.list.mat.encode.pop), size = nrow(edge.list.mat.encode.pop), replace = TRUE)])
    })

    compression.pop <- lapply(iter, function(i){
      memCompress(bits.pop[[i]], "gzip")
    })

    kcomp.pop <- lapply(iter, function(i){
      length(compression.pop[[i]])
    })

  }else{

    edge.list.mat.encode.pop2 <- edge.list.mat.encode.pop

    edge.list.mat.encode.pop2 <- cbind(edge.list.mat.encode.pop2, NA)

    for(i in 1:nrow(edge.list.mat.encode.pop2)){
      edge.list.mat.encode.pop2[i,3] <- mat.encode.pop[edge.list.mat.encode.pop2[i,1],edge.list.mat.encode.pop2[i,2]]
    }

    iter <- 1:1000

    bits.pop <- lapply(iter, function(i){
      toString(edge.list.mat.encode.pop2[sample(1:nrow(edge.list.mat.encode.pop2), size = nrow(edge.list.mat.encode.pop2), replace = TRUE),3])
    })

    compression.pop <- lapply(iter, function(i){
      memCompress(bits.pop[[i]], "gzip")
    })

    kcomp.pop <- lapply(iter, function(i){
      length(compression.pop[[i]])
    })

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
