#%%%%%%%%%%%%%%%%%
# DEVELOPMENT ----
#%%%%%%%%%%%%%%%%%

#' @noRd
#'
# Polytomous IRT parameters
# Updated 20.03.2022
poly.irt <- function(loadings, data)
{
  # Check for standardized loadings
  if(any(class(loadings) == "NetLoads")){
    loadings <- loadings$std
  }

  # Convert loadings to matrix
  loadings <- as.matrix(loadings)

  # Unique variance
  s <- sqrt(1 - rowSums(loadings^2))

  # Estimated discrimination parameters
  est_a <- sweep(loadings, 1, s, "/")

  # Estimated threshold parameters
  thresholds <- lavaan::lavCor(
    as.data.frame(data), # must be data frame
    ordered = colnames(data),
    output = "thresholds"
  )

  # Separate thresholds
  threshs <- vector("list", length = ncol(data))
  names(threshs) <- colnames(data)

  # Get thresholds for data
  for(i in 1:ncol(data)){

    # Number of categories
    cats <- length(unique(na.omit(data[,i]))) - 1

    # Initialize vector
    thresh <- numeric(cats)

    # Loop through thresholds
    for(j in 1:cats){
      thresh[j] <- thresholds[
        grep(
          paste(colnames(data)[i], "|t", j, sep = ""), # set exact match
          names(thresholds),
          fixed = TRUE
        )
      ]
    }

    # Insert into thresholds
    threshs[[i]] <- thresh

  }

  # Estimate location parameters
  ## All variables have the same number of categories
  if(all(unlist(lapply(threshs, length)))){

    threshs <- t(simplify2array(threshs))

    est_d <- sweep(-threshs, 1, s, "/")
  }else{## Some variables have different numbers of categories

    est_d <- lapply(threshs, function(x, s){
      -x/s
    }, s = s)

  }

  # Results list
  results <- list()
  results$discrimination <- est_a * 1.702
  results$location <- est_d

  return(results)

}

#%%%%%%%%%%%%%
# hierEGA ----
#%%%%%%%%%%%%%

#' Oblimin rotation (to be consistent in \code{\link[EGAnet]{net.loads}})
#' @noRd
# Updated 28.11.2022 -- Marcos
oblimin_rotate <- function(
    loadings,
    n.rotations = 10, # defaults to 10
    nfactors, # number of factors
    maxit = 1e4, # iterations
    eps = 1e-5, # convergence
    gamma = 0, # gamma in oblimin
    rotate = "oblimin" # rotation
)
{
  
  # Check for number of factors
  if(missing(nfactors)){
    nfactors <- ncol(loadings)
  }
  
  # Check for unidimensional structure (28.11.2022)
  if(nfactors == 1){
    
    # Return results
    return(
      list(
        loadings = loadings,
        Phi = 1
      )
    )
    
  }
  
  x <- list()
  fs <- vector(length = n.rotations)
  
  for(i in 1:n.rotations) {
    
    X <- replicate(nfactors, rnorm(nfactors))
    Q <- qr.Q(qr(X)) # Random orthogonal matrix (initial value)
    x[[i]] <- GPArotation::GPFoblq(A = loadings, method = rotate, 
                                   Tmat = Q, maxit = maxit, eps = eps) # Fit
    fs[i] <- oblimin(x[[i]]$loadings, gamma) # Fit values
    
  }

  # Results
  index <- which.min(fs) # Index pertaining to the minimum fit value
  loadings <- x[[i]]$loadings # Select the corresponding loadings
  Phi <- x[[i]]$Phi # Factor correlation matrix
  
  # Return results
  return(
    list(
      loadings = loadings,
      Phi = Phi
    )
  )
  
}

#' Oblimin rotation
#' @noRd
# Updated 21.11.2022 -- Marcos
oblimin <- function(L, gamma = 0) {
  
  # Fit value for oblimin
  # L = rotated loading matrix
  # gamma = fixed parameter
  
  nfactors <- ncol(L)
  p <- nrow(L)
  I <- diag(p)
  gC <- matrix(gamma/p, p, p)
  IgC <- I - gC
  
  N <- matrix(1, nfactors, nfactors)
  diag(N) <- 0
  L2 <- L*L
  f <- sum(diag(t(L2) %*% (IgC %*% L2 %*% N)))/4
  return(f)
  
}

#' Factor scores
#' @noRd
# Updated 28.11.2022 -- Marcos
fscores <- function(S, loadings, Phi, Shat, scores = NULL, method = "Thurstone") {
  
  # S = empirical correlation matrix
  # loadings = rotated loading matrix
  # Phi = estimated factor correlation matrix
  # Shat = model correlation matrix
  # scores = raw scores
  # method = factor.scores method
  
  if(is.null(scores)) stop("Please, provide the matrix of observed scores")
  
  uniquenesses <- 1 - diag(Shat)
  invS <- solve(S)
  n <- nrow(scores)
  p <- nrow(loadings)
  q <- ncol(loadings)
  z <- scale(scores)
  
  Lambda <- loadings
  Phi <- Phi
  LP <- Lambda %*% Phi # Correlations between factors and items
  
  # Find the weights:
  
  if(method == "regression" | method == "Thurstone") {
    
    weights <- solve(S, LP)
    
  } else if(method == "tenBerge") {
    
    SVD <- svd(Phi)
    Phi12 <- SVD$u %*% diag(sqrt(SVD$d)) %*% t(SVD$v)
    SVD <- svd(S)
    R12 <- SVD$u %*% diag(1/sqrt(SVD$d)) %*% t(SVD$v)
    L <- Lambda %*% Phi12
    SVD <- svd(t(L) %*% invS %*% L)
    LRL12 <- SVD$u %*% diag(1/sqrt(SVD$d)) %*% t(SVD$v)
    C <- R12 %*% L %*% LRL12
    weights <- R12 %*% C %*% Phi12
    
  } else if(method == "Bartlett") {
    
    U <- c(fit$efa$uniquenesses)
    U2 <- diag(1/(U*U))
    weights <- U2 %*% Lambda %*% solve(t(Lambda) %*% U2 %*% Lambda)
    
  } else if(method == "Harman") {
    
    weights <- solve(fit$efa$Rhat) %*% Lambda
    
  }
  
  fs <- z %*% weights # Factor scores
  
  # Validity coefficients:
  # invL <- diag(1/apply(fs, MARGIN = 2, FUN = sd))
  C <- t(weights) %*% S %*% weights
  # invL was updated (28.11.2022) to ensure matrix
  # when nfactors = 1
  invL <- diag(sqrt(diag(C)), ncol(loadings)) # Standard deviations of the factor scores
  validity_univocality <- t(LP) %*% weights %*% invL
  validity <- matrix(diag(validity_univocality), nrow = 1)
  rownames(validity) <- ""
  univocality <- validity_univocality
  diag(univocality) <- NA
  
  # Accuracy:
  accuracy <- stats::cor(fs)
  
  # Standard errors for factor scores:
  r <- matrix(diag(invS), ncol = 1)
  Rj <- matrix(1-c(validity^2), nrow = 1)
  se <- sqrt(r %*% Rj / (n-p-1))
  se <- matrix(se, nrow = p, ncol = q)
  
  colnames(fs) <- colnames(weights) <- colnames(validity) <-
    colnames(univocality) <- rownames(univocality) <-
    colnames(accuracy) <- rownames(accuracy) <- colnames(se) <-
    paste("F", sprintf(paste("%0", nchar(q), "d", sep = ""), 1:q), sep = "")
  
  result <- list(fscores = fs, weights = weights, validity = validity,
                 univocality = univocality, accuracy = accuracy, se = se)
  
  return(result)
  
}

#' Exploratory factor analysis
#' @noRd
# Updated 21.11.2022 -- Marcos
efa <- function(data, nfactors, fm = "minres", rotate = "oblimin",
                n.rotations = 10, maxit = 1e4, factor.scores = "Thurstone",
                gamma = 0, eps = 1e-5) {
  
  # Function to perform EFA
  # data = raw data
  # nfactors = number of factors to extract
  # fm = factor extraction method
  # rotate = rotation method (only oblimin is implemented right now)
  # n.rotations = number of rotations to perform with different initial values
  # maxit = maximum number of iterations for the rotation to converge
  # factor.scores = factor scores' method
  # gamma = fixed parameter for oblimin
  # eps = stop criteria for the rotation (relative tolerance value)
  
  if(rotate != "oblimin") stop("oblimin is the only rotation currently implemented.")
  
  # Compute the correlation matrix
  S <- qgraph::cor_auto(data, verbose = FALSE)
  
  # Factor extraction:
  fa <- psych::fa(S, nfactors = nfactors, fm = fm, rotate = "none")
  # x <- list()
  # fs <- vector(length = n.rotations)
  # 
  # for(i in 1:n.rotations) {
  #   
  #   X <- replicate(nfactors, rnorm(nfactors))
  #   Q <- qr.Q(qr(X)) # Random orthogonal matrix (initial value)
  #   x[[i]] <- GPArotation::GPFoblq(A = fa$loadings, method = rotate, 
  #                                  Tmat = Q, maxit = maxit, eps = eps) # Fit
  #   fs[i] <- oblimin(x[[i]]$loadings, gamma) # Fit values
  #   
  # }
  # 
  # index <- which.min(fs) # Index pertaining to the minimum fit value
  # loadings <- x[[i]]$loadings # Select the corresponding loadings
  # Phi <- x[[i]]$Phi # Factor correlation matrix
  
  # Perform oblimin rotation
  rotation <- oblimin_rotate(
    loadings = fa$loadings,
    n.rotations = n.rotations,
    nfactors = nfactors,
    maxit = maxit,
    eps = eps,
    gamma = gamma,
    rotate = rotate
  )
  loadings <- rotation$loadings # Select the corresponding loadings
  Phi <- rotation$Phi # Factor correlation matrix
  Shat <- fa$model # Model correlation matrix
  uniquenesses <- fa$uniquenesses
  
  # Compute the factor scores
  fs <- fscores(S = S, loadings = loadings, Phi = Phi, 
                Shat = Shat, scores = data, method = factor.scores)
  
  result <- list(loadings = loadings, Phi = Phi, Shat = Shat,
                 uniquenesses = uniquenesses, factor.scores = fs$fscores)
  
}

#%%%%%%%%%%%%%
# LOUVAIN ----
#%%%%%%%%%%%%%

#' @noRd
# Re-index membership function
# Updated 04.05.2023
reindex_comm <- function(wc) {
  
  # Initialize index array
  index_array <- numeric(length(wc))
  
  # Initialize count
  count <- 0
  
  # Loop over to determine index correspondence
  for(i in seq_along(wc)){
    
    # Check if index already exists
    if(index_array[wc[i]] == 0){
      
      # Increase count
      count <- count + 1
      
      # Set index array
      index_array[wc[i]] <- count
      
    }
    
  }
  
  # Update membership
  for(i in seq_along(wc)){
    
    # Check value in index array
    wc[i] <- index_array[wc[i]]
    
  }
  
  # Return re-indexed membership
  return(wc)
  
}

# Lancichinetti & Fortunato (2012)
#' @noRd
# Consensus Clustering
# Updated 17.11.2022
consensus_clustering <- function(
    network, corr,
    order = c("lower", "higher"),
    consensus.iter,
    resolution = 1,
    type = c(
      "all", "highest_modularity",
      "most_common", "iterative",
      "lowest_tefi"
    )
)
{
  # Check for type
  if(missing(type)){
    type <- "all"
  }else{
    type <- match.arg(type, several.ok = TRUE)
  }
  
  # Obtain network names
  network_names <- colnames(network)
  
  # Check for empty network
  if(sum(network) == 0){
    
    # Return individual communities
    wc <- 1:ncol(network)
    
    # Assign names
    names(wc) <- network_names
    
    # Set up results
    results <- list()
    results$highest_modularity <- wc
    results$most_common <- wc
    results$iterative <- wc
    results$lowest_tefi <- wc
    results$summary_table <- "Empty network. No general factors found."
    
    # Return consensus
    return(results)
    
    
  }
  
  # NEED TO FIX DISCONNECTED NODES!!!!

  # Apply Louvain
  communities <- lapply(1:consensus.iter, function(j){

    # igraph output
    output <- signed.louvain(network)
    
    # Obtain order
    if(order == "lower"){
      wc <- output$memberships[1,]
    }else if(order == "higher"){
      wc <- output$memberships[nrow(output$memberships),]
    }

    # Return
    return(wc)

  })

  # Simplify to a matrix
  wc_matrix <- t(simplify2array(communities, higher = FALSE))

  # Make data frame
  df <- as.data.frame(wc_matrix)

  # Obtain duplicate indices
  dupe_ind <- duplicated(df)

  # Rows for non-duplicates
  non_dupes <- data.frame(df[!dupe_ind,])

  # Rows for duplicates
  dupes <- data.frame(df[dupe_ind,])

  # Match duplicates with non-duplicates
  dupe_count <- table(
    match(
      data.frame(t(dupes)), data.frame(t(non_dupes))
    ))

  # Set up summary table
  summary_table <- data.frame(
    N_Dimensions = apply(non_dupes, 1, function(x){
      length(na.omit(unique(x)))
    }),
    Proportion = as.matrix(count(wc_matrix) / nrow(wc_matrix))
  )
  
  # Change column names of non_dupes
  if(!is.null(colnames(network))){
    colnames(non_dupes) <- colnames(network)
  }
  
  # Attach non-duplicate solutions
  summary_table <- cbind(summary_table, non_dupes)
  
  if(type == "all" | type == "most_common" | type == "iterative"){
    
    # Obtain max proportion
    wc_proportion <- unlist(summary_table[
      which.max(summary_table[,"Proportion"]),
      -c(1,2)
    ])
    
  }
  
  if(type == "all" | type == "lowest_tefi"){
    
    # Compute TEFI
    TEFI <- apply(non_dupes, 1, function(x){
      tefi(abs(corr), x)$VN.Entropy.Fit
    })
    
    # Add TEFI to summary table
    summary_table$TEFI <- TEFI
    
    # Obtain minimum TEFI
    wc_tefi <- unlist(summary_table[
      which.min(summary_table[,"TEFI"]),
      -c(1:2)
    ])
    
    # Remove TEFI
    wc_tefi <- wc_tefi[-length(wc_tefi)]
    
  }
  
  if(type == "all" | type == "highest_modularity"){
    
    # Compute modularity
    modularities <- apply(non_dupes, 1, modularity, network, 1)
    
    # Add modularities to summary table
    summary_table$Modularity <- modularities
    
    # Ensure descending order
    summary_table <- summary_table[order(summary_table[,"Modularity"], decreasing = TRUE),]
    row.names(summary_table) <- NULL
    
    # Obtain max modularity
    wc_modularity <- unlist(summary_table[
      which.max(summary_table[,"Modularity"]),
      -c(1:2)
    ])
    
    # Check for TEFI
    if("TEFI" %in% names(wc_modularity)){
      wc_modularity <- wc_modularity[-c(
        length(wc_modularity) - 1,
        length(wc_modularity)
      )]
    }else{
      wc_modularity <- wc_modularity[-c(
        length(wc_modularity)
      )]
    }
    
  }

  # Remove row names
  row.names(summary_table) <- NULL

  # Traditional consensus clustering
  
  if(type == "all" | type == "iterative"){
    
    # Binary check function
    binary <- function(b_matrix){
      all(b_matrix == 0 | b_matrix == 1)
    }
    
    # Initialize count for homogenizing membership
    iter <- 1
    
    # Set up while loop
    while(!binary(network)){
      
      if(iter != 1){
        
        # Apply Louvain
        communities <- lapply(1:consensus.iter, function(j){
          
          # igraph output
          output <- signed.louvain(network)
          
          # Obtain order
          if(order == "lower"){
            wc <- output$memberships[1,]
          }else if(order == "higher"){
            wc <- output$memberships[nrow(output$memberships),]
          }
          
          # Return
          return(wc)
          
        })
        
        # Simplify to a matrix
        wc_matrix <- t(simplify2array(communities, higher = FALSE))
        
      }else{
        
        # Check for non-unique memberships
        if(length(wc_proportion) != length(na.omit(unique(wc_proportion)))){
          
          # Homogenize memberships
          wc_matrix <- t(homogenize.membership(
            target.wc = wc_proportion,
            convert.wc = t(wc_matrix)
          ))
          
        }
        
        # ^^^ checks for whether all variables are in individual
        # communities
        #
        # current workaround for higher order dimensions with
        # singleton dimensions
        
      }
      
      # Get indices for matrix
      d_matrix <- matrix(0, nrow = ncol(wc_matrix), ncol = ncol(wc_matrix))
      
      # Obtain combinations for lower
      combinations <- cbind(
        rep(1:ncol(wc_matrix), times = ncol(wc_matrix)),
        rep(1:ncol(wc_matrix), each = ncol(wc_matrix))
      )
      
      # Fill lower order matrix
      for(i in 1:nrow(combinations)){
        
        # Get indices
        index1 <- combinations[i,1]
        index2 <- combinations[i,2]
        
        d_matrix[index1, index2] <- mean(wc_matrix[,index1] == wc_matrix[,index2], na.rm = TRUE)
        
      }
      
      # Set values less than threshold to zero
      d_matrix <- ifelse(d_matrix <= 0.30, 0, d_matrix)
      
      # Start over
      network <- d_matrix
      
      # Increase count
      iter <- iter + 1
      
    }
    
    # Check for same dimensions
    if(sum(network) == ncol(network)){
      
      # Set to all unique
      wc <- 1:ncol(network)
      
    }else{
      
      # Obtain memberships
      wc <- signed.louvain(network)$memberships
      
      # Obtain order
      if(order == "lower"){
        wc <- wc[1,]
      }else if(order == "higher"){
        wc <- wc[nrow(wc),]
      }
      
      # Ensure vector
      wc <- as.vector(wc)
      
    }
    
    # Assign names
    names(wc) <- network_names
    
    # Assign to traditional
    wc_traditional <- wc
    
  }

  # Set up results
  results <- list()
  
  # Set up results
  if(exists("wc_modularity", envir = environment())){
    results$highest_modularity <- wc_modularity
  }
  
  if(exists("wc_proportion", envir = environment())){
    results$most_common <- wc_proportion
  }
  
  if(exists("wc_traditional", envir = environment())){
    results$iterative <- wc_traditional
  }
  
  if(exists("wc_tefi", envir = environment())){
    results$lowest_tefi <- wc_tefi
  }
  
  results$summary_table <- summary_table
  
  # Return consensus
  return(results)
}

# Lancichinetti & Fortunato (2012)
#' @noRd
# Most Common Consensus Clustering
# Updated 22.07.2022
most_common_consensus <- function(
    network,
    order = c("lower", "higher"),
    consensus.iter,
    resolution = 1
)
{
  
  # Check for empty network
  if(sum(network) == 0){
    
    # Return individual communities
    wc <- 1:ncol(network)
    
    # Assign names
    names(wc) <- colnames(network)
    
    # Set up results
    results <- list()
    results$most_common <- wc
    results$summary_table <- "Empty network. No general factors found."
    
    # Return consensus
    return(results)
    
    
  }
  
  # Apply Louvain
  communities <- lapply(1:consensus.iter, function(j){
    
    # igraph output
    output <- signed.louvain(network)
    
    # Obtain order
    if(order == "lower"){
      wc <- output$memberships[1,]
    }else if(order == "higher"){
      wc <- output$memberships[nrow(output$memberships),]
    }
    
    # Obtain text for progress
    progress_text <- paste0(
      "\r Consensus iteration ", formatC(
        j, digits = digits(consensus.iter) - 1,
        format = "d", flag = "0"
      ), " of ", consensus.iter,
      " complete."
    )
    
    # Print progress in message text
    cat(
      colortext(text = progress_text, defaults = "message")
    )
    
    # Close cat
    cat("\r")
    
    # Return
    return(wc)
    
  })
  
  # Simplify to a matrix
  df <- as.data.frame(
    t(simplify2array(communities, higher = FALSE))
  )
  
  # Remove communities
  rm(communities); gc(verbose = FALSE);
  
  # Obtain duplicate indices
  dupe_ind <- duplicated(df)
  
  # Rows for non-duplicates
  non_dupes <- data.frame(df[!dupe_ind,])
  
  # Rows for duplicates
  dupes <- data.frame(df[dupe_ind,])
  
  # Remove data frame 
  rm(df); gc(verbose = FALSE);
  
  # Match duplicates with non-duplicates
  dupe_count <- table(
    match(
      data.frame(t(dupes)), data.frame(t(non_dupes))
    ))
  
  # Remove dupes
  rm(dupes); gc(verbose = FALSE);
  
  # Obtain counts
  counts <- rep(1, nrow(non_dupes))
  counts[as.numeric(names(dupe_count))] <- counts[as.numeric(names(dupe_count))] + dupe_count
  
  # Change column names of non_dupes
  if(!is.null(colnames(network))){
    colnames(non_dupes) <- colnames(network)
  }
  
  # Set up summary table
  summary_table <- data.frame(
    N_Dimensions = apply(non_dupes, 1, function(x){
      length(na.omit(unique(x)))
    }),
    Proportion = as.matrix(counts / consensus.iter)
  )
  
  # Attach non-duplicate solutions
  summary_table <- cbind(summary_table, non_dupes)
  
  # Obtain max proportion
  wc_proportion <- unlist(summary_table[
    which.max(summary_table[,"Proportion"]),
    -c(1:2)
  ])
  
  # Set up results
  results <- list()
  results$most_common <- wc_proportion
  results$summary_table <- summary_table
  
  # Return consensus
  return(results)
}

# Lancichinetti & Fortunato (2012)
#' @noRd
# Most Common Consensus Clustering
# Updated 08.02.2023
most_common_consensus_top5 <- function(
    network,
    order = c("lower", "higher"),
    consensus.iter,
    resolution = 1
)
{
  
  # Obtain network names
  network_names <- colnames(network)
  
  # Check for empty network
  if(sum(network) == 0){
    
    # Return individual communities
    wc <- 1:ncol(network)
    
    # Assign names
    names(wc) <- network_names
    
    # Set up results
    results <- list()
    results$highest_modularity <- wc
    results$most_common <- wc
    results$iterative <- wc
    results$lowest_tefi <- wc
    results$summary_table <- "Empty network. No general factors found."
    
    # Return consensus
    return(results)
    
    
  }
  
  # Convert network to igraph
  igraph_network <- suppressWarnings(
    convert2igraph(abs(network))
  )
  
  # Ensure all nodes are included in igraph
  if(igraph::vcount(igraph_network) != ncol(network)){
    
    igraph_network <- igraph::add.vertices(
      igraph_network,
      nv = ncol(network) -
        igraph::vcount(igraph_network)
    )
    
  }
  
  # Apply Louvain
  communities <- lapply(1:consensus.iter, function(j, resolution){
    
    # igraph output
    output <- igraph::cluster_louvain(igraph_network, resolution = resolution)
    
    # Obtain memberships
    wc <- output$memberships
    
    # Check for no rows
    if(nrow(wc) == 0){
      wc <- output$membership
    }else{
      
      # Obtain order
      if(order == "lower"){
        wc <- wc[1,]
      }else if(order == "higher"){
        wc <- wc[nrow(wc),]
      }
      
    }
    
    # Obtain text for progress
    progress_text <- paste0(
      "\r Consensus iteration ", formatC(
        j, digits = digits(consensus.iter) - 1,
        format = "d", flag = "0"
      ), " of ", consensus.iter,
      " complete."
    )
    
    # Print progress in message text
    cat(
      colortext(text = progress_text, defaults = "message")
    )
    
    # Close cat
    cat("\r")
    
    # Return
    return(wc)
    
  }, resolution = resolution)
  
  # Simplify to a matrix
  wc_matrix <- t(simplify2array(communities, higher = FALSE))
  
  # Make data frame
  df <- as.data.frame(wc_matrix)
  
  # Obtain duplicate indices
  dupe_ind <- duplicated(df)
  
  # Rows for non-duplicates
  non_dupes <- data.frame(df[!dupe_ind,])
  
  # Rows for duplicates
  dupes <- data.frame(df[dupe_ind,])
  
  # Match duplicates with non-duplicates
  dupe_count <- table(
    match(
      data.frame(t(dupes)), data.frame(t(non_dupes))
    ))
  
  # Change column names of non_dupes
  if(!is.null(colnames(network))){
    colnames(non_dupes) <- colnames(network)
  }
  
  # Set up summary table
  summary_table <- data.frame(
    N_Dimensions = apply(non_dupes, 1, function(x){
      length(na.omit(unique(x)))
    }),
    Proportion = as.matrix(count(wc_matrix) / nrow(wc_matrix))
  )
  
  # Attach non-duplicate solutions
  summary_table <- cbind(summary_table, non_dupes)
  
  # Ordered table
  summary_table <- summary_table[
    order(summary_table$Proportion, decreasing = TRUE),
  ]
  
  # Check for number of rows in summary table
  if(nrow(summary_table) > 5){
    
    # Obtain top five solutions
    top5 <- summary_table[1:5,]
    
  }else if(nrow(summary_table) < 3){
    
    # Obtain max proportion
    wc_proportion <- unlist(summary_table[
      which.max(summary_table[,"Proportion"]),
      -c(1:2)
    ])
    
    # Set up results
    results <- list()
    results$most_common <- wc_proportion
    results$summary_table <- summary_table
    
  }else{
    
    # Obtain top solutions
    top5 <- summary_table
    
  }
  
  # Obtain proportions
  top5_proportions <- top5$Proportion
  
  # Obtain solutions (ensures matrix)
  top5_solutions <- as.matrix(
    top5[,-c(1,2)],
    ncol = ncol(top5) - 2
  )
  
  # Compute NMI over
  nmi_matrix <- matrix(
    0, nrow = length(top5_proportions),
    ncol = length(top5_proportions)
  )
  
  # Loop over to get NMI
  for(i in 1:nrow(nmi_matrix))
    for(j in 1:ncol(nmi_matrix)){
      
      # Lower triangle
      if(i > j){
        
        # Populate NMI matrix
        nmi_matrix[i,j] <- igraph::compare(
          top5_solutions[i,], top5_solutions[j,],
          method = "nmi"
        )
        
      }
      
    }
  
  # Make symmetric
  nmi_matrix <- nmi_matrix + t(nmi_matrix); diag(nmi_matrix) <- 1;
  
  # Matrix multiply by proportions
  wc_proportion <- top5_solutions[
    which.max(as.vector(nmi_matrix %*% top5_proportions)),
  ]
  
  # Set up results
  results <- list()
  results$most_common <- wc_proportion
  results$summary_table <- summary_table
  
  # Return consensus
  return(results)
}

# Lancichinetti & Fortunato (2012)
#' @noRd
# Most Common Consensus Clustering with TEFI Selection
# Updated 31.01.2023
most_common_tefi <- function(
    network, correlation,
    order = c("lower", "higher"),
    consensus.iter,
    resolution = 1
)
{
  
  # Obtain network names
  network_names <- colnames(network)
  
  # Check for empty network
  if(sum(network, na.rm = TRUE) == 0){
    
    # Return individual communities
    wc <- 1:ncol(network)
    
    # Assign names
    names(wc) <- network_names
    
    # Set up results
    results <- list()
    results$highest_modularity <- wc
    results$most_common <- wc
    results$iterative <- wc
    results$lowest_tefi <- wc
    results$summary_table <- "Empty network. No general factors found."
    
    # Return consensus
    return(results)
    
    
  }
  
  # Convert network to igraph
  igraph_network <- suppressWarnings(
    convert2igraph(abs(network))
  )
  
  # Ensure all nodes are included in igraph
  if(igraph::vcount(igraph_network) != ncol(network)){
    
    igraph_network <- igraph::add.vertices(
      igraph_network,
      nv = ncol(network) -
        igraph::vcount(igraph_network)
    )
    
  }
  
  # Apply Louvain
  communities <- lapply(1:consensus.iter, function(j, resolution){
    
    # igraph output
    output <- igraph::cluster_louvain(igraph_network, resolution = resolution)
    
    # Obtain memberships
    wc <- output$memberships
    
    # Check for no rows
    if(nrow(wc) == 0){
      wc <- output$membership
    }else{
      
      # Obtain order
      if(order == "lower"){
        wc <- wc[1,]
      }else if(order == "higher"){
        wc <- wc[nrow(wc),]
      }
      
    }
    
    # Return
    return(wc)
    
  }, resolution = resolution)
  
  # Simplify to a matrix
  wc_matrix <- t(simplify2array(communities, higher = FALSE))
  
  # Make data frame
  df <- as.data.frame(wc_matrix)
  
  # Obtain duplicate indices
  dupe_ind <- duplicated(df)
  
  # Rows for non-duplicates
  non_dupes <- data.frame(df[!dupe_ind,])
  
  # Rows for duplicates
  dupes <- data.frame(df[dupe_ind,])
  
  # Match duplicates with non-duplicates
  dupe_count <- table(
    match(
      data.frame(t(dupes)), data.frame(t(non_dupes))
    ))
  
  # Change column names of non_dupes
  if(!is.null(colnames(network))){
    colnames(non_dupes) <- colnames(network)
  }
  
  # Set up summary table
  summary_table <- data.frame(
    N_Dimensions = apply(non_dupes, 1, function(x){
      length(na.omit(unique(x)))
    }),
    Proportion = as.matrix(count(wc_matrix) / nrow(wc_matrix))
  )
  
  # Attach non-duplicate solutions
  summary_table <- cbind(summary_table, non_dupes)
  
  # Obtain total proportions
  ## Unique number of dimensions
  unique_dimensions <- unique(summary_table$N_Dimensions)
  
  ## Initialize total proportions
  total_proportions <- numeric(
    length(unique_dimensions)
  )
  names(total_proportions) <- unique_dimensions
  
  # Loop over unique dimensions to collect total proportions
  for(i in seq_along(total_proportions)){
    total_proportions[i] <- sum(
      summary_table$Proportion[
        summary_table$N_Dimensions == unique_dimensions[i]
      ],
      na.rm = TRUE
    )
  }
  
  # Maximum proportions (including tie)
  maximum_proportions <- max(total_proportions, na.rm = TRUE)
  
  # Indices for maximum proportions
  index_max_proportions <- total_proportions %in% maximum_proportions
  
  # Obtain dimensions
  maximum_dimensions <- as.numeric(names(total_proportions[index_max_proportions]))
  
  # Get solution indices
  solution_indices <- summary_table$N_Dimensions %in% maximum_dimensions
  
  # Get solutions
  solutions <- summary_table[
    solution_indices, -which(
      colnames(summary_table) %in% c("N_Dimensions", "Proportion")
    )
  ]
  
  # Check for more than one solution
  if(nrow(solutions) > 1){
    
    # Initialize TEFI vector
    tefi_vector <- numeric(nrow(solutions))
    
    # Populate TEFI vector
    for(i in seq_along(tefi_vector)){
      
      tefi_vector[i] <- tefi(
        data = abs(correlation),
        structure = unlist(solutions[i,])
      )$VN.Entropy.Fit
      
    }
    
    # Minimum TEFI
    minimum_tefi <- which.min(tefi_vector)
    
    # Final solution
    wc_proportion <- unlist(solutions[minimum_tefi,])
    
  }else{
    
    # Set most common
    wc_proportion <- unlist(solutions)
    
  }
  
  # Set up results
  results <- list()
  results$most_common <- wc_proportion
  results$summary_table <- summary_table
  
  # Return consensus
  return(results)
  
}

#%%%%%%%%%%%%%%%%%%%%
# NETWORKTOOLBOX ----
#%%%%%%%%%%%%%%%%%%%%

#' @noRd
# Cohen's d
# Updated 01.08.2021
d <- function(samp1, samp2)
{
  # Remove NAs
  samp1 <- samp1[!is.na(samp1)]
  samp2 <- samp2[!is.na(samp2)]
  
  # Means
  m1 <- mean(samp1, na.rm = TRUE)
  m2 <- mean(samp2, na.rm = TRUE)
  
  # Numerator
  num <- m1 - m2
  
  # Degrees of freedom
  df1 <- length(samp1) - 1
  df2 <- length(samp2) - 1
  
  # Variance
  var1 <- var(samp1, na.rm = TRUE)
  var2 <- var(samp2, na.rm = TRUE)
  
  # Denominator
  denom <- sqrt(
    ((df1 * var1) + (df2 * var2)) / (df1 + df2)
  )
  
  return(abs(num / denom))
  
}

# adapt.a
#' @noRd
#' @importFrom stats qchisq qf t.test var
# Adaptive Alpha
# Updated 01.08.2022
adapt.a <- function (test = c("anova","chisq","cor","one.sample","two.sample","paired"),
                     ref.n = NULL, n = NULL, alpha = .05, power = .80,
                     efxize = c("small","medium","large"), groups = NULL, df = NULL)
{
  
  # Need a test
  if(missing(test)){
    stop("test must be selected")
  }else{test <- match.arg(test)}
  
  # Assign medium effect size
  if(missing(efxize)){
    efxize <- "medium"
    message("No effect size selected. Medium effect size computed.")
  }else{efxize <- efxize}
  
  # ANOVA
  if(test == "anova"){
    
    # Check for groups
    if(is.null(groups)){
      stop("ANOVA is selected. Number of groups must be set")
    }
    
    # Set effect size
    efxize <- switch(
      efxize,
      "small" = 0.10,
      "medium" = 0.25,
      "large" = 0.40
    )
    
    # Determine reference sample size
    if(is.null(ref.n)){
      ref.n <- pwr::pwr.anova.test(f=efxize,power=power,sig.level=alpha,k=groups)$n
      message("ref.n is observations per group")
    }
    
    # Numerator
    num <- sqrt(ref.n*(log(ref.n)+qchisq((1-alpha),1)))
    
  }else if(test == "chisq"){ # Chi-square
    
    # Needs degrees of freedom
    if(is.null(df)){
      stop("Chi-square is selected. Degrees of freedom must be set")
    }
    
    # Set effect size
    efxize <- switch(
      efxize,
      "small" = 0.10,
      "medium" = 0.30,
      "large" = 0.50
    )

    # Determine reference sample size
    if(is.null(ref.n)){
      ref.n <- pwr::pwr.chisq.test(w=efxize,df=df,power=power,sig.level=alpha)$N
    }
    # Numerator
    num <- sqrt(ref.n*(log(ref.n)+qchisq((1-alpha),1)))
    
  }else if(test == "cor"){ # Correlation
    
    # Set effect size
    efxize <- switch(
      efxize,
      "small" = 0.10,
      "medium" = 0.30,
      "large" = 0.50
    )
    
    # Determine reference sample size
    if(is.null(ref.n)){
      ref.n <- pwr::pwr.r.test(r=efxize,power=power,sig.level=alpha)$n
    }
    
    # Numerator
    num <- sqrt(ref.n*(log(ref.n)+qchisq((1-alpha),1)))
    
  }else if(any(c("one.sample", "two.sample", "paired") %in% test)){# t-test
    
    # Set effect size
    efxize <- switch(
      efxize,
      "small" = 0.20,
      "medium" = 0.50,
      "large" = 0.80
    )
    
    # Determine reference sample size
    if(is.null(ref.n)){
      ref.n <- pwr::pwr.t.test(d=efxize,power=power,sig.level=alpha,type=test)$n
    }
    
    # Numerator
    num <- sqrt(ref.n*(log(ref.n)+qchisq((1-alpha),1)))
    
  }else{stop("test does not exist")}
  
  # Denominator
  denom <- (sqrt(n*(log(n)+qchisq((1-alpha),1))))
  
  # Adjusted alpha calculation
  adj.a <- alpha*num/denom
  
  # Critical values
  if(test == "anova"){
    
    critical.f <- function (groups, n, a)
    {
      df1 <- groups - 1
      df2 <- n - groups
      cvf <- qf(a, df1, df2, lower.tail = FALSE)
      return(cvf)
    }
    
    cv <- critical.f(groups, n, adj.a)
    
  }else if(test == "chisq"){
    
    critical.chi <- function (df, a)
    {
      cvchi <- qchisq(a, df, lower.tail = FALSE)
      return(cvchi)
    }
    
    cv <- critical.chi(df, adj.a)
    
  }else if(test == "cor"){
    
    critical.r <- function (n, a)
    {
      df <- n - 2
      critical.t <- qt( a/2, df, lower.tail = FALSE )
      cvr <- sqrt( (critical.t^2) / ( (critical.t^2) + df ) )
      return(cvr)
    }
    
    cv <- critical.r(n, adj.a)
    
  }else if(any(c("one.sample", "two.sample", "paired") %in% test)){
    
    critical.t <- function (n, a)
    {
      df <- n - 2
      cvt <- qt( a/2, df, lower.tail = FALSE )
      return(cvt)
    }
    
    cv <- critical.t(n, adj.a)
    
  }
  
  # Output
  output <- list(
    adapt.a = adj.a, crit.value = cv,
    orig.a = alpha, ref.n = ref.n,
    exp.n = n, power = power,
    efxize = efxize
  )
  # Check for ANOVA or Chi-square
  if(test == "anova"){
    output$groups <- groups
    output$df <- c((groups - 1), (n - groups))
    
  }else if(test=="chisq"){
    output$df <- df
  }
  # Add test
  output$test <- test
  
  return(output)
}

# binarize
#' @noRd
# Make adjacency matrix unweighted
# Updated 30.12.2021
binarize <- function (A)
{
  A <- as.matrix(A)
  bin <- ifelse(A!=0,1,0)
  row.names(bin) <- row.names(A)
  colnames(bin) <- colnames(A)

  return(bin)
}

# conn
#' @noRd
# Network Connectivity
# Updated 30.12.2021
conn <- function (A)
{
  diag(A)<-0

  weights<-0
  wc<-0
  B<-A[lower.tri(A)]
  for(i in 1:length(B))
    if (B[i]!=0)
    {
      wc <- wc+1
      weights[wc] <- B[i]
    }
  tot<-sum(weights)
  mea<-mean(weights)
  s<-sd(weights)

  possible<-sum(ifelse(A!=0,1,0)/2)
  den<-possible/((ncol(A)^2-ncol(A))/2)

  return(list(weights=weights,mean=mea,sd=s,total=tot,density=den))
}

# clustcoeff
#' @noRd
# Clustering Coefficient
# Updated 30.12.2021
clustcoeff <- function (A, weighted = FALSE)
{
  diag(A) <- 0

  if(nrow(A)!=ncol(A))
  {stop("Input not an adjacency matrix")}
  if(!weighted)
  {n<-ncol(A)
  A<-ifelse(A!=0,1,0)
  C<-matrix(0,nrow=n,ncol=1)

  for(i in 1:n)
  {
    v<-which(A[i,]!=0)
    k<-length(v)
    if(k >= 2)
    {
      S<-A[v,v]
      C[i]<-sum(S)/(k^2-k)
    }}

  C <- round(as.vector(C),3)
  names(C)<-colnames(A)

  CCi<-C
  CC <- mean(C)

  }else{
    K<-colSums(A!=0)
    m<-A^(1/3)
    cyc<-diag(m%*%m%*%m)
    K[cyc==0]<-Inf
    C <- as.vector(round(cyc/(K*(K-1)),3))
    names(C)<-colnames(A)
    CCi<-C
    CC<-mean(C)
  }

  return(list(CC=CC,CCi=CCi))
}

# degree
#' @noRd
# Degree
# Updated 30.12.2021
degree <- function (A)
{
  if(nrow(A)!=ncol(A))
  {stop("Input not an adjacency matrix")}

  A <- as.matrix(A)

  A <- binarize(A)

  if(isSymmetric(A, check.attributes = FALSE))
  {
    Deg <- as.vector(colSums(A))
    names(Deg) <- colnames(A)
    return(Deg)
  }else
  {
    #In-degree
    inDeg <- as.vector(colSums(A, na.rm = TRUE))
    names(inDeg) <- colnames(A)

    #Out-degree
    outDeg <- as.vector(rowSums(A, na.rm = TRUE))
    names(outDeg) <- colnames(A)

    #Relative influence
    relinf <- as.vector((outDeg-inDeg)/(outDeg+inDeg))
    names(relinf) <- colnames(A)

    #Reciprocal degree
    for(i in 1:nrow(A))
      for(j in 1:ncol(A))
      {
        A[i,j] <- ifelse(A[i,j] == 1 & A[j,i] == 1, 1, 0)
        A[j,i] <- ifelse(A[i,j] == 1 & A[j,i] == 1, 1, 0)
      }

    recipDeg <- colSums(A, na.rm = TRUE)
    names(recipDeg) <- colnames(A)

    if(all(relinf<.001))
    {Deg <- as.vector(inDeg)
    names(Deg) <- colnames(A)
    return(Deg)
    }else{return(list(inDegree=inDeg,
                      outDegree=outDeg,
                      recipDegree=recipDeg,
                      relInf=relinf))}
  }
}

# distance
#' @noRd
# Distance between nodes
# Updated 30.12.2021
distance<-function (A, weighted = FALSE)
{
  if(nrow(A)!=ncol(A))
  {stop("Input not an adjacency matrix")}
  if(!weighted)
  {B<-ifelse(A!=0,1,0)
  l<-1
  Lpath<-B
  D<-B
  Idx<-matrix(TRUE,nrow=nrow(B),ncol=ncol(B))
  while(any(Idx))
  {
    l<-l+1
    Lpath<-(as.matrix(Lpath))%*%(as.matrix(B))
    for(e in 1:nrow(Lpath))
      for(w in 1:ncol(Lpath))
        Idx[e,w]<-(Lpath[e,w]!=0&&(D[e,w]==0))
    D[Idx]<-l
  }

  D[!D]<-Inf
  diag(D)<-0
  }else if(weighted){
    G<-ifelse(1/A==Inf,0,1/A)

    if(any(G==-Inf))
    {G<-ifelse(G==-Inf,0,G)}

    if(any(!G==t(G)))
    {if(max(abs(G-t(G)))<1e-10)
    {G<-(G+G)/2}}

    n<-ncol(G)
    D<-matrix(Inf,nrow=n,ncol=n)
    diag(D)<-0
    B<-matrix(0,nrow=n,ncol=n)

    for(u in 1:n)
    {
      S<-matrix(TRUE,nrow=n,ncol=1)
      L1<-G
      V<-u
      while(TRUE)
      {
        S[V]<-0
        L1[,V]<-0
        for(v in V)
        {
          W<-which(L1[v,]!=0)
          d<-apply(rbind(D[u,W],(D[u,v]+L1[v,W])),2,min)
          wi<-apply(rbind(D[u,W],(D[u,v]+L1[v,W])),2,which.min)
          D[u,W]<-d
          ind<-W[wi==2]
          B[u,ind]<-B[u,v]+1
        }

        minD<-suppressWarnings(min(D[u,S==TRUE]))
        if(length(minD)==0||is.infinite(minD)){break}

        V<-which(D[u,]==minD)
      }
    }
  }

  D<-ifelse(D==Inf,0,D)

  colnames(D)<-colnames(A)
  row.names(D)<-colnames(A)
  return(D)
}

# lattnet
#' @noRd
# Generate lattice network
# Updated 30.12.2021
lattnet <- function (nodes, edges)
{
  dlat<-matrix(0,nrow=nodes,ncol=nodes)
  lat<-matrix(0,nrow=nodes,ncol=nodes)

  balance <- sum(lat) - edges

  count <- 0

  while(sign(balance) == -1){

    if(count == 0){

      for(i in 1:nodes){

        if(i != nodes){
          dlat[i, (i + 1)] <- 1
        }
      }

    }else{

      for(i in 1:nodes){

        if(i < (nodes - count)){
          dlat[i, (i + (count + 1))] <- 1
        }

      }

    }

    count <- count + 1

    balance <- sum(dlat) - edges

  }

  over <- sum(dlat) - edges

  if(over != 0){

    rp <- sample(which(dlat==1), over, replace = FALSE)

    dlat[rp] <- 0

  }

  lat <- dlat + t(dlat)

  return(lat)
}

# LoGo
#' @noRd
# Local-Global Sparse Inverse Covariance Matrix
# Updated 30.12.2021
LoGo <- function (cormat, cliques, separators,
                  partial = TRUE, ...)
{
  #covariance matrix
  standardize <- TRUE


  S <- cormat

  if(missing(separators))
  {separators<-NULL}

  if(missing(cliques))
  {cliques<-NULL}


  if(is.null(separators)&is.null(cliques))
  {
    tmfg<-TMFG(cormat)
    separators<-tmfg$separators
    cliques<-tmfg$cliques
  }

  n<-ncol(S)
  Jlogo<-matrix(0,nrow=n,ncol=n)

  if(!is.list(cliques)&!is.list(separators))
  {
    for(i in 1:nrow(cliques))
    {
      v<-cliques[i,]
      Jlogo[v,v]<-Jlogo[v,v]+solve(S[v,v])
    }

    for(i in 1:nrow(separators))
    {
      v<-separators[i,]
      Jlogo[v,v]<-Jlogo[v,v]-solve(S[v,v])
    }
  }else{
    for(i in 1:length(cliques))
    {
      v<-cliques[[i]]
      Jlogo[v,v]<-Jlogo[v,v]+solve(S[v,v])
    }

    for(i in 1:length(separators))
    {
      v<-separators[[i]]
      Jlogo[v,v]<-Jlogo[v,v]-solve(S[v,v])
    }
  }

  if(partial)
  {
    Jlogo<-(-cov2cor(Jlogo))
    if(any(is.na(Jlogo)))
    {Jlogo <- ifelse(is.na(Jlogo),0,Jlogo)}
    diag(Jlogo)<-0
  }

  colnames(Jlogo)<-colnames(data)
  row.names(Jlogo)<-colnames(data)

  if(!isSymmetric(Jlogo))
  {Jlogo[lower.tri(Jlogo)] <- Jlogo[upper.tri(Jlogo)]}

  return(logo=Jlogo)
}

# pathlengths
#' @noRd
# Path Lengths
# Updated 30.12.2021
pathlengths <- function (A, weighted = FALSE)
{
  if(nrow(A)!=ncol(A))
  {stop("Input not an adjacency matrix")}
  if(!weighted)
  {D<-distance(A,weighted=FALSE)}else if(weighted){D<-distance(A,weighted=TRUE)}
  n<-nrow(D)
  for(i in 1:ncol(D))
    for(j in 1:nrow(D))
      if(is.infinite(D[j,i]))
      {D[j,i]<-0}
  if(any(colSums(D)==0))
  {
    newD <- D
    newD[,(which(colSums(D)==0))] <- rep(NA,length(which(colSums(D)==0)))
  }else{newD <- D}

  aspli<-colSums(newD*(newD!=0))/(ncol(newD)-1)
  aspl<-mean(aspli,na.rm=TRUE)

  Emat<-(D*(D!=0))

  ecc<-matrix(nrow=nrow(Emat),ncol=1)

  for(i in 1:nrow(Emat))
  {ecc[i,]<-max(Emat[i,])}

  d<-max(ecc)

  ecc <- as.vector(ecc)
  names(ecc) <- colnames(A)

  aspli <- as.vector(aspli)
  names(aspli) <- colnames(A)

  return(list(ASPL=aspl,ASPLi=aspli,ecc=ecc,diameter=d))
}

# pwr.r.test
#' @noRd
#' @importFrom stats uniroot
# Power for correlations from pwr 1.3.0
# Updated 30.12.2021
pwr.r.test <- function (n = NULL, r = NULL, sig.level = 0.05, power = NULL,
                        alternative = c("two.sided", "less", "greater"))
{
  if (sum(sapply(list(n, r, power, sig.level), is.null)) !=
      1)
    stop("exactly one of n, r, power, and sig.level must be NULL")
  if (!is.null(r) && is.character(r))
    r <- cohen.ES(test = "r", size = r)$effect.size
  if (!is.null(sig.level) && !is.numeric(sig.level) || any(0 >
                                                           sig.level | sig.level > 1))
    stop(sQuote("sig.level"), " must be numeric in [0, 1]")
  if (!is.null(power) && !is.numeric(power) || any(0 > power |
                                                   power > 1))
    stop(sQuote("power"), " must be numeric in [0, 1]")
  if (!is.null(n) && any(n < 4))
    stop("number of observations must be at least 4")
  alternative <- match.arg(alternative)
  tside <- switch(alternative, less = 1, two.sided = 2, greater = 3)
  if (tside == 2 && !is.null(r))
    r <- abs(r)
  if (tside == 3) {
    p.body <- quote({
      ttt <- qt(sig.level, df = n - 2, lower = FALSE)
      rc <- sqrt(ttt^2/(ttt^2 + n - 2))
      zr <- atanh(r) + r/(2 * (n - 1))
      zrc <- atanh(rc)
      pnorm((zr - zrc) * sqrt(n - 3))
    })
  }
  if (tside == 1) {
    p.body <- quote({
      r <- -r
      ttt <- qt(sig.level, df = n - 2, lower = FALSE)
      rc <- sqrt(ttt^2/(ttt^2 + n - 2))
      zr <- atanh(r) + r/(2 * (n - 1))
      zrc <- atanh(rc)
      pnorm((zr - zrc) * sqrt(n - 3))
    })
  }
  if (tside == 2) {
    p.body <- quote({
      ttt <- qt(sig.level/2, df = n - 2, lower = FALSE)
      rc <- sqrt(ttt^2/(ttt^2 + n - 2))
      zr <- atanh(r) + r/(2 * (n - 1))
      zrc <- atanh(rc)
      pnorm((zr - zrc) * sqrt(n - 3)) + pnorm((-zr - zrc) *
                                                sqrt(n - 3))
    })
  }
  if (is.null(power))
    power <- eval(p.body)
  else if (is.null(n))
    n <- uniroot(function(n) eval(p.body) - power, c(4 +
                                                       1e-10, 1e+09))$root
  else if (is.null(r)) {
    if (tside == 2)
      r <- uniroot(function(r) eval(p.body) - power, c(1e-10,
                                                       1 - 1e-10))$root
    else r <- uniroot(function(r) eval(p.body) - power, c(-1 +
                                                            1e-10, 1 - 1e-10))$root
  }
  else if (is.null(sig.level))
    sig.level <- uniroot(function(sig.level) eval(p.body) -
                           power, c(1e-10, 1 - 1e-10))$root
  else stop("internal error")
  METHOD <- "approximate correlation power calculation (arctangh transformation)"
  structure(list(n = n, r = r, sig.level = sig.level, power = power,
                 alternative = alternative, method = METHOD), class = "power.htest")
}

# cohen.ES
#' @noRd
# Effect sizes for correlations from pwr 1.3.0
# Updated 30.12.2021
cohen.ES <- function(test=c("p","t","r","anov","chisq","f2"),size=c("small","medium","large")){
  test <- match.arg(test)
  size <- match.arg(size)
  ntest <- switch(test, p = 1, t = 2,r=3,anov=4,chisq=5,f2=6)
  if(ntest==1){
    ES<-switch(size,small=0.2,medium=0.5,large=0.8)
  }
  if(ntest==2){
    ES<-switch(size,small=0.2,medium=0.5,large=0.8)
  }
  if(ntest==3){
    ES<-switch(size,small=0.1,medium=0.3,large=0.5)
  }
  if(ntest==4){
    ES<-switch(size,small=0.1,medium=0.25,large=0.4)
  }
  if(ntest==5){
    ES<-switch(size,small=0.1,medium=0.3,large=0.5)
  }
  if(ntest==6){
    ES<-switch(size,small=0.02,medium=0.15,large=0.35)
  }

  METHOD <- "Conventional effect size from Cohen (1982)"
  structure(list(test = test,size=size,effect.size=ES,
                 method = METHOD), class = "power.htest")
}

# randnet
#' @noRd
# Generate random network
# Updated 02.07.2022
randnet <- function (nodes = NULL, edges = NULL, A = NULL)
{
  if(is.null(A)){
    
    # Initialize matrix
    mat <- matrix(1, nrow = nodes, ncol = nodes)

    # Set diagonal to zero
    diag(mat) <- 0

    # Indices of upper diagonal
    ind <- ifelse(upper.tri(mat) == TRUE, 1, 0)
    i <- which(ind == 1)

    # Sample zeros and ones
    rp <- sample(length(i))
    # Get indices
    irp <- i[rp]

    # Initialize random matrix
    rand <- matrix(0, nrow = nodes, ncol = nodes)

    # Insert edges
    rand[irp[1:edges]] <- 1

    # Make symmetric
    rand <- rand + t(rand)

  }else{

    # Make diagonal of network zero
    diag(A) <- 0

    # Compute degree
    degrees <- degree(A)
    
    # Identify disconnected nodes
    if(any(degrees == 0)){
      disconnected <- which(degrees == 0)
      degrees <- degrees[-disconnected]
    }

    # Get degrees based on directed or undirected
    # Use igraph
    if(is.list(degrees)){
      rand <- as.matrix(igraph::as_adj(igraph::sample_degseq(out.deg = degrees$outDegree, in.deg = degrees$inDegree, method = "vl")))
    }else{
      rand <- as.matrix(igraph::as_adj(igraph::sample_degseq(out.deg = degrees, method = "vl")))
    }
    
    # Add back disconnected nodes
    if(exists("disconnected", envir = environment())){
      
      # New random matrix
      new_rand <- matrix(0, nrow = ncol(A), ncol = ncol(A))
      
      # Insert old random matrix into new random matrix
      new_rand[-disconnected, -disconnected] <- rand
      
      # Copy new random matrix into old
      rand <- new_rand
      
    }
    
    # Apply back row and column names
    row.names(rand) <- colnames(A)
    colnames(rand) <- colnames(A)
    
  }

  return(rand)
}

# strength
#' @noRd
# Node Strength
# Updated 30.12.2021
strength <- function (A, absolute = TRUE)
{
  if(is.vector(A))
  {return(0)
  }else if(nrow(A)!=ncol(A))
  {stop("Input not an adjacency matrix")}

  if(absolute)
  {A <- abs(A)}
  A <- as.matrix(A)

  if(isSymmetric(A, check.attributes = FALSE))
  {
    Str <- round(as.vector(colSums(A)),2)
    names(Str) <- colnames(A)
    return(Str)
  }else{
    #In-strength
    inStr <- as.vector(colSums(A))
    names(inStr) <- colnames(A)
    #Out-strength
    outStr <- as.vector(rowSums(A))
    names(outStr) <- colnames(A)
    #Relative influence
    relinf <- as.vector((outStr-inStr)/(outStr+inStr))
    names(relinf) <- colnames(A)

    if(all(relinf<.001))
    {Str <- round(as.vector(colSums(A)),2)
    names(Str) <- colnames(A)
    return(Str)
    }else{return(list(inStrength=inStr,outStrength=outStr,relInf=relinf))}
  }
}

# smallworldness
#' @noRd
# Small-worldness measures
# Updated 30.12.2021
smallworldness <- function (A, iter = 100, progBar = FALSE, method = c("HG","rand","TJHBL"))
{
  if(missing(method))
  {method<-"TJHBL"
  }else{method<-match.arg(method)}

  mat<-matrix(0,nrow=nrow(A),ncol=ncol(A)) #Initialize bootstrap matrix
  asamps<-matrix(0,nrow=iter) #Initialize sample matrix
  csamps<-matrix(0,nrow=iter) #Initialize sample matrix
  if(progBar)
  {pb <- txtProgressBar(max=iter, style = 3)}
  for(i in 1:iter) #Generate array of bootstrapped samples
  {
    rand<-randnet(A = A)
    if(method=="TJHBL")
    {latt<-lattnet(ncol(A),sum(ifelse(A!=0,1,0))/2)}
    asamps[i,]<-pathlengths(rand)$ASPL
    if(method=="rand")
    {csamps[i,]<-clustcoeff(rand)$CC
    }else if(method=="HG"){csamps[i,]<-transitivity(rand)
    }else if(method=="TJHBL"){csamps[i,]<-clustcoeff(latt)$CC}else{stop("Method not available")}
    if(progBar)
    {setTxtProgressBar(pb, i)}
  }
  if(progBar)
  {close(pb)}

  nodes<-ncol(A)
  ASPL<-pathlengths(A)$ASPL
  CC<-clustcoeff(A)$CC
  trans<-transitivity(A)
  rASPL<-mean(asamps)

  if(method=="rand")
  {rCC<-mean(csamps)
  swm<-(CC/rCC)/(ASPL/rASPL)
  lrCCt<-rCC
  }else if(method=="HG")
  {rtrans<-mean(csamps)
  swm<-(trans/rtrans)/(ASPL/rASPL)
  lrCCt<-rtrans
  }else if(method=="TJHBL")
  {lCC<-mean(csamps)
  swm<-(rASPL/ASPL)-(CC/lCC)
  lrCCt<-lCC
  }

  return(list(swm=swm, rASPL=rASPL, lrCCt=lrCCt))
}

# transitivity
#' @noRd
# Transitivity
# Updated 30.12.2021
transitivity <- function (A, weighted = FALSE)
{
  if(!weighted)
  {
    A<-ifelse(A!=0,1,0)
    trans<-sum(diag(A%*%A%*%A))/((sum(A%*%A))-sum(diag(A%*%A)))
  }else if(weighted){
    K<-colSums(ifelse(A!=0,1,0))
    W<-A^(1/3)
    cyc<-diag(W%*%W%*%W)
    trans<-sum(cyc)/sum(K*(K-1))
  }

  return(trans)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# WEIGHTED TOPOLOGICAL OVERLAP ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# From wTO 1.6.3
#' @noRd
#'
# Weighted topological overlap
wTO <- function (A, sign = c("abs", "sign"))
{
  if (sign %in% c("abs", "absolute")) {
    A = abs(A)
  }
  A_TF = as.data.frame(subset(A, select = row.names(A)))
  C = as.matrix(A) %*% t(A)
  W = C + A_TF
  K = matrix(NA, nrow(A_TF), ncol(A_TF))
  KI = rowSums(abs(A), na.rm = T)
  for (ii in 1:nrow(A_TF)) {
    for (jj in 1:ncol(A_TF)) {
      K[ii, jj] = min(KI[ii], KI[jj])
    }
  }
  WTO = round(W/(K + 1 - abs(A_TF)), 3)
  return(WTO)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%
# NETWORK DESCRIPTIVES ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%

# From WGCNA version 1.70-3
#' @noRd
#' @importFrom graphics hist
#' @importFrom stats lm
#'
#Scale-free fit index
#Updated 12.05.2021
scaleFreeFitIndex=function(k,nBreaks=10, removeFirst = FALSE)
{
  discretized.k = cut(k, nBreaks)
  dk = tapply(k, discretized.k, mean)
  p.dk = as.vector(tapply(k, discretized.k, length)/length(k))
  breaks1 = seq(from = min(k), to = max(k),
                length = nBreaks + 1)
  hist1 = hist(k, breaks = breaks1, plot = FALSE, right = TRUE)
  dk2 = hist1$mids
  dk = ifelse(is.na(dk), dk2, dk)
  dk = ifelse(dk == 0, dk2, dk)
  p.dk = ifelse(is.na(p.dk), 0, p.dk)
  log.dk = as.vector(log10(dk))
  if (removeFirst) {
    p.dk = p.dk[-1]
    log.dk = log.dk[-1]
  }
  log.p.dk= as.numeric(log10(p.dk + 1e-09))
  lm1 = try(lm(log.p.dk ~ log.dk));
  if (inherits(lm1, "try-error")) browser();
  lm2 = lm(log.p.dk ~ log.dk + I(10^log.dk))
  datout=data.frame(Rsquared.SFT=summary(lm1)$r.squared,
                    slope.SFT=summary(lm1)$coefficients[2, 1],
                    truncatedExponentialAdjRsquared= summary(lm2)$adj.r.squared)
  datout
}

#%%%%%%%%%%%%%%%%%
# ENTROPY FIT ----
#%%%%%%%%%%%%%%%%%

#' @noRd
# Mimics count from plyr::count
# Updated 24.02.2022
count <- function(data)
{
  # Make data frame
  df <- as.data.frame(data)

  # Obtain duplicate indices
  dupe_ind <- duplicated(df)

  # Rows for non-duplicates
  non_dupes <- data.frame(df[!dupe_ind,])

  # Rows for duplicates
  dupes <- data.frame(df[dupe_ind,])

  # Match duplicates with non-duplicates
  dupe_count <- table(
    match(
      data.frame(t(dupes)), data.frame(t(non_dupes))
    ))

  # Obtain counts
  counts <- rep(1, nrow(non_dupes))
  counts[as.numeric(names(dupe_count))] <- counts[as.numeric(names(dupe_count))] + dupe_count

  return(counts)
  
}

#%%%%%%%%%
# EGA ----
#%%%%%%%%%

#' @noRd
# Unidimensionality check
# Updated 18.10.2022
unidimensionality.check <- function(
    data,
    n,
    corr,
    correlation, 
    uni.method = c("louvain", "LE", "expand"),
    model, model.args,
    algorithm, algorithm.args,
    consensus.method, consensus.iter
)
{
  
  # Make unidimensional method lowercase
  uni.method <- tolower(uni.method)
  
  # Perform algorithm
  if(uni.method == "louvain"){
    
    # Most common consensus with Louvain
    wc <- most_common_consensus(
      network = abs(correlation),
      order = "higher",
      consensus.iter = consensus.iter,
      resolution = 0.95
    )$most_common
    
  }else if(uni.method == "le"){
    
    
    # Try Leading Eigenvalue
    wc <- try(
      igraph::cluster_leading_eigen(
        convert2igraph(abs(correlation))
      )$membership,
      silent = TRUE
    )
    
    # If error, then use Louvain
    if(any(class(wc) == "try-error")){
      
      # Most common consensus with Louvain
      wc <- most_common_consensus(
        network = abs(correlation),
        order = "higher",
        consensus.iter = consensus.iter,
        resolution = 0.95
      )$most_common
      
      warning("Error occurred in Leading Eigenvalue algorithm. Using \"louvain\" for unidimensional check")
      
    }
    
  }else if(uni.method == "expand"){
    
    # Check for {igraph} algorithm
    if(is.function(algorithm)){
      
      # Identify whether algorithm is Spinglass
      if("spins" %in% methods::formalArgs(algorithm)){
        
        # Check for whether data are a correlation matrix
        if(isSymmetric(unname(as.matrix(data)))){
          data <- MASS_mvrnorm(n = n, mu = rep(0, ncol(data)), Sigma = data)
        }
        
        # Simulate data from unidimensional factor model
        sim.data <- sim.func(data = data, nvar = 4, nfact = 1, load = .70)
        
        ## Compute correlation matrix
        correlation <- suppressMessages(
          switch(corr,
                 "cor_auto" = qgraph::cor_auto(sim.data, forcePD = TRUE),
                 "pearson" = cor(sim.data, use = "pairwise.complete.obs"),
                 "spearman" = cor(sim.data, method = "spearman", use = "pairwise.complete.obs")
          )
        )
        
      }else{## Expand correlation matrix
        correlation <- expand.corr(correlation)
      }
      
    }else{## Expand correlation matrix
      correlation <- expand.corr(correlation)
    }
    
    # Unidimensional result
    wc <- EGA.estimate(
      data = correlation, n = n,
      model = model, model.args = model.args,
      algorithm = algorithm, algorithm.args = algorithm.args,
      consensus.method = consensus.method, consensus.iter = consensus.iter
    )$wc
    
  }
  
  # Assign names
  names(wc) <- colnames(correlation)
  
  # Collect dimensions
  n.dim <- length(na.omit(unique(wc)))
  
  # Populate results
  results <- list(
    wc = wc,
    n.dim = n.dim
  )
  
  # Return results
  return(results)
  
}

#%%%%%%%%%%%%%%
# PLOTTING ----
#%%%%%%%%%%%%%%

#' @noRd
# Sub-routine for compare.EGA.plots
# Updated 05.06.2021
compare.EGA <- function(ega.object1, ega.object2)
{
  # Plots
  plot1 <- ega.object1
  plot2 <- ega.object2

  # Reorder node coordinates for plot2
  plot2$data <- plot2$data[row.names(plot1$data),]

  # Reorder edge coordinates for plot2
  for(i in 1:nrow(plot2$layers[[1]]$data)){

    plot2$layers[[1]]$data$X1[i] <- which(plot2$layers[[1]]$data$X1[i] == plot2$data$x)
    plot2$layers[[1]]$data$X2[i] <- which(plot2$layers[[1]]$data$X2[i] == plot2$data$x)
    plot2$layers[[1]]$data$Y1[i] <- which(plot2$layers[[1]]$data$Y1[i] == plot2$data$y)
    plot2$layers[[1]]$data$Y2[i] <- which(plot2$layers[[1]]$data$Y2[i] == plot2$data$y)

  }

  # Reassign edge coordinates based on plot1
  plot2$layers[[1]]$data$X1 <- plot1$data$x[plot2$layers[[1]]$data$X1] # X1
  plot2$layers[[1]]$data$X2 <- plot1$data$x[plot2$layers[[1]]$data$X2] # X2
  plot2$layers[[1]]$data$Y1 <- plot1$data$y[plot2$layers[[1]]$data$Y1] # Y1
  plot2$layers[[1]]$data$Y2 <- plot1$data$y[plot2$layers[[1]]$data$Y2] # Y2

  # Assign coordinates of plot1 to plot2
  plot2$data$x <- plot1$data$x
  plot2$data$y <- plot1$data$y

  # Make plot list
  plots.net <- list()
  plots.net[[1]] <- plot1
  plots.net[[2]] <- plot2

  # Return plot list
  return(plots.net)

}

#' @noRd
# Defaults for GGally plotting
# For plots and methods
# Updated 28.07.2021
GGally.args <- function(plot.args)
{
  default.args <- formals(GGally::ggnet2)
  ega.default.args <- list(node.size = 12, edge.size = 8,
                           alpha = 0.5, label.size = 5,
                           edge.alpha = 0.5, layout.exp = 0.2)
  default.args[names(ega.default.args)]  <- ega.default.args
  default.args <- default.args[-length(default.args)]


  if("node.alpha" %in% names(plot.args)){
    plot.args$alpha <- plot.args$node.alpha
    plot.args$node.alpha <- NULL
  }

  if("vsize" %in% names(plot.args)){
    plot.args$node.size <- plot.args$vsize
    plot.args$vsize <- NULL
  }

  if("legend.names" %in% names(plot.args)){
    legend.names <- plot.args$legend.names
    plot.args$legend.names <- NULL
  }

  if(!"color.palette" %in% names(plot.args)){
    default.args$color.palette <- "polychrome"
  }

  if("color.palette" %in% names(plot.args)){

    if(length(plot.args$color.palette) == 1){
      
      if(tolower(plot.args$color.palette) == "greyscale" |
         tolower(plot.args$color.palette) == "grayscale" | 
         tolower(plot.args$color.palette) == "colorblind"){
        plot.args$edge.color <- c("#293132", "grey25")
        plot.args$edge.lty <- c("solid", "dashed")
      }
      
    }

  }

  if(!"edge.color" %in% names(plot.args)){
    plot.args$edge.color <- c("darkgreen", "red")
  }else if(length(plot.args$edge.color) != 2){
    stop("Two colors needed for 'edge.color'")
  }

  if(!"edge.lty" %in% names(plot.args)){
    plot.args$edge.lty <- c("solid", "solid")
  }else if(length(plot.args$edge.lty) != 2){
    stop("Two line types needed for 'edge.lty'")
  }


  if(any(names(plot.args) %in% names(default.args))){
    target.args <- plot.args[which(names(plot.args) %in% names(default.args))]
    default.args[names(target.args)] <- target.args
  }

  plot.args <- default.args

  return(plot.args)
}

#' @importFrom graphics text
#' @noRd
# Color sorting for EGA palettes
# For EGA_color_palette
# Updated 17.12.2020
color.sort <- function (wc)
{
  unlist(lapply(sort(wc, na.last = TRUE), function(x, uniq){

    if(is.na(x)){
      NA
    }else{
      which(x == uniq)
    }
  }, uniq = sort(unique(wc), na.last = TRUE)))
}

#' @noRd
# Rescale edges for GGally
# For plots
# Updated 17.01.2021
rescale.edges <- function (network, size)
{
  # Set rescaling
  scaling <- seq(0, 1, .000001) * size
  names(scaling) <- seq(0, 1, .000001)

  # Vectorize edges
  edges <- round(as.vector(as.matrix(network)), 5)

  # Get absolute edges
  abs.edges <- abs(edges)

  # Get edge signs
  signs.edges <- sign(edges)

  # Rescale edges
  rescaled.edges <- unname(scaling[as.character(abs.edges)])

  return(rescaled.edges)
}

#' @noRd
# Compare plots fix
# Updated 28.03.2022
compare.plot.fix.EGA <- function(object.list,
                                 plot.args = list()){


  ## Original plot arguments
  original.plot.args <- plot.args

  ## Check for input plot arguments
  if("legend.names" %in% names(plot.args)){
    legend.names <- plot.args$legend.names
  }
  plot.args <- GGally.args(plot.args)
  color.palette <- plot.args$color.palette

  ## Initialize plot list
  ega.plots <- list()

  # Loop through object list
  for(i in 1:length(object.list)){

    if(is(object.list[[i]], "EGA")){
      x <- object.list[[i]]
    }else if(is(object.list[[i]], "bootEGA")){
      x <- list(
        network = object.list[[i]]$typicalGraph$graph,
        wc = object.list[[i]]$typicalGraph$wc
      )
    }else if(is(object.list[[i]], "dynEGA")){
      x <- object.list[[i]]$dynEGA
    }

    ### Plot ###
    if(i != 1){
      ## Reset plot arguments
      plot.args <- original.plot.args
      
      ## Check for input plot arguments
      if("legend.names" %in% names(plot.args)){
        legend.names <- plot.args$legend.names
      }
      plot.args <- GGally.args(plot.args)
      color.palette <- plot.args$color.palette
    }
    
    # Insignificant values (keeps ggnet2 from erroring out)
    x$network <- ifelse(abs(as.matrix(x$network)) <= .00001, 0, as.matrix(x$network))
    
    if(exists("legend.names")){
      for(l in 1:length(unique(legend.names))){
        x$wc[x$wc == l] <- legend.names[l]
      }
    }
    
    # Reorder network and communities
    if(i == 1){
      x$network <- x$network[order(x$wc), order(x$wc)]
      x$wc <- x$wc[order(x$wc)]
      fix.order.wc <- names(x$wc)
    }else{
      x$network <- x$network[fix.order.wc, fix.order.wc]
      x$wc <- x$wc[fix.order.wc]
    }
    
    # weighted  network
    network1 <- network::network(x$network,
                                 ignore.eval = FALSE,
                                 names.eval = "weights",
                                 directed = FALSE)
    
    network::set.vertex.attribute(network1, attrname= "Communities", value = x$wc)
    network::set.vertex.attribute(network1, attrname= "Names", value = network::network.vertex.names(network1))
    network::set.edge.attribute(network1, "color", ifelse(network::get.edge.value(network1, "weights") > 0, plot.args$edge.color[1], plot.args$edge.color[2]))
    network::set.edge.attribute(network1, "line", ifelse(network::get.edge.value(network1, "weights") > 0, plot.args$edge.lty[1], plot.args$edge.lty[2]))
    network::set.edge.value(network1,attrname="AbsWeights",value=abs(x$network))
    network::set.edge.value(network1,attrname="ScaledWeights",
                            value=matrix(rescale.edges(x$network, plot.args$edge.size),
                                         nrow = nrow(x$network),
                                         ncol = ncol(x$network)))
    
    # Layout "Spring"
    graph1 <- convert2igraph(x$network)
    edge.list <- igraph::as_edgelist(graph1)
    layout.spring <- qgraph::qgraph.layout.fruchtermanreingold(edgelist = edge.list,
                                                               weights =
                                                                 abs(igraph::E(graph1)$weight/max(abs(igraph::E(graph1)$weight)))^2,
                                                               vcount = ncol(x$network))
    
    
    set.seed(1234)
    plot.args$net <- network1
    plot.args$node.color <- "Communities"
    plot.args$node.alpha <- plot.args$alpha
    plot.args$node.shape <- plot.args$shape
    node.size <- plot.args$node.size
    plot.args$node.size <- 0
    plot.args$color.palette <- NULL
    plot.args$palette <- NULL
    plot.args$edge.color <- "color"
    plot.args$edge.lty <- "line"
    plot.args$edge.size <- "ScaledWeights"
    
    lower <- abs(x$network[lower.tri(x$network)])
    non.zero <- sqrt(lower[lower != 0])
    
    plot.args$edge.alpha <- non.zero
    plot.args$mode <- layout.spring
    plot.args$label <- colnames(x$network)
    plot.args$node.label <- rep("", ncol(x$network))
    if(plot.args$label.size == "max_size/2"){plot.args$label.size <- plot.args$size/2}
    if(plot.args$edge.label.size == "max_size/2"){plot.args$edge.label.size <- plot.args$size/2}
    
    palette <- color_palette_EGA(color.palette, as.numeric(factor(x$wc)))
    palette <- ifelse(is.na(palette), "white", palette)
    
    ega.plot <- suppressWarnings(
      suppressMessages(
        do.call(GGally::ggnet2, plot.args) +
          ggplot2::theme(legend.title = ggplot2::element_blank())
      )
    )
    
    name <- colnames(x$network)
    
    name.split <- lapply(name, function(x){
      unlist(strsplit(x, split = " "))
    })
    
    name <- unlist(
      lapply(name.split, function(x){
        
        len <- length(x)
        
        if(len > 1){
          
          add.line <- round(len / 2)
          
          paste(
            paste(x[1:add.line], collapse = " "),
            paste(x[(add.line+1):length(x)], collapse = " "),
            sep = "\n"
          )
          
        }else{x}
        
      })
    )
    
    # Border color
    if(all(color.palette == "grayscale" |
           color.palette == "greyscale" |
           color.palette == "colorblind")){
      border.color <- ifelse(palette == "white", "white", "black")
    }else{border.color <- palette}
    
    # Custom nodes: transparent insides and dark borders
    ega.plot <- ega.plot +
      ggplot2::geom_point(ggplot2::aes(color = "color"), size = node.size,
                          color = border.color,
                          shape = 1, stroke = 1.5, alpha = .8) +
      ggplot2::geom_point(ggplot2::aes(color = "color"), size = node.size + .5,
                          color = palette,
                          shape = 19, alpha = plot.args$alpha) +
      ggplot2::geom_text(ggplot2::aes(label = name), color = "black", size = plot.args$label.size) +
      ggplot2::guides(
        color = ggplot2::guide_legend(override.aes = list(
          color = unique(palette),
          size = node.size,
          alpha = as.numeric(names(which.max(table(plot.args$alpha)))),
          stroke = 1.5
        ))
      )

    set.seed(NULL)

    # Return to object.list
    ega.plots[[i]] <- ega.plot

    # Reset plot.args
    plot.args <- original.plot.args

  }

  return(ega.plots)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# MULTI-FUNCTION SUB-ROUTINES ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @noRd
#'
# Zero-order correlations to partial correlations
# Updated 08.10.2022
cor2pcor <- function(corr)
{
  # Obtain inverse covariance
  inv <- solve(corr)
  
  # Convert to partial correlations
  pcor <- -cov2cor(inv)
  
  # Set diagonal to zero
  diag(pcor) <- 0
  
  # Return
  return(pcor)
  
}

#' @noRd
#'
# Partial correlations to zero-order correlations
# Updated 08.10.2022
pcor2cor <- function(pcor)
{
  
  # Obtain inverse covariance
  inv <- pcor2inv(pcor)
  
  # Convert to correlations
  corr <- cov2cor(inv)
  
  # Return
  return(corr)
  
}

#' @noRd
#'
# Partial correlations to inverse covariance matrix
# Updated 08.10.2022
pcor2inv <- function(pcor)
{
  
  # Set diagonal to negative 1
  diag(pcor) <- -1
  
  # Obtain inverse covariance
  inv <- solve(-pcor)
  
  # Return
  return(inv)
  
}

#' @noRd
#'
# Custom progress bar for timing
# Updated 02.08.2022
custom_progress <- function (
    i, max, start_time,
    caps = "|", progress = "+"
) {
  
  # Calculate percent complete
  percent <- i / max * 100
  
  # Check for calculating
  if(is.character(start_time)){
    
    # Update progress
    cat(
      sprintf(
        paste0("\r", caps, "%-49s", caps, " %s"),
        paste(rep(progress, percent / 2), collapse = ""),
        # Add percentage
        timing <- paste0(
          round(percent), "% ~",
          start_time
        )
      )
    )
    
  }else{
    
    # Obtain end time
    end_time <- Sys.time()
    
    # Time difference
    time_difference <- as.numeric(
      difftime(
        time1 = end_time,
        time2 = start_time,
        units = "sec"
      )
    )
    
    # Obtain time to finish based on remaining computations
    time_multiple <- (max - i) / i
    
    # Multiple time difference by remaining computations
    seconds <- time_multiple * time_difference
    if(i == max){
      seconds <- time_difference
    }
    
    # Obtain remaining minutes
    minutes <- as.numeric(seconds / 60)
    
    # Remaining seconds
    seconds <- round(60 * minutes %% 1)
    
    # Floor remaining minutes
    minutes <- floor(minutes)
    
    # Set timing with seconds
    timing <- paste0(
      formatC(
        seconds, digits = 1,
        flag = "0", format = "d"
      ), "s"
    )
    
    # Set timing with minutes
    if(minutes != 0){
      timing <- paste(
        paste0(minutes, "m"),
        timing
      )
    }
    
    # If max, then end
    if(i == max){
      
      # Add percentage
      timing <- paste0(
        round(percent), "% elapsed=",
        timing
      )
      
      # Update progress
      cat(
        sprintf(
          paste0("\r", caps, "%-49s", caps, " %s"),
          paste(rep(progress, percent / 2), collapse = ""),
          timing
        )
      )
      
      cat("\n")
      
    }else{
      
      # Add percentage
      timing <- paste0(
        round(percent), "% ~",
        timing
      )
      
      # Update progress
      cat(
        sprintf(
          paste0("\r", caps, "%-49s", caps, " %s        "),
          paste(rep(progress, percent / 2), collapse = ""),
          timing
        )
      )
      
    }
    
  }
  
}

#' Proportion table
#'
#' @param boot.mat Matrix.
#' A matrix of bootstrapped memberships
#'
#' @return Matrix of proportions based on dimensions
#'
#' @noRd
#'
# Proportion Table
# Updated 26.02.2021
proportion.table <- function (boot.mat)
{
  # Get maximum number of dimensions
  max.dim <- max(boot.mat, na.rm = TRUE)

  # Set up table
  tab <- matrix(0, nrow = nrow(boot.mat), ncol = max.dim)
  colnames(tab) <- 1:max.dim # assign column names

  # Check if there are row names
  if(!is.null(row.names(boot.mat))){
    row.names(tab) <- row.names(boot.mat)
  }

  # Loop through maximum dimensions
  for(i in 1:max.dim)
  {tab[,i] <- apply(boot.mat, 1, function(x){mean(x == i, na.rm = TRUE)})}

  return(tab)
}


#' @noRd
# Normalize DNN weights
# Updated 24.03.2021
min.max <- function(vec)
{
  exp.min <- exp(0) / exp(1)
  exp.max <- exp(1) / exp(0)

  return((vec - exp.min) / (exp.max - exp.min))
}

#' @noRd
# Custom range min-max
# Updated 26.05.2021
custom.min.max <- function(vec, ran)
{
  a <- ran[1]
  b <- ran[2]

  return((b - a) * ((vec - min(vec)) / (max(vec) - min (vec))) + a)
}

#' @noRd
# New NA function
# Updated 24.10.2021
is.NA <- function(x){

  # Regular check
  reg.na <- is.na(x)

  # Character check
  char.na <- x == "NA"

  # Get any NA
  return(as.logical(reg.na + char.na))

}

#' @noRd
# Function to create long results from list
# Updated 23.12.2021
long_results <- function(results_list){

  # Create long results
  rows <- unlist(lapply(results_list, nrow))
  end <- cumsum(rows)
  start <- (end + 1) - rows

  # Initialize matrix
  res_long <- matrix(
    ncol = ncol(results_list[[1]]),
    nrow = max(end)
  )
  colnames(res_long) <- colnames(results_list[[1]])

  # Loop through to populate
  for(i in seq_along(results_list)){
    res_long[start[i]:end[i],] <- as.matrix(results_list[[i]])
  }

  return(res_long)

}

#' @noRd
# Function to obtain arguments
# Updated 27.07.2022
obtain.arguments <- function(FUN, FUN.args)
{
  
  # Obtain formal arguments
  FUN.formals <- formals(FUN)
  
  # Check for input arguments
  if(length(FUN.args) != 0){
    
    ## Check for matching arguments
    if(any(names(FUN.args) %in% names(FUN.formals))){
      
      replace.args <- FUN.args[na.omit(match(names(FUN.formals), names(FUN.args)))]
      
      FUN.formals[names(replace.args)] <- replace.args
    }
    
  }
  
  # Remove ellipses
  if("..." %in% names(FUN.formals)){
    FUN.formals[which(names(FUN.formals) == "...")] <- NULL
  }
  
  # Return agrumnets
  return(FUN.formals)
  
}

#%%%%%%%%%%%%%%%%%%%%%%%%
# DATA GENERATION ----
#%%%%%%%%%%%%%%%%%%%%%%%%

#' @noRd
# Expand correlation matrix (unidimensional EGA)
# For EGA
# Updated 04.02.2021
expand.corr <- function(corr)
{
  # Number of variables
  nvar <- 4

  # Create additional correlations
  add.corr <- matrix(.50, nrow = nvar, ncol = nvar)
  diag(add.corr) <- 1

  # Initialize new correlation matrix
  new.corr <- matrix(0, nrow = nvar + nrow(corr), ncol = nvar + ncol(corr))

  # Input correlations into new correlation matrix
  new.corr[1:nvar, 1:nvar] <- add.corr
  new.corr[(nvar+1):nrow(new.corr), (nvar+1):ncol(new.corr)] <- corr

  # Add names
  colnames(new.corr) <- c(paste("SIM",1:nvar, sep = ""), colnames(corr))
  row.names(new.corr) <- colnames(new.corr)

  return(new.corr)
}

#' A sub-routine to simulate data for \code{\link[EGAnet]{EGA}}
#'
#' @param data Data for number of cases
#'
#' @param nvar Number of variables
#'
#' @param nfact Number of factors
#'
#' @param load Magnitude of loadings
#'
#' @return Simulated data
#'
#' @importFrom utils data
#'
#' @author Loli Dolores Nieto and Luis Eduardo Garrido
#'
#' @noRd
#'
# Simulate data function
# Updated 04.12.2020
sim.func <- function(data, nvar, nfact, load, n = 500)
{

  # Check for unidimensional structure
  ## Set up data simulation
  if(missing(data)){
    data <- matrix(nrow = n, ncol = 0)
  }else{
    n <- nrow(data)
  }
  corf <- 0
  J <- nvar*nfact
  sdcross = 0

  ## GENERATE SAMPLE DATA MATRIX
  check.eig <- TRUE
  check.com <- TRUE

  while(check.eig == TRUE|check.com == TRUE)
  {
    SATF = matrix(0, J, nfact)

    for(j in 1:nfact)
    {
      SATF[(j*nvar-nvar+1):(j*nvar),j]<-runif(nvar, load-.10, load+.10)

      if(nfact>1)
      {
        CROSS.L <- apply(as.matrix(SATF[(j*nvar-nvar+1+2):(j*nvar),-c(j)]), 2, function(x) rnorm((nvar-2), 0, sdcross))

        SATF[(j*nvar-nvar+1+2):(j*nvar),-c(j)] <- CROSS.L
      }
    }

    #SATF # Population factor loading matrix with cross-loadings and marker items

    FCOR      = matrix(corf, nfact, nfact); diag(FCOR)<-1 ## Factor correlation matrix
    R         = SATF%*%FCOR%*%t(SATF)                          ## Rr
    check.com = any(diag(R) > .90)                                  ## Check communalities values
    diag(R)   = 1                                                                    ## Insert ones in the diagonal of Rr
    #R                                                                                       ## Rp
    check.eig = any(eigen(R)$values <= 0)                      ## Check eigenvalues
  }

  U = chol(R)                                                                       ## Cholesky decomposition of Rp
  Z = MASS_mvrnorm(n, mu = rep(0, J), Sigma = diag(J))                                  ## Obtain sample matrix of continuous variables
  X = Z%*%U
  colnames(X) <- paste0("X", 1:ncol(X))

  data.sim <- cbind(X, data)

  return(data.sim)
}

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## DATA CATEGORIZATION FUNCTION
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Script: Function to categorize continuous data
## % DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## First, likely values for the skewness of the variables must be specified
## according to a set of values ranging from -2 to 2 in increments of 0.5.
## This is done with the skew.values argument. Second, for the selected
## number of response categories to simulate, an object with the thresholds
## to produce different skewness is generated. Then, the skewness for each
## variable is randomly assigned and each variable is categorized according
## to the thresholds previously defined. The thresholds used to categorize
## the data set are those defined by Garrido, Abad, & Ponsoda (2011, 2013).
## %% ARGUMENTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## data:        data set of continuous variables
## ncat:        number of response categories
## skew.values: a vector with several of the following values: -2, -1.5,
##              -1, -0.5, 0, 0.5, 1, 1.5, or 2. It can also be a positive
##              integer with a single value.
## %% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## data.cat:  categorized data
## SIM.SKEW:  simulated skewness for each variable
## TRUE.SKEW: real generated skewness for each variable
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Updated 04.12.2020
# Authors: Maria Dolores Nieto and Luis Eduardo Garrido
categorize<-function(data, ncat, skew.values){

  data.cat<-matrix(0, nrow = nrow(data), ncol = ncol(data))
  SIM.SKEW  = rep(0, ncol(data))
  TRUE.SKEW = rep(0, ncol(data))

  if(ncat==2){# SKEW = -2     -1.5     -1     -0.5  0   0.5     1     1.5     2
    thres<-matrix(c(-1.0518,-0.8416,-0.5936,-0.3088,0,0.3088,0.5936,0.8416,1.0518), ncol = 9)
    colnames(thres)<-seq(-2,2,0.5)
    SKEW.TABLE<-skew.values

    for(j in 1:ncol(data)){

      SKEW<-sample(SKEW.TABLE, 1)
      COLUMN<-which(colnames(thres)==SKEW)
      data.cat[,j][data[,j] < thres[1,COLUMN]]<-1
      data.cat[,j][data[,j] >=thres[1,COLUMN]]<-2
      SIM.SKEW[j]  = SKEW
      TRUE.SKEW[j] = psych::describe(data.cat[,j])$skew
    }
  }
  if(ncat==3){
    thres.2  <-c(0.8518,  1.3754)  ## SKEW = 2
    thres.1.5<-c(0.6131,  1.1969)  ## SKEW = 1.5
    thres.1  <-c(0.3195,  0.9921)  ## SKEW = 1.0
    thres.0.5<-c(-0.0236, 0.7256)  ## SKEW = 0.5
    thres.0  <-c(-1.0000, 1.0000)  ## SKEW = 0
    thres    <-matrix(c(sort(-thres.2), sort(-thres.1.5), sort(-thres.1), sort(-thres.0.5),
                        (thres.0),
                        thres.0.5, thres.1, thres.1.5, thres.2),
                      ncol = 9)
    colnames(thres)<-seq(-2,2,0.5)
    SKEW.TABLE<-skew.values

    for(j in 1:ncol(data)){
      SKEW<-sample(SKEW.TABLE, 1)
      COLUMN<-which(colnames(thres)==SKEW)
      data.cat[,j][ data[,j]< thres[1,COLUMN]]<-1
      data.cat[,j][(data[,j]>=thres[1,COLUMN])&(data[,j]<thres[2,COLUMN])]<-2
      data.cat[,j][ data[,j]>=thres[2,COLUMN]]<-3
      SIM.SKEW[j]  = SKEW
      TRUE.SKEW[j] = psych::describe(data.cat[,j])$skew
    }
  }
  if(ncat==4){
    thres.2  <-c(0.7515,  1.1341, 1.5980) ## SKEW = 2
    thres.1.5<-c(0.4945,  0.9299, 1.4359) ## SKEW = 1.5
    thres.1  <-c(0.1678,  0.6873, 1.2513) ## SKEW = 1.0
    thres.0.5<-c(-0.2057, 0.3706, 0.9809) ## SKEW = 0.5
    thres.0  <-c(-1.5000, 0.0000, 1.5000) ## SKEW = 0
    thres    <-matrix(c(sort(-thres.2), sort(-thres.1.5), sort(-thres.1), sort(-thres.0.5),
                        (thres.0),
                        thres.0.5, thres.1, thres.1.5, thres.2),
                      ncol = 9)
    colnames(thres)<-seq(-2,2,0.5)
    SKEW.TABLE<-skew.values

    for(j in 1:ncol(data)){
      SKEW<-sample(SKEW.TABLE, 1)
      COLUMN<-which(colnames(thres)==SKEW)
      data.cat[,j][ data[,j]< thres[1,COLUMN]]<-1
      data.cat[,j][(data[,j]>=thres[1,COLUMN])&(data[,j]<thres[2,COLUMN])]<-2
      data.cat[,j][(data[,j]>=thres[2,COLUMN])&(data[,j]<thres[3,COLUMN])]<-3
      data.cat[,j][ data[,j]>=thres[3,COLUMN]]<-4
      SIM.SKEW[j]  = SKEW
      TRUE.SKEW[j] = psych::describe(data.cat[,j])$skew
    }
  }
  if(ncat==5){
    thres.2  <-c(0.6792,  1.0043,  1.3441, 1.7703) ## SKEW = 2
    thres.1.5<-c(0.4071,  0.7827,  1.1596, 1.6186) ## SKEW = 1.5
    thres.1  <-c(0.0502,  0.5117,  0.9432, 1.4462) ## SKEW = 1.0
    thres.0.5<-c(-0.3414, 0.1642,  0.6257, 1.1645) ## SKEW = 0.5
    thres.0  <-c(-1.8000, -0.6000, 0.6000, 1.8000) ## SKEW = 0
    thres    <-matrix(c(sort(-thres.2), sort(-thres.1.5), sort(-thres.1), sort(-thres.0.5),
                        (thres.0),
                        thres.0.5, thres.1, thres.1.5, thres.2),
                      ncol = 9)
    colnames(thres)<-seq(-2,2,0.5)
    SKEW.TABLE<-skew.values

    for(j in 1:ncol(data)){
      SKEW<-sample(SKEW.TABLE, 1)
      COLUMN<-which(colnames(thres)==SKEW)
      data.cat[,j][ data[,j]< thres[1,COLUMN]]<-1
      data.cat[,j][(data[,j]>=thres[1,COLUMN])&(data[,j]<thres[2,COLUMN])]<-2
      data.cat[,j][(data[,j]>=thres[2,COLUMN])&(data[,j]<thres[3,COLUMN])]<-3
      data.cat[,j][(data[,j]>=thres[3,COLUMN])&(data[,j]<thres[4,COLUMN])]<-4
      data.cat[,j][ data[,j]>=thres[4,COLUMN]]<-5
      SIM.SKEW[j]  = SKEW
      TRUE.SKEW[j] = psych::describe(data.cat[,j])$skew
    }
  }

  RESULTS<-list(data.cat  = data.cat,
                SIM.SKEW  = SIM.SKEW,
                TRUE.SKEW = TRUE.SKEW)

  return(RESULTS)
}


#%%%%%%%%%%%%%
# GNN LCT ----
#%%%%%%%%%%%%%

#' @noRd
# Put data into pytorch format
# Updated 31.12.2021
torch_format <- function(data, ...)
{
  # Get arguments for EGAnet
  args <- list(...)
  args$data <- data
  args$plot.EGA <- FALSE

  # Estimate EGA
  ega <- suppressWarnings(
    suppressMessages(
      do.call(
        EGA, args
      )
    )
  )

  # Estimate graph
  graph <- abs(ega$network)

  # Round values to 5 decimal places
  graph <- round(graph, 5)

  # Obtain node attributes

  # Convert to igraph
  g <- convert2igraph(graph)

  # Obtain largest component
  components <- igraph::components(g)$membership == 1
  include <- components & !is.na(ega$wc)

  # Re-obtain igraph
  g <- convert2igraph(abs(graph[include, include]))

  # Check if dimensions = 1
  if(ega$n.dim == 1){

    # Loadings
    ## Network
    network_dom <- abs(as.vector(as.matrix(net.loads(ega)$std))[include]) # Absolute loadings
    network_cross <- rep(0, length(network_dom))

    ## Factor
    factor_loads <- suppressMessages(
      suppressWarnings(
        psych::fa(
          ega$correlation,
          nfactors = 1,
          n.obs = nrow(data)
        )$loadings[,1]
      )
    )
    factor_dom <- abs(as.vector(factor_loads)[include]) # Absolute loadings
    factor_cross <- rep(0, length(factor_dom))

  }else{

    # Loadings
    ## Network
    network_loads <- abs(net.loads(ega)$std) # Absolute loadings
    network_loads <- network_loads[names(ega$wc[include]),]
    network_loads <- network_loads[,order(colnames(network_loads))]
    network_dom <- unlist(lapply(1:nrow(network_loads), function(i){
      network_loads[i, ega$wc[include][i]]
    }))
    network_cross <- unlist(lapply(1:nrow(network_loads), function(i){
      sum(network_loads[i, -ega$wc[include][i]])
    }))

    ## Factor
    factor_loads <- suppressMessages(
      suppressWarnings(
        psych::fa(
          ega$correlation,
          nfactors = ncol(network_loads),
          n.obs = nrow(data)
        )$loadings[,1:ncol(network_loads)]
      )
    )
    factor_loads <- abs(factor_loads[names(ega$wc[include]),]) # Absolute loadings
    factor_dom <- unlist(lapply(1:nrow(factor_loads), function(i){
      factor_loads[i,which.max(factor_loads[i,])]
    }))
    factor_cross <- unlist(lapply(1:nrow(factor_loads), function(i){
      sum(factor_loads[i,-which.max(factor_loads[i,])])
    }))

  }

  # Network measures
  q <- max(igraph::cluster_louvain(g)$modularity, na.rm = TRUE)

  node_attributes <- round(cbind(
    q,
    factor_cross,
    network_cross
  ), 5)

  ## Graph attributes
  # aspl <- pathlengths(graph)$ASPL
  # cc <- clustcoeff(graph)$CC
  # q <- max(igraph::cluster_louvain(convert2igraph(abs(graph)))$modularity)
  # graph_attributes <- round(c(aspl, cc, q), 5)
  # names(graph_attributes) <- c("aspl", "cc", "q")

  # Remove non-dimensional nodes
  graph <- graph[include, include]

  # Make graph sparse
  graph <- sparsify(graph)

  # Remove zero weights
  graph[,"weight"] <- ifelse(graph[,"weight"] == 0, NA, graph[,"weight"])
  graph <- na.omit(graph)
  attr(graph, "na.action") <- NULL

  # Obtain remaining nodes
  node_label <- sort(unique(as.vector(graph[,c("from", "to")])))

  # Set up re-numbering for nodes
  node_num <- seq_along(node_label) # set to zero for Python
  names(node_num) <- node_label

  # Replace node numbers to be sequential
  graph[,"from"] <- node_num[as.character(graph[,"from"])]
  graph[,"to"] <- node_num[as.character(graph[,"to"])]

  # Remove node_attributes for disconnected nodes
  node_attributes <- node_attributes[node_label,]

  # Set graph indicator
  graph_indicator <- rep(1, length(node_label))

  # Return results
  results <- list()
  results$A <- graph[,c("from", "to")]
  results$edge_weight <- graph[,"weight"]
  results$node_attributes <- node_attributes
  # results$graph_attributes <- graph_attributes
  results$graph_indicator <- graph_indicator

  return(results)
}

#' @noRd
# Sparse matrix
# Updated 31.12.2021
sparsify <- function(network){

  # Number of nodes
  nodes <- ncol(network)

  # Sparse matrix (with graph label)
  from <- rep(1:nodes, each = nodes)
  to <- rep(1:nodes, times = nodes)
  weight <- as.vector(network)
  sparse_matrix <- cbind(
    from, to, weight
  )

  return(sparse_matrix)

}

#' @noRd
# Normalize function
# Updated 31.12.2021
normalize <- function(values, min_desired, max_desired){

  min_actual <- min(values, na.rm = TRUE)
  num <- values - min_actual
  denom <- max(values, na.rm = TRUE) - min_actual
  mult  <- max_desired - min_desired
  add <- min_desired

  return((num / denom) * mult + add)

}

#%%%%%%%%%
# LCT ----
#%%%%%%%%%

#' @noRd
# Dynamic organization
# Updated 26.05.2021
dyn.org <- function(data, gen)
{
  for(i in 1:ncol(data)){

    gen[,i] <- custom.min.max(
      sort(gen[,i])[rank(data[,i], ties.method = "first")],
      range(data[,i])
    )

  }

  return(gen)
}

#' A sub-routine to compute the deep learning neural network model
#' weights for \code{\link[EGAnet]{LCT}}
#'
#' @param loads Matrix of loadings
#'
#' @param weights Weights from specific model (see \code{\link[EGAnet]{dnn.weights}})
#'
#' @return A prediction or probability of the specified model
#'
#' @noRd
#'
# DNN weights function
# Updated 09.08.2020
dnn.model.weights <- function (loads, weights)
{
  wb <- seq(1, length(weights), 3)

  for(i in wb[-length(wb)])
  {
    if(i == 1)
    {
      input <- as.vector(t(weights[[i]]) %*% as.matrix(loads)) + weights[[(i+1)]]
      input <- ifelse(input > 0, input, input * weights[[i+2]])
      layer <- input
    }else{
      layer <- as.vector(t(weights[[i]]) %*% as.matrix(layer)) + weights[[(i+1)]]
      layer <- ifelse(layer > 0, layer, layer * weights[[i+2]])
    }
  }

  # Output
  output <- as.vector(t(weights[[wb[length(wb)]]]) %*% as.matrix(layer)) + weights[[length(weights)]]

  # Sigmoid activation function
  sigmoid <- function(x)
  {exp(x) / (exp(x) + 1)}

  # Prediction
  prediction <- sigmoid(output)

  return(prediction)
}

#' A sub-routine to predict the model for \code{\link[EGAnet]{LCT}}
#'
#' @param loads Matrix of loadings
#'
#' @return The model prediction
#'
#' @importFrom utils data
#'
#' @noRd
#'
# DNN prediction function
# Updated 30.03.2021
dnn.predict <- function (loads)
{
  # Load deep learning neural network weights
  dnn.weights <- get(data("dnn.weights", envir = environment()))

  # Compute ratios
  ## Small ratio
  small.ratio <- exp(loads[1]) / exp(loads[6])
  ## Moderate ratio
  moderate.ratio <- exp(loads[2]) / exp(loads[7])
  ## Large ratio
  large.ratio <- exp(loads[3]) / exp(loads[8])
  ## Dominant ratio
  dominant.ratio <- exp(loads[4]) / exp(loads[9])
  ## Cross ratio
  cross.ratio <- exp(loads[5]) / exp(loads[10])

  small.ratio <- min.max(small.ratio)
  moderate.ratio <- min.max(moderate.ratio)
  large.ratio <- min.max(large.ratio)
  dominant.ratio <- min.max(dominant.ratio)
  cross.ratio <- min.max(cross.ratio)

  # Factor versus network model
  f_n <- vector("numeric", length = 3)

  # Check for low correlation factor versus network model
  f_n[1] <- dnn.model.weights(c(loads, dominant.ratio), dnn.weights$lf_n_weights)

  # Check for high correlation with variables greater than factors versus network model
  f_n[2] <- dnn.model.weights(c(loads, dominant.ratio), dnn.weights$hvgf_n_weights)

  # Check for high correlation factor versus network model
  f_n[3] <- dnn.model.weights(c(loads, dominant.ratio), dnn.weights$hvlf_n_weights)

  # Check for factor model
  ifelse(any(f_n >= .50), return(1), return(2))

}

#%%%%%%%%%
# UVA ----
#%%%%%%%%%

#' @noRd
# Redundancy Processing
# Updated 04.05.2022
redundancy.process <- function(data, cormat, n, model, method, type, sig, plot.redundancy, plot.args)
{

  if(method == "wto"){

    if(model == "glasso"){

      for(i in c(0.50, 0.25, 0))
      {
        net <- EBICglasso.qgraph(data = cormat, n = n, gamma = i)

        if(all(colSums(net)!=0))
        {break}
      }

    }else if(model == "tmfg"){

      net <- TMFG(cormat)$A

    }else{

      stop(paste(model, "does not exist as an option for the argument 'model'"))

    }

    tom <- wTO(net, sign = "sign")

  }else if(method == "pcor"){

    tom <- -cov2cor(solve(cormat))

  }else{tom <- cormat}


  # Ensure column names
  if(is.null(colnames(tom))){
    colnames(tom) <- paste("V", 1:ncol(tom), sep = "")
    row.names(tom) <- colnames(tom)
  }

  # Number of variables; diagonal zero; absolute values
  vars <- ncol(tom); diag(tom) <- 0; tom <- abs(tom)

  # Lower triangle
  lower <- tom[lower.tri(tom, diag = FALSE)]

  # Names lower triangle
  name1 <- colnames(tom); name2 <- name1

  # Initialize name matrix
  name.matrix <- tom

  # Generate name matrix
  for(i in 1:vars)
    for(j in 1:vars){
      name.matrix[i,j] <- paste(name1[j], name2[i], sep = "--")
    }

  # Get lower triangle of names
  names(lower) <- name.matrix[lower.tri(name.matrix, diag = FALSE)]

  # Obtain positive values only
  pos.vals <- na.omit(ifelse(lower == 0, NA, lower)); attr(pos.vals, "na.action") <- NULL

  # Get redundant pairings
  if(type == "threshold"){## Threshold values
    redund <- pos.vals[which(pos.vals >= sig)]
    aic <- NULL
    g.dist <- NULL
  }else{## Determine distribution
    
    stop('"alpha" and "adapt" have been removed from UVA. Please use "threshold"')

    # # Distributions, initialize AIC vector
    # distr <- c("norm", "gamma")
    # aic <- numeric(length(distr))
    # names(aic) <- c("normal", "gamma")
    # 
    # ## Obtain distribution
    # for(i in 1:length(distr)){
    #   capture.output(
    #     aic[i] <- fitdistrplus::fitdist(pos.vals, distr[i], method="mle")$aic
    #   )
    # }
    # 
    # ## Obtain parameters
    # g.dist <- suppressWarnings(
    #   fitdistrplus::fitdist(
    #     pos.vals, distr[which.min(aic)]
    #   )$estimate
    # )
    # 
    # # Estimate p-values
    # pval <- switch(names(aic)[which.min(aic)],
    # 
    #                normal = 1 - unlist(lapply(pos.vals, # positive values
    #                                           pnorm, # probability in normal distribution
    #                                           mean = g.dist["mean"], #mean of normal
    #                                           sd = g.dist["sd"]) #standard deviation of normal
    #                ),
    # 
    #                gamma = 1 - unlist(lapply(pos.vals, # positive values
    #                                          pgamma, # probability in gamma distribution
    #                                          shape = g.dist["shape"], # shape of gamma
    #                                          rate = g.dist["rate"]) # rate of gamma
    #                ),
    # )
    # 
    # # Check if using adaptive alpha
    # if(type == "adapt"){
    #   sig <- adapt.a("cor", alpha = sig, n = length(pos.vals), efxize = "medium")$adapt.a
    # }
    # 
    # # Get redundant pairings
    # redund <- pos.vals[which(pval <= sig)]

  }

  # Check for redundant pairings
  if(length(redund) == 0){
    message("No redundant variables identified.")
    res.list <- NA
  }else{

    # Create result matrix
    split.res <- unlist(strsplit(names(redund), split = "--"))
    res.mat <- t(simplify2array(sapply(names(redund), strsplit, split = "--")))

    # Initialize result list
    res.list <- list()

    # Initialize count
    count <- 0

    while(nrow(res.mat)!=0){
      # Increase count
      count <- count + 1

      # Get variable counts
      var.counts <- sort(table(split.res), decreasing = TRUE)

      if(!all(var.counts==1)){
        # Identify targets
        target <- which(res.mat == names(var.counts[1]), arr.ind = TRUE)[,"row"]

        # Insert values into list
        res.list[[names(var.counts[1])]] <- setdiff(unlist(strsplit(names(target),split="--")),names(var.counts[1]))

        # Remove rows from result matrix
        res.mat <- res.mat[-target,]

        # Remove variables from split result
        split.res <- as.vector(res.mat)

        # Force matrix
        if(is.vector(res.mat))
        {res.mat <- t(as.matrix(res.mat))}

      }else{
        for(i in 1:nrow(res.mat))
        {res.list[[res.mat[i,1]]] <- unname(res.mat[i,2])}

        res.mat <- res.mat[-c(1:nrow(res.mat)),]
      }
    }

  }

  # Revert tom to matrix
  tom <- as.matrix(tom)

  # Check for plot
  if(plot.redundancy){

    # Initialize plot matrix
    plot.mat <- matrix(0, nrow = nrow(tom), ncol = ncol(tom))
    colnames(plot.mat) <- colnames(tom)
    row.names(plot.mat) <- colnames(tom)

    for(i in 1:length(res.list)){
      plot.mat[names(res.list)[i],res.list[[i]]] <- tom[names(res.list)[i],res.list[[i]]]
      plot.mat[res.list[[i]],names(res.list)[i]] <- tom[res.list[[i]],names(res.list)[i]]
    }

    if(any(colSums(plot.mat) == 0)){

      rm.mat <- which(colSums(plot.mat) == 0)
      plot.mat <- plot.mat[-rm.mat, -rm.mat]

    }

    plot.args$title <- switch(method,
                              "wto" = "Weighted\nTopological\nOverlap",
                              "pcor" = "Partial\nCorrelation",
                              "cor" = "Zero-order\nCorrelation",
                              "irt" = "IRT\nCorrelated\nResiudals"
    )

    if(ncol(plot.mat) <= 2){
      warning("No plot was produced because there are only two redundant variables")
    }else{

      # Global suppress warnings (need a better workaround)
      warn <- options("warn")[[1]]
      options(warn = -1)
      suppressMessages(
        net.plot <- redund.plot(plot.mat, plot.args)
      )
      options(warn = warn)
    }

  }

  # Get redundancy descriptives
  desc <- redund.desc(pos.vals = pos.vals, method = method, type = type, sig = sig,
                      aic = aic, g.dist = g.dist)

  # Add p-values
  # if(type != "threshold"){
  #   desc$basic <- cbind(round(sig, 5), desc$basic)
  #   colnames(desc$basic)[1] <- "Sig"
  #   desc$centralTendency <- cbind(round(pval[row.names(desc$centralTendency)], 5),
  #                                 desc$centralTendency)
  #   colnames(desc$centralTendency)[1] <- "p-value"
  # }

  # Results list
  res <- list()
  res$redundant <- res.list
  res$data <- data
  res$correlation <- cormat
  res$weights <- tom
  if(exists("net")){res$network <- net}
  if(exists("net.plot")){res$plot <- net.plot}
  res$descriptives <- desc
  res$method <- method
  res$model <- model
  res$type <- type
  # if(type != "threshold"){res$distribution <- names(aic)[which.min(aic)]}

  class(res) <- "node.redundant"

  return(res)

}

#' @noRd
# Redundancy Naming
# Updated 13.12.2020
redund.names <- function(node.redundant.obj, key)
{
  # Check for node.redundant object class
  if(is(node.redundant.obj) != "node.redundant"){
    stop("A 'node.redundant' object must be used as input")
  }

  # Obtain and remove data from node redundant object
  data <- node.redundant.obj$data

  # Check that columns match key
  if(ncol(data) != length(as.vector(key))){
    stop("Number of columns in data does not match the length of 'key'")
  }

  # Names of node.redundant object
  nr.names <- names(node.redundant.obj$redundant)

  # Key names
  key.names <- colnames(data)

  # Key change
  key.chn <- key

  for(i in 1:length(nr.names)){

    # Target redundant node
    target.r <- match(names(node.redundant.obj$redundant)[i],key.names)

    # Replace item name with description
    names(node.redundant.obj$redundant)[i] <- as.character(key.chn[target.r])

    # Target other nodes
    target.o <- match(node.redundant.obj$redundant[[i]],key.names)

    # Replace item names with description
    node.redundant.obj$redundant[[i]] <- as.character(key.chn[target.o])

  }

  # Create key code
  names(key) <- colnames(data)
  node.redundant.obj$key <- key

  return(node.redundant.obj)
}

#' @noRd
# Redundancy Descriptives
# Updated 13.12.2020
redund.desc <- function(pos.vals, method, type, sig, aic, g.dist)
{
  # Initialize descriptives matrix
  desc <- matrix(0, nrow = 1, ncol = 9)

  # Row name
  row.names(desc) <- switch(method,
                            "wto" = "wTO",
                            "pcor"= "pcor",
                            "cor" = "cor"
  )

  colnames(desc) <- c("Mean", "SD", "Median", "MAD", "3*MAD", "6*MAD", "Minimum", "Maximum", "Critical Value")

  desc[,"Mean"] <- mean(pos.vals, na.rm = TRUE)
  desc[,"SD"] <- sd(pos.vals, na.rm = TRUE)
  desc[,"Median"] <- median(pos.vals, na.rm = TRUE)
  desc[,"MAD"] <- mad(pos.vals, constant = 1, na.rm = TRUE)
  desc[,"3*MAD"] <- mad(pos.vals, constant = 1, na.rm = TRUE) * 3
  desc[,"6*MAD"] <- mad(pos.vals, constant = 1, na.rm = TRUE) * 6
  desc[,"Minimum"] <- range(pos.vals, na.rm = TRUE)[1]
  desc[,"Maximum"] <- range(pos.vals, na.rm = TRUE)[2]

  # Critical value
  if(type == "threshold"){
    desc[,"Critical Value"] <- sig
  }else{

    desc[,"Critical Value"] <- switch(names(aic)[which.min(aic)],

                                      normal = qnorm(sig, #significance
                                                     mean = g.dist["mean"], #mean of normal
                                                     sd = g.dist["sd"], #sd of normal
                                                     lower.tail = FALSE),

                                      gamma = qgamma(sig, #significance
                                                     shape = g.dist["shape"], #shape of gamma
                                                     rate = g.dist["rate"], #rate of gamma
                                                     lower.tail = FALSE),

    )

  }

  # Organize positive values output
  ordered.pos <- sort(pos.vals, decreasing = TRUE)
  sd.from.mean <- (ordered.pos - mean(pos.vals, na.rm = TRUE)) / sd(ordered.pos, na.rm = TRUE)
  mad.from.median <- (ordered.pos - median(pos.vals, na.rm = TRUE)) / mad(ordered.pos, constant = 1, na.rm = TRUE)
  pos.output <- round(cbind(ordered.pos, sd.from.mean, mad.from.median), 3)

  colnames(pos.output)[1] <- switch(method,
                                    "wto" = "wTO",
                                    "pcor"= "pcor",
                                    "cor" = "cor",
                                    "irt" = "IRT"
  )

  colnames(pos.output)[2:3] <- c("SD from Mean", "MAD from Median")

  res.desc <- list()
  res.desc$basic <- round(desc, 3)
  res.desc$centralTendency <- pos.output

  return(res.desc)

}

#' @noRd
# Redundancy Plot
# Updated 13.12.2020
redund.plot <- function(plot.matrix, plot.args, plot.reduce = FALSE)
{
  # Convert to plot.mat
  plot.mat <- plot.matrix

  # weighted  network
  network1 <- network::network(plot.mat,
                               ignore.eval = FALSE,
                               names.eval = "weights",
                               directed = FALSE)

  if(isTRUE(plot.reduce)){
    wc <- c("Target", rep("Possible", ncol(plot.mat)-1))
    network::set.vertex.attribute(network1, attrname= "Communities", value = wc)
  }else{
    wc <- rep(plot.args$title, ncol(plot.mat))
    network::set.vertex.attribute(network1, attrname= "Communities", value = wc)
  }

  network::set.vertex.attribute(network1, attrname= "Names", value = network::network.vertex.names(network1))
  network::set.edge.attribute(network1, "color", ifelse(network::get.edge.value(network1, "weights") > 0, plot.args$edge.color[1], plot.args$edge.color[2]))
  network::set.edge.attribute(network1, "line", ifelse(network::get.edge.value(network1, "weights") > 0, plot.args$edge.lty[1], plot.args$edge.lty[2]))
  network::set.edge.value(network1,attrname="AbsWeights",value=abs(plot.mat))
  network::set.edge.value(network1,attrname="ScaledWeights",
                          value=matrix(rescale.edges(plot.mat, 5),
                                       nrow = nrow(plot.mat),
                                       ncol = ncol(plot.mat)))

  # Layout "Spring"
  graph1 <- convert2igraph(plot.mat)
  edge.list <- igraph::as_edgelist(graph1)
  layout.spring <- qgraph::qgraph.layout.fruchtermanreingold(edgelist = edge.list,
                                                             weights =
                                                               abs(igraph::E(graph1)$weight/max(abs(igraph::E(graph1)$weight)))^2,
                                                             vcount = ncol(plot.mat))


  lower <- abs(plot.mat[lower.tri(plot.mat)])
  non.zero <- sqrt(lower[lower != 0])

  set.seed(1234)
  plot.args$net <- network1
  plot.args$node.color <- "Communities"
  plot.args$node.alpha <- plot.args$alpha
  plot.args$node.shape <- plot.args$shape
  node.size <- plot.args$node.size
  plot.args$node.size <- 0
  color.palette <- plot.args$color.palette
  plot.args$color.palette <- NULL
  plot.args$palette <- NULL
  plot.args$edge.lty <- "line"
  plot.args$edge.color <- "color"
  plot.args$edge.size <- "ScaledWeights"

  lower <- abs(plot.mat[lower.tri(plot.mat)])
  non.zero <- sqrt(lower[lower != 0])

  plot.args$edge.alpha <- non.zero
  plot.args$mode <- layout.spring
  plot.args$label <- colnames(plot.mat)
  plot.args$node.label <- rep("", ncol(plot.mat))
  if(plot.args$label.size == "max_size/2"){plot.args$label.size <- plot.args$size/2}
  if(plot.args$edge.label.size == "max_size/2"){plot.args$edge.label.size <- plot.args$size/2}

  redund.net <- suppressMessages(
    do.call(GGally::ggnet2, plot.args) +
      ggplot2::theme(legend.title = ggplot2::element_blank()) +
      ggplot2::scale_color_manual(values = color_palette_EGA(color.palette, na.omit(as.numeric(factor(wc)))),
                                  breaks = sort(as.numeric(factor(wc)))) +
      ggplot2::guides(
        color = ggplot2::guide_legend(override.aes = list(
          size = node.size,
          alpha = plot.args$alpha,
          stroke = 1.5
        ))
      )
  )

  redund.net <- redund.net + ggplot2::annotate("text", x = -Inf, y = Inf,
                                               hjust = 0, vjust = 1,
                                               label = plot.args$title, size = 5.5)

  set.seed(NULL)

  name <- colnames(plot.mat)

  # Custom nodes: transparent insides and dark borders
  redund.net <- redund.net +
    ggplot2::geom_point(ggplot2::aes(color = "color"), size = node.size,
                        color = color_palette_EGA(color.palette, na.omit(as.numeric(factor(wc))), sorted = FALSE),
                        shape = 1, stroke = 1.5, alpha = .8) +
    ggplot2::geom_point(ggplot2::aes(color = "color"), size = node.size + .5,
                        color = color_palette_EGA(color.palette, na.omit(as.numeric(factor(wc))), sorted = FALSE),
                        shape = 19, alpha = plot.args$alpha) +
    ggplot2::geom_text(ggplot2::aes(label = name), color = "black", size = plot.args$label.size)

  return(redund.net)
}

#' @importFrom graphics text
#' @noRd
# Redundancy Reduction
# Updated 04.05.2022
redund.reduce <- function(node.redundant.obj, reduce.method, plot.args, lavaan.args, corr)
{
  # Check for node.redundant object class
  if(is(node.redundant.obj) != "node.redundant"){
    stop("A 'node.redundant' object must be used as input")
  }

  # Redundant list
  redund <- node.redundant.obj$redundant

  # Copied data
  new.data <- node.redundant.obj$data

  # Weights
  if("network" %in% names(node.redundant.obj)){
    weights <- as.matrix(node.redundant.obj$network)
  }else if("correlation" %in% names(node.redundant.obj)){
    weights <- as.matrix(node.redundant.obj$correlation)
  }else{weights <- as.matrix(node.redundant.obj$weights)}

  # Track merged items
  merged <- list()

  # Track changed names
  name.chn <- vector("character")

  # Initialize count
  count <- 0

  # Get key
  if("key" %in% names(node.redundant.obj)){
    key <- node.redundant.obj$key
    names(key) <- names(node.redundant.obj$key)
  }else{
    key <- colnames(node.redundant.obj$data)
    names(key) <- key
  }

  # Line break
  linebreak()

  # Previous state
  prev.state <- list(redund)
  prev.state.data <- list(new.data)

  # Loop through named node redundant list
  while(length(redund) != 0){

    # Targeting redundancy
    target.item <- names(redund)[1]

    # Potential redundancies
    pot <- redund[[1]]

    # Construct potential redundancies and menu
    poss <- redundancy.menu(redund = redund, reduce.method = reduce.method,
                            pot = pot, target.item = target.item,
                            weights = weights, plot.args = plot.args,
                            key = key, node.redundant.obj = node.redundant.obj)

    # Get input
    input <- input.check(poss, type = "redund")

    # Check if going back is necessary
    if(any(input == "b")){

      # Let user know they can't go back
      if(is.null(unlist(prev.state[length(prev.state) - 1], recursive = FALSE))){
        message("\nCannot go back. This is the start of the redundancy list")
      }else{

        # Renew redund list
        redund <- unlist(prev.state[length(prev.state) - 1], recursive = FALSE)

        # Remove previous state
        prev.state <- prev.state[-length(prev.state)]

        # Renew new data
        new.data <- as.data.frame(prev.state.data[length(prev.state.data) - 1])

        # Remove previous state data
        prev.state.data <- prev.state.data[-length(prev.state.data)]

        # Reduce count
        count <- count - 1

        # Remove merge
        merged <- merged[-length(merged)]

      }

    }else if(all(input != "0")){

      # Convert to numeric
      re.items <- as.numeric(unlist(strsplit(unlist(strsplit(input, split = " ")), split = ",")))

      # Items to combine with target
      comb <- poss[re.items]

      # Index items
      idx <- names(key)[match(comb, key)]

      # Target index
      tar.idx <- names(key)[match(target.item, key)]

      # Update merged list
      count <- count + 1
      merged[[count]] <- c(key[tar.idx], key[idx])

      # Combine into target index
      if(reduce.method == "latent"){

        # Latent variable
        ## create model
        mod <- paste(paste("comb =~ ",sep=""), paste(colnames(new.data[,c(tar.idx, idx)]), collapse = " + "))

        # Replace arguments
        lavaan.args$model <- mod
        lavaan.args$data <- new.data
        
        ## Get default estimator
        categories <- apply(new.data[,c(tar.idx, idx)], 2, function(x){
          length(unique(x))
        })

        # Check categories
        if(sum(categories < 6) > 1){# Not all continuous
          lavaan.args$estimator <- "WLSMV"
          lavaan.args$missing <- "pairwise"
          lavaan.args$ordered <- TRUE
        }else{# All can be considered continuous
          lavaan.args$estimator <- "MLR"
          lavaan.args$missing <- "fiml"
        }

        ## get CFA function from lavaan
        FUN <- lavaan::cfa

        ## fit model
        fit <- suppressWarnings(
          suppressMessages(
            do.call(what = "FUN", args = as.list(lavaan.args))
          )
        )

        ## identify cases
        cases <- lavaan::inspect(fit, "case.idx")

        ## compute latent variable score
        latent <- as.numeric(lavaan::lavPredict(fit))

        ## check for missing cases and handle
        if(length(cases) != nrow(new.data))
        {
          new.vec <- as.vector(matrix(NA, nrow = nrow(new.data), ncol = 1))
          new.vec[cases] <- latent
        }else{new.vec <- latent}

        ## check for reverse scoring/labelling
        corrs <- as.matrix(cor(cbind(latent,new.data[,c(tar.idx, idx)]), use = "complete.obs")[1,-1])
        row.names(corrs) <- c(key[tar.idx], key[idx])
        colnames(corrs) <- "latent"

        if(any(sign(corrs)==-1))
        {
          message("Some variables are reverse coded (negative correlations with latent variable were found). Correlations with latent variable:")
          print(corrs)

          input2 <- "o"

          while(input2 != "y" && input2 != "n")
          {input2 <- readline("Reverse code for positive labelling (y/n): ")}

          if(input2 == "y")
          {new.vec <- -new.vec}
        }

        # input new vector
        new.data[,tar.idx] <- new.vec

        # Ask for new label
        lab <- readline(prompt = "New label for latent variable (no quotations): ")
        name.chn[count] <- lab
        col.idx <- match(tar.idx, colnames(new.data))
        colnames(new.data)[col.idx] <- lab

        message(paste("\nNew LATENT variable called '", lab,"' was created. Redundant variables were REMOVED", sep = ""))

      }else if(reduce.method == "remove" | reduce.method == "sum"){

        target.key <- c(tar.idx, idx)
        target.data <- new.data[,target.key]

        means <- round(colMeans(target.data, na.rm = TRUE), 2)
        sds <- round(apply(target.data, 2, sd, na.rm = TRUE), 2)
        ranges <- round(apply(target.data, 2, range, na.rm = TRUE), 2)

        if(ncol(target.data) > 2){

          # Corrected item-total correlations
          cor.corr <- round(item.total(target.data, corr), 2)

          ## Use information utility?
          categories <- apply(target.data, 2, function(x){
            length(unique(x))
          })

          ## Check categories
          if(all(categories <= 7)){

            ## Use
            #util <- round(info.util(target.data), 2)

            tab <- cbind(#util,
              cor.corr, means, sds, t(ranges))
            colnames(tab) <- c(#"Utility Gain",
              "Item-Total r", "Mean", "SD", "Low", "High")
          }else{

            tab <- cbind(cor.corr, means, sds, t(ranges))
            colnames(tab) <- c("Item-Total r", "Mean", "SD", "Low", "High")

          }

        }else{
          tab <- cbind(means, sds, t(ranges))
          colnames(tab) <- c("Mean", "SD", "Low", "High")
        }
        row.names(tab) <- c("0 (Target)", 1:length(comb))

        tab[,1:(ncol(tab) - 2)] <- matrix(sprintf("%.2f", tab[,1:(ncol(tab) - 2)]), nrow = nrow(tab), ncol = ncol(tab) - 2)

        if(!isSymmetric(new.data)){
          gridExtra::grid.arrange(gridExtra::tableGrob(tab))
        }

        # Input check
        new.input <- input.check(poss = c(target.item, comb), type = "remove")

        # All variables
        ind <- names(key[match(c(target.item, comb), key)])

        # All except selected variable
        idx <- ind[-(as.numeric(new.input)+1)]

        # Input selected
        merged[[count]] <- key[idx]

        # Name merged input
        name.chn[count] <- key[setdiff(ind, idx)]

        # Message user
        message(paste("\nKEPT '", key[ind[as.numeric(new.input) + 1]],"' and REMOVED all others", sep = ""))

      }

      if(reduce.method != "sum"){

        # Remove redundant variables from data
        rm.idx <- match(idx, colnames(new.data))

        if(isSymmetric(new.data)){
          new.data <- new.data[-rm.idx, -rm.idx]
        }else{
          new.data <- new.data[,-rm.idx]
        }

      }

      # Update previous state data
      prev.state.data[length(prev.state.data) + 1] <- list(new.data)

      # Remove target item
      redund[[1]] <- NULL

      # Make ind
      if(reduce.method == "latent" | reduce.method == "sum"){ind <- idx}

      # Remove variables within future options
      for(i in 1:length(ind)){

        ## Get elements from other redundancy lists
        elements <- sapply(redund, function(x){key[ind[i]] %in% x})

        ## If there are any
        if(any(elements)){

          ### Target elements
          target.elem <- which(elements)

          for(j in 1:length(target.elem)){

            #### Target list
            list.target <- redund[[target.elem[j]]]

            #### Remove from target list
            redund[[target.elem[j]]] <- redund[[target.elem[j]]][-which(key[ind[i]] == list.target)]

          }

        }

      }

      # Remove empty list objects
      rm.list <- which(
        unlist(lapply(redund, function(x){
          length(x)
        })) == 0
      )

      if(length(rm.list) != 0){
        redund <- redund[-rm.list]
      }

      # Remove object names
      if(any(key[ind] %in% names(redund))){

        ## Get target names
        name.targets <- which(key[ind] %in% names(redund))

        ## Loop through
        for(i in 1:length(name.targets)){

          # If there is only one item left, then remove from list
          if(length(redund[[key[ind][name.targets[i]]]]) == 1){
            redund[key[ind][name.targets[i]]] <- NULL
          }else{# Otherwise, replace name with first element and remove first element from list
            names(redund[key[ind]][name.targets[i]]) <- redund[[key[ind][name.targets[i]]]][1]
            redund[[key[ind][name.targets[i]]]][1] <- NA
            redund[[key[ind][name.targets[i]]]] <- na.omit(redund[[key[ind][name.targets[i]]]])
            redund[key[ind][name.targets[i]]] <- NULL
          }

        }

      }

      # Update previous state
      prev.state[length(prev.state) + 1] <- list(redund)

    }else{

      # Update previous state data
      prev.state.data[length(prev.state.data) + 1] <- list(new.data)

      redund[[1]] <- NULL

      # Update previous state
      prev.state[length(prev.state) + 1] <- list(redund)

    }

    if(!is.null(input)){
      linebreak()
      input <- NULL
    }

  }

  # Transform merged list to matrix
  if(length(merged) != 0){

    # Number of rows for matrix
    m.rows <- max(unlist(lapply(merged, length)))

    # Initialize merged matrix
    m.mat <- matrix("", nrow = m.rows, ncol = length(merged))

    # Input into merged matrix
    for(i in 1:length(merged))
    {
      diff <- m.rows - length(merged[[i]])

      m.mat[,i] <- c(merged[[i]], rep("", diff))
    }

    colnames(m.mat) <- name.chn

  }

  # Replace column names for item names not changed
  if(any(colnames(new.data) %in% names(key))){

    # Target names
    target.names <- which(colnames(new.data) %in% names(key))

    # new.data names
    new.data.names <- colnames(new.data)[target.names]

    # Insert into new data
    colnames(new.data)[target.names] <- key[new.data.names]

  }

  # Check if 'm.mat' exists
  if(!exists("m.mat")){
    m.mat <- NULL
  }else{
    m.mat <- t(m.mat)

    if(reduce.method == "latent"){
      colnames(m.mat) <- c("Target", paste("Redundancy_", 1:(ncol(m.mat)-1), sep = ""))
    }else{
      colnames(m.mat) <- c(paste("Redundancy_", 1:ncol(m.mat), sep = ""))
    }
  }

  # Check for "sum"
  if(reduce.method == "sum"){

    # Reinstate new.data
    new.data <- as.data.frame(node.redundant.obj$data)

    # Collapse across rows
    for(i in 1:nrow(m.mat)){

      # Collapse
      collapse <- row.names(m.mat)[i]

      # Redundant
      redunds <- m.mat[i,]
      redunds <- redunds[redunds != ""]

      # Obtain columns that exist in data
      extant_cols <- c(collapse, redunds)
      extant_cols <- extant_cols[extant_cols %in% key]
      extant_cols <- names(key)[match(extant_cols, key)]

      # Collapse and insert into matrix
      new.data[[collapse]] <- rowSums(new.data[,extant_cols], na.rm = TRUE)

      # Remove redundant terms
      new.data <- new.data[,-match(extant_cols, colnames(new.data))]
    }

    # Convert new.data back to matrix (for symmetric check)
    new.data <- as.matrix(new.data)

  }



  # Initialize results list
  res <- list()
  res$data <- new.data
  res$merged <- m.mat

  return(res)

}

#' @noRd
# Redundancy Reduction (Automated)
# Updated 13.09.2022
redund.reduce.auto <- function(node.redundant.obj,
                               reduce.method, lavaan.args, corr)
{
  # Check for node.redundant object class
  if(is(node.redundant.obj) != "node.redundant"){
    stop("A 'node.redundant' object must be used as input")
  }

  # Redundant list
  redund <- node.redundant.obj$redundant

  # Copied data
  new.data <- node.redundant.obj$data

  # Track merged items
  merged <- list()

  # Track changed names
  name.chn <- vector("character")

  # Initialize count
  count <- 0

  # Get key
  if("key" %in% names(node.redundant.obj)){

    key <- node.redundant.obj$key
    names(key) <- names(node.redundant.obj$key)

  }else{

    key <- colnames(node.redundant.obj$data)
    names(key) <- key

  }

  # Line break
  #linebreak()

  # Previous state
  prev.state <- list(redund)
  prev.state.data <- list(new.data)

  # Loop through named node redundant list
  while(length(redund) != 0){

    # Targeting redundancy
    target.item <- names(redund)[1]

    # Potential redundancies
    pot <- redund[[1]]

    # Get input
    input <- as.character(1:length(pot))

    # Auto
    if(all(input != "0")){

      # Convert to numeric
      re.items <- as.numeric(unlist(strsplit(unlist(strsplit(input, split = " ")), split = ",")))

      # Items to combine with target
      comb <- pot

      # Index items
      idx <- names(key)[match(comb, key)]

      # Target index
      tar.idx <- names(key)[match(target.item, key)]

      # Update merged list
      count <- count + 1
      merged[[count]] <- c(key[tar.idx], key[idx])

      # Combine into target index
      if(reduce.method == "latent"){

        # Latent variable
        ## create model
        mod <- paste(paste("comb =~ ",sep=""), paste(colnames(new.data[,c(tar.idx, idx)]), collapse = " + "))

        # Replace arguments
        lavaan.args$model <- mod
        lavaan.args$data <- new.data
        
        ## Get default estimator
        categories <- apply(new.data[,c(tar.idx, idx)], 2, function(x){
          length(unique(x))
        })

        # Check categories
        if(sum(categories < 6) > 1){# Not all continuous
          lavaan.args$estimator <- "WLSMV"
          lavaan.args$missing <- "pairwise"
          lavaan.args$ordered <- TRUE
        }else{# All can be considered continuous
          lavaan.args$estimator <- "MLR"
          lavaan.args$missing <- "fiml"
          lavaan.args$ordered <- FALSE
        }

        ## get CFA function from lavaan
        FUN <- lavaan::cfa

        ## fit model
        fit <- suppressWarnings(
          suppressMessages(
            do.call(what = "FUN", args = as.list(lavaan.args))
          )
        )

        ## identify cases
        cases <- lavaan::inspect(fit, "case.idx")

        ## compute latent variable score
        latent <- as.numeric(lavaan::lavPredict(fit))

        ## check for missing cases and handle
        if(length(cases) != nrow(new.data)){
          new.vec <- as.vector(matrix(NA, nrow = nrow(new.data), ncol = 1))
          new.vec[cases] <- latent
        }else{new.vec <- latent}

        ## check for reverse scoring/labelling
        corrs <- as.matrix(cor(cbind(latent,new.data[,c(tar.idx, idx)]), use = "complete.obs")[1,-1])
        row.names(corrs) <- c(key[tar.idx], key[idx])
        colnames(corrs) <- "latent"

        # input new vector
        new.data[,tar.idx] <- new.vec

        # Ask for new label
        lab <- paste("LV_", count, sep = "")
        name.chn[count] <- lab
        col.idx <- match(tar.idx, colnames(new.data))
        colnames(new.data)[col.idx] <- lab

      }

      # Remove redundant variables from data
      rm.idx <- match(idx, colnames(new.data))

      if(isSymmetric(new.data)){
        new.data <- new.data[-rm.idx, -rm.idx]
      }else{
        new.data <- new.data[,-rm.idx]
      }

      # Update previous state data
      prev.state.data[length(prev.state.data) + 1] <- list(new.data)

      # Remove target item
      redund[[1]] <- NULL

      # Make ind
      if(reduce.method == "latent" | reduce.method == "sum"){ind <- idx}

      # Remove variables within future options
      for(i in 1:length(ind)){

        ## Get elements from other redundancy lists
        elements <- sapply(redund, function(x){key[ind[i]] %in% x})

        ## If there are any
        if(any(elements)){

          ### Target elements
          target.elem <- which(elements)

          for(j in 1:length(target.elem)){

            #### Target list
            list.target <- redund[[target.elem[j]]]

            #### Remove from target list
            redund[[target.elem[j]]] <- redund[[target.elem[j]]][-which(key[ind[i]] == list.target)]

          }

        }

      }

      # Remove empty list objects
      rm.list <- which(
        unlist(lapply(redund, function(x){
          length(x)
        })) == 0
      )

      if(length(rm.list) != 0){
        redund <- redund[-rm.list]
      }

      # Remove object names
      if(any(key[ind] %in% names(redund))){

        ## Get target names
        name.targets <- which(key[ind] %in% names(redund))

        ## Loop through
        for(i in 1:length(name.targets)){

          # If there is only one item left, then remove from list
          if(length(redund[[key[ind][name.targets[i]]]]) == 1){
            redund[key[ind][name.targets[i]]] <- NULL
          }else{# Otherwise, replace name with first element and remove first element from list
            names(redund[key[ind]][name.targets[i]]) <- redund[[key[ind][name.targets[i]]]][1]
            redund[[key[ind][name.targets[i]]]][1] <- NA
            redund[[key[ind][name.targets[i]]]] <- na.omit(redund[[key[ind][name.targets[i]]]])
            redund[key[ind][name.targets[i]]] <- NULL
          }

        }

      }

      # Update previous state
      prev.state[length(prev.state) + 1] <- list(redund)

    }else{

      # Update previous state data
      prev.state.data[length(prev.state.data) + 1] <- list(new.data)

      redund[[1]] <- NULL

      # Update previous state
      prev.state[length(prev.state) + 1] <- list(redund)

    }

    if(!is.null(input)){
      #linebreak()
      input <- NULL
    }

    # Artificial pause for smoothness of experience
    #Sys.sleep(1)

  }

  # Transform merged list to matrix
  if(length(merged) != 0){

    # Number of rows for matrix
    m.rows <- max(unlist(lapply(merged, length)))

    # Initialize merged matrix
    m.mat <- matrix("", nrow = m.rows, ncol = length(merged))

    # Input into merged matrix
    for(i in 1:length(merged))
    {
      diff <- m.rows - length(merged[[i]])

      m.mat[,i] <- c(merged[[i]], rep("", diff))
    }

    # Set names for redundancy matrix
    if(length(name.chn) != 0){
      colnames(m.mat) <- name.chn
    }else{
      colnames(m.mat) <- paste("LV", 1:length(merged), sep = "_")
    }

  }

  # Replace column names for item names not changed
  if(any(colnames(new.data) %in% names(key))){

    # Target names
    target.names <- which(colnames(new.data) %in% names(key))

    # new.data names
    new.data.names <- colnames(new.data)[target.names]

    # Insert into new data
    colnames(new.data)[target.names] <- key[new.data.names]

  }

  # Check if 'm.mat' exists
  if(!exists("m.mat")){
    m.mat <- NULL
  }else{
    m.mat <- t(m.mat)

    if(reduce.method == "latent" | reduce.method == "sum"){
      colnames(m.mat) <- c("Target", paste("Redundancy_", 1:(ncol(m.mat)-1), sep = ""))
    }else if(reduce.method == "remove"){
      colnames(m.mat) <- c(paste("Redundancy_", 1:ncol(m.mat), sep = ""))
    }
  }

  # Check for "sum"
  if(reduce.method == "sum"){

    # Reinstate new.data
    new.data <- node.redundant.obj$data

    # Get key
    if("key" %in% names(node.redundant.obj)){

      key <- node.redundant.obj$key
      names(key) <- names(node.redundant.obj$key)

    }else{

      key <- colnames(node.redundant.obj$data)
      names(key) <- key

    }

    # Replace column names for item names not changed
    if(any(colnames(new.data) %in% names(key))){

      # Target names
      target.names <- which(colnames(new.data) %in% names(key))

      # new.data names
      new.data.names <- colnames(new.data)[target.names]

      # Insert into new data
      colnames(new.data)[target.names] <- key[new.data.names]

    }

    # Convert to data frame
    new.data <- as.data.frame(new.data)

    # Collapse across rows
    for(i in 1:nrow(m.mat)){

      # Collapse
      collapse <- row.names(m.mat)[i]

      # Redundant
      redunds <- m.mat[i,]
      redunds <- redunds[redunds != ""]

      # Obtain columns that exist in data
      # extant_cols <- c(collapse, redunds)
      # extant_cols[extant_cols %in% key] <- names(key)[match(extant_cols[extant_cols %in% key], key)]

      # Collapse and insert into matrix
      new.data[[collapse]] <- rowSums(new.data[,redunds], na.rm = TRUE)

      # Remove redundant terms
      new.data <- new.data[,-match(redunds, colnames(new.data))]
    }

    # Convert new.data back to matrix (for symmetric check)
    new.data <- as.matrix(new.data)

  }

  # Initialize results list
  res <- list()
  res$data <- new.data
  res$merged <- m.mat

  return(res)

}

#' @noRd
# Redundancy Adhoc Reduction (Automated)
# Updated 20.07.2022
redund.adhoc.auto <- function(node.redundant.obj,
                              node.redundant.reduced,
                              node.redundant.original,
                              reduce.method, lavaan.args, corr)
{

  # Redundant list
  redund <- node.redundant.obj$redundant

  # Check for overlaps
  overlap <- unlist(
    lapply(redund, function(x, name){
      any(name %in% x)
    }, name = names(redund))
  )

  # Merge overlaps
  if(length(which(overlap)) != 0){

    for(i in which(overlap)){
      redund[[i]] <- unique(unlist(c(redund[[i]], redund[names(redund) %in% redund[[i]]])))

      redund <- redund[-which(names(redund) %in% redund[[i]])]
    }

  }
  
  # Check for other overlaps
  for(i in 1:length(redund)){
    
    # Identify any overlap
    target <- any(redund[[i]] %in% unlist(redund[-i]))
    
    # Remove latter overlap
    if(isTRUE(target)){
      
      # Obtain matched target
      matched <- unlist(redund[-i])[match(redund[[i]], unlist(redund[-i]))]
      matched <- matched[!is.na(matched)]
      
      # Remove from each list
      for(j in 1:length(matched)){
        
        # Target list
        target_list <- redund[[names(matched)[j]]]
        target_list[which(target_list == matched[j])] <- NA
        
        # Return target list
        redund[[names(matched)[j]]] <- na.omit(target_list)
      }
      
    }
    
  }
  
  # Remove empty lists
  lengths <- unlist(lapply(redund, length))
  if(any(lengths == 0)){
    redund <- redund[-which(lengths == 0)]
  }

  # Copied data
  new.data <- node.redundant.original$data

  # Get reduced
  reduced.merged <- node.redundant.reduced$merged

  # Make merged
  merged <- lapply(apply(reduced.merged, 1, as.list), function(x){unname(unlist(x))})

  # Set counter
  count <- max(as.numeric(gsub("LV_", "", names(merged))))
  original_count <- count

  # Update merged
  for(i in 1:length(redund)){

    # Target variables
    target <- c(names(redund)[i], redund[[i]])

    # Expand latent variables
    lv <- target[grep("LV_", target)]

    # If any latent variables
    if(length(lv) != 0){

      # Create new latent variable
      new_LV <- unname(c(unlist(merged[lv]), target[-grep("LV_", target)]))

      # Insert into merged
      count <- count + 1
      merged[[paste("LV_", count, sep = "")]] <- new_LV

      # Remove old latent variables
      merged[lv] <- NULL

    }else{

      # Create new latent variable
      new_LV <- target

      # Insert into merged
      count <- count + 1
      merged[[paste("LV_", count, sep = "")]] <- new_LV

      # Remove old variables
      merged[names(redund)[i]] <- NULL

    }
  }

  # Remove all ""
  merged <- lapply(merged, function(x){
    x <- na.omit(ifelse(x == "", NA, x))
    attr(x, which = "na.action") <- NULL
    return(x)
  })

  # Reorder longest to shortest
  merged <- merged[order(unlist(lapply(merged, length)), decreasing = TRUE)]

  # Get key
  if("key" %in% names(node.redundant.original)){
    key <- node.redundant.original$key
    names(key) <- names(node.redundant.original$key)
  }else{
    key <- colnames(node.redundant.original$data)
    names(key) <- key
  }

  # Remove missing redundancies
  merged <- lapply(merged, function(x){
    if(length(x) == 0){
      NULL
    }else{x}
  })

  nulls <- unlist(lapply(merged, is.null))

  if(any(nulls)){
    merged <- merged[!nulls]
  }

  # Loop through to make new variables
  for(i in 1:length(merged)){

    # Combine into target index
    if(reduce.method == "latent"){

      # Get indices from key
      idx <- names(key[match(merged[[i]], key)])

      # Create model
      mod <- paste(paste("comb =~ ",sep=""), paste(colnames(new.data[,idx]), collapse = " + "))

      # Replace arguments
      lavaan.args$model <- mod
      lavaan.args$data <- new.data
      ## Get default estimator
      categories <- apply(new.data[,idx], 2, function(x){
        length(unique(x))
      })

      # Check categories
      if(sum(categories < 6) > 1){# Not all continuous
        lavaan.args$estimator <- "WLSMV"
        lavaan.args$missing <- "pairwise"
        lavaan.args$ordered <- TRUE
      }else{# All can be considered continuous
        lavaan.args$estimator <- "MLR"
        lavaan.args$missing <- "fiml"
        lavaan.args$ordered <- FALSE
      }

      ## get CFA function from lavaan
      FUN <- lavaan::cfa

      ## fit model
      fit <- suppressWarnings(
        suppressMessages(
          do.call(what = "FUN", args = as.list(lavaan.args))
        )
      )

      ## identify cases
      cases <- lavaan::inspect(fit, "case.idx")

      ## compute latent variable score
      latent <- as.numeric(lavaan::lavPredict(fit))

      ## check for missing cases and handle
      if(length(cases) != nrow(new.data)){
        new.vec <- as.vector(matrix(NA, nrow = nrow(new.data), ncol = 1))
        new.vec[cases] <- latent
      }else{new.vec <- latent}

      # Tack on latent variable
      ## Must come first!
      new.data <- cbind(new.data, new.vec)

      # Remove variables (ensure matrix)
      new.data <- as.matrix(new.data[,-match(idx, colnames(new.data))])

      # Rename latent variable
      colnames(new.data)[ncol(new.data)] <- names(merged)[i]

    }

  }

  # Transform merged list to matrix
  if(length(merged) != 0){

    # Number of rows for matrix
    m.rows <- max(unlist(lapply(merged, length)))

    # Initialize merged matrix
    m.mat <- matrix("", nrow = m.rows, ncol = length(merged))

    # Input into merged matrix
    for(i in 1:length(merged))
    {
      diff <- m.rows - length(merged[[i]])

      m.mat[,i] <- c(merged[[i]], rep("", diff))
    }

    # Set names for redundancy matrix
    colnames(m.mat) <- names(merged)

  }

  # Replace column names for item names not changed
  if(any(colnames(new.data) %in% names(key))){

    # Target names
    target.names <- which(colnames(new.data) %in% names(key))

    # new.data names
    new.data.names <- colnames(new.data)[target.names]

    # Insert into new data
    colnames(new.data)[target.names] <- key[new.data.names]

  }

  # Check if 'm.mat' exists
  if(!exists("m.mat")){
    m.mat <- NULL
  }else{
    m.mat <- t(m.mat)

    if(reduce.method == "latent" | reduce.method == "sum"){
      colnames(m.mat) <- c("Target", paste("Redundancy_", 1:(ncol(m.mat)-1), sep = ""))
    }else if(reduce.method == "remove"){
      colnames(m.mat) <- c(paste("Redundancy_", 1:ncol(m.mat), sep = ""))
    }
  }

  # Check for "sum"
  if(reduce.method == "sum"){

    # Reinstate new.data
    new.data <- as.data.frame(new.data)

    # Collapse across rows
    for(i in 1:nrow(m.mat)){

      # Collapse
      collapse <- row.names(m.mat)[i]

      # Redundant
      redunds <- m.mat[i,]
      redunds <- redunds[redunds != ""]

      # Obtain columns that exist in data
      # extant_cols <- c(collapse, redunds)
      # extant_cols[extant_cols %in% key] <- names(key)[match(extant_cols[extant_cols %in% key], key)]

      # Collapse and insert into matrix
      new.data[[collapse]] <- rowSums(new.data[,redunds], na.rm = TRUE)

      # Remove redundant terms
      new.data <- new.data[,-match(redunds, colnames(new.data))]
    }

    # Convert new.data back to matrix (for symmetric check)
    new.data <- as.matrix(new.data)

  }

  # Initialize results list
  res <- list()
  res$data <- new.data
  res$merged <- m.mat

  return(res)

}

#' @noRd
# Item-total correlations
# For UVA
# Updated 15.02.2021
item.total <- function (data.sub, corr)
{
  # Get correlations
  corrs <- switch(corr,
                  "cor_auto" = suppressMessages(qgraph::cor_auto(data.sub)),
                  "pearson" = suppressMessages(cor(data.sub, use = "pairwise.complete.obs")),
                  "spearman" = suppressMessages(cor(data.sub, method = "spearman", use = "pairwise.complete.obs"))
  )

  # Check for negatives (reverse if so)
  for(i in 1:nrow(corrs)){

    if(any(sign(corrs[i,-i]) == -1)){
      data.sub[,i] <- (max(data.sub[,i]) + min(data.sub[,i])) - data.sub[,i]

      if(any(sign(corrs[i,-i]) == 1)){
        target <- colnames(corrs)[(sign(corrs[i,]) == 1)[-i]]

        for(j in 1:length(target)){
          data.sub[,target[j]] <- (max(data.sub[,target[j]]) + min(data.sub[,target[j]])) - data.sub[,target[j]]
        }

      }

      corrs <- switch(corr,
                      "cor_auto" = suppressMessages(qgraph::cor_auto(data.sub)),
                      "pearson" = suppressMessages(cor(data.sub, use = "pairwise.complete.obs")),
                      "spearman" = suppressMessages(cor(data.sub, method = "spearman", use = "pairwise.complete.obs"))
      )
    }

  }

  # Compute corrected item-total correlations
  cor.corr <- numeric(ncol(data.sub))

  # Loop through
  for(i in 1:ncol(data.sub)){

    cor.corr[i] <- switch(corr,
                          "cor_auto" = suppressMessages(qgraph::cor_auto(cbind(data.sub[,i], rowSums(data.sub[,-i]))))[1,2],
                          "pearson" = suppressMessages(cor(cbind(data.sub[,i],
                                                                 rowSums(data.sub[,-i])),
                                                           method = "pearson",
                                                           use = "pairwise.complete.obs"))[1,2],
                          "spearman" = suppressMessages(cor(cbind(data.sub[,i],
                                                                  rowSums(data.sub[,-i])),
                                                            method = "spearman",
                                                            use = "pairwise.complete.obs"))[1,2]
    )

  }

  return(cor.corr)

}

#' @noRd
# Menu for redundancy
# For UVA
# Updated 15.02.2021
redundancy.menu <- function (redund, reduce.method, pot, target.item, weights,
                             plot.args, key, node.redundant.obj)
{
  # Configure into list
  pot <- list(pot)
  names(pot) <- target.item
  possible <- unname(unlist(pot))

  # Loop through potential redundancies
  count2 <- 1

  # Check names
  for(i in 1:length(possible)){

    if(possible[i] %in% names(redund)){
      count2 <- count2 + 1
      pot[count2] <- redund[possible[i]]
      names(pot)[count2] <- possible[i]
    }

  }

  # Check elements
  for(i in 1:length(possible)){
    elements <- redund[sapply(redund, function(x){possible[i] %in% x})]
    pot[(count2 + 1):(count2 + length(elements))] <- elements
    names(pot)[(count2 + 1):(count2 + length(elements))] <- names(elements)
    count2 <- count2 + length(elements)
  }

  # Get unique lists
  pot <- redund[unique(names(pot))]

  # Possible options
  poss <- unique(c(unname(unlist(pot)), names(pot)[-1]))

  # Organize plot of redundancy connections
  mat <- matrix(0, nrow = length(poss) + 1, ncol = length(poss) + 1)
  colnames(mat) <- c(paste("Target"), 1:length(poss))
  row.names(mat) <- colnames(mat)

  mat["Target",paste(1:length(unlist(pot[[1]])))] <- weights[names(key[match(target.item, key)]),names(key[match(unlist(pot[[1]]),key)])]
  mat[paste(1:length(unlist(pot[[1]]))),"Target"] <- weights[names(key[match(target.item, key)]),names(key[match(unlist(pot[[1]]),key)])]

  if(length(pot) != 1)
  {
    # Remove first element
    ext <- pot[-1]

    # Loop through rest of extended
    for(i in 1:length(ext))
    {
      # Target extended
      target.ext <- ext[[i]]

      # Loop through target
      for(j in 1:length(target.ext))
      {
        # Single out each element
        single <- target.ext[[j]]

        # Get element in possible redundancies
        elem <- match(names(ext)[i], poss)

        # Get elements redundant with element
        red.elem <- match(single, poss)

        # Put into matrix
        mat[paste(elem),paste(red.elem)] <- weights[names(key[match(poss[elem], key)]),names(key[match(poss[red.elem],key)])]
        mat[paste(red.elem),paste(elem)] <- weights[names(key[match(poss[elem], key)]),names(key[match(poss[red.elem],key)])]
      }
    }
  }

  # Print target and potential options
  cat(paste("Target variable: '", target.item, "'", sep = ""))
  cat("\n\nPotential redundancies:\n\n")
  if(reduce.method == "latent"){
    cat("0. Do not combine with any")
  }else if(reduce.method == "remove" | reduce.method == "sum"){
    cat("0. None")
  }

  cat(paste("\n", 1:length(poss), ". ", "'", poss, "'", sep = ""),"\n")

  # Plot
  if(node.redundant.obj$model == "tmfg"){

    plot.args$title <- "Zero-order Correlations"

  }else{

    plot.args$title <- switch(node.redundant.obj$method,
                              "wto" = "Regularized Partial Correlations",
                              "pcor" = "Partial Correlations",
                              "cor" = "Zero-order Correlations",
                              "irt" = "Correlated Residuals"
    )

  }

  if(length(poss) > 1){
    # Global suppress warnings (need a better work around)
    warn <- options("warn")[[1]]
    options(warn = -1)
    plot(redund.plot(plot.matrix = mat, plot.args = plot.args, plot.reduce = TRUE))
    options(warn = warn)
  }else{

    if(node.redundant.obj$model == "tmfg"){

      plot.args$title <- "Zero-order Correlation"

    }else{

      plot.args$title <- switch(node.redundant.obj$method,
                                "wto" = "Regularized Partial Correlation",
                                "pcor" = "Partial Correlation",
                                "cor" = "Zero-order Correlation",
                                "irt" = "Correlated Residuals"
      )

    }

    par(mar = c(0,0,0,0))
    plot(c(0, 1), c(0, 1), ann = FALSE, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
    text(x = 0.5, y = 0.5, paste("There is only one redundant variable with the target variable.\nTheir ",
                                 tolower(plot.args$title), " = ", round(mat[1,2], 3),
                                 sep = ""),
         cex = 1.6, col = "black")
    par(mar = c(5, 4, 4, 2) + 0.1)

  }

  message("\nPress 'B' to go back")

  return(poss)
}

#' @noRd
# Input check for redundancy
# For UVA
# Updated 21.12.2020
input.check <- function (poss, type = c("redund", "remove"))
{

  if(type == "redund"){

    message("\nEnter numbers of variables redundant with the target variable (separate by commas)")
    input <- readline(prompt = "Selection: ")

    # Redo input check
    re.input <- in.check(input, poss = poss)

    while(re.input){
      # Print message to try again
      message("Inappropriate input. Try again. 'B' can be used to go back.\n")

      # Get input
      message("Enter numbers of variables redundant with the target variable (separate by commas)")
      input <- readline(prompt = "Selection: ")

      # Redo input check
      re.input <- in.check(input, poss)
    }

  }else if (type == "remove"){

    cat(
      paste("\n0. ", "'", poss[1], "'", " (Target)", sep = ""),
      paste("\n", 1:(length(poss) - 1), ". ", "'", poss[-1], "'", sep = ""),
      "\n\n"
    )

    input <- readline(prompt = "Select variable to KEEP: ")

    # Redo input check
    re.input <- in.check(input, poss = poss)

    while(re.input)
    {
      # Print message to try again
      message("Inappropriate input. Try again.\n")

      # Get input
      input <- readline(prompt = "Select variable to KEEP: ")

      # Redo input check
      re.input <- in.check(input, poss)
    }

  }

  return(input)
}

#' @noRd
# Input check for redundancy
# For UVA
# Updated 21.12.2020
in.check <- function(input, poss)
{
  if(tolower(input) == "b"){
    ret.val <- FALSE
  }else{
    inp <- suppressWarnings(as.numeric(unlist(strsplit(unlist(strsplit(input, split = " ")), split = ","))))

    ret.val <- FALSE

    if(any(is.na(inp)))
    {ret.val <- TRUE}

    if(length(inp) == 0)
    {ret.val <- TRUE}

    if(length(setdiff(inp, 0:length(poss))) != 0)
    {ret.val <- TRUE}
  }

  return(ret.val)
}

#' @noRd
# Creates a nice looking line break
# For UVA
# Updated 26.02.2021
linebreak <- function(){cat("\n", colortext(paste(rep("-", getOption("width")), collapse = ""), defaults = "message"), "\n\n")}

#' @noRd
# Changes names that have leading characters that are numbers
# This creates an error in lavaan's formulas
# Numbers are moved to the end of the name
# For UVA
# Updated 24.03.2021
lavaan.formula.names <- function (data){

  # Original column names
  original.names <- colnames(data)

  # Move leading numeric values to the end of the variable name
  colnames(data) <- unlist(
    lapply(strsplit(colnames(data),
                    split = ""), function(x){

                      ind <- grepl("[[:digit:]]", x)

                      if(isTRUE(rle(ind)$values[1])){
                        rm.ind <- x[1:rle(ind)$lengths[1]]
                        new.name <- paste(
                          paste(x[-c(1:rle(ind)$lengths[1])], collapse = ""),
                          paste(rm.ind, collapse = ""),
                          sep = "_"
                        )
                        return(new.name)
                      }else{return(paste(x, collapse = ""))}

                    })
  )

  # Message user so they know
  if(any(colnames(data) != original.names)){

    message(paste("Some variable names begin with a number.",
                  "This creates errors in 'lavaan' formulas.",
                  styletext("\nAll numbers have been moved to the end of the variable name.", defaults = "bold"),
                  "\nTo avoid this message, make sure all variable names begin with a character.",
                  sep = "\n"))

  }

  return(data)

}

#%%%%%%%%%%%%%%%%%%%
# itemStability ----
#%%%%%%%%%%%%%%%%%%%

#' @noRd
# Converts memberships to numeric
# For itemStability
# Updated 26.02.2021
numeric.membership <- function(membership){

  # Get unique membership
  unique.membership <- unique(membership)

  # Check for named memberships
  if(all(is.character(unique.membership))){

    # Assign numeric values to memberships
    membership.numbers <- membership

    for(i in 1:length(unique.membership)){
      membership.numbers[which(membership.numbers == unique.membership[i])] <- i
    }

  }else{membership.numbers <- membership}

  # Make sure the memberships are numeric
  membership.numbers <- as.numeric(membership.numbers)

  # Name elements
  names(membership.numbers) <- names(membership)

  return(membership.numbers)

}

#' Computes the mode
#'
#' @param v Numeric vector.
#' Vector of values to find mode in
#'
#' @param fin.vec Alphanumeric vector.
#' Vector of current state of \code{v}
#'
#' @return The mode of a vector
#'
#' @noRd
#'
# Mode
# Updated 15.06.2020
mode <- function(v, fin.vec)
{
  #unique values
  uniqv <- unique(v)

  #find mode
  uniq.val <- uniqv[which.max(tabulate(match(v, uniqv)))]

  #do not overwrite already identified dimension
  while(uniq.val %in% fin.vec)
  {
    #remove unique value
    uniqv <- uniqv[-which(uniq.val==uniqv)]

    if(length(uniqv)==0)
    {
      uniq.val <- NA
      break
    }

    #find mode
    uniq.val <- uniqv[which.max(tabulate(match(v, uniqv)))]
  }

  return(uniq.val)
}

#' Converts membership vector into a target membership vector
#'
#' @param target.wc Numeric vector.
#' Target membership vector
#'
#' @param convert.wc Numeric vector, matrix, or data frame.
#' Memberships to be converted to the target membership vector
#'
#' @return Vector or matrix of homogenized memberships
#'
#' @noRd
#'
# Homogenize Membership
# For itemStability
# Updated 01.03.2022
homogenize.membership <- function (target.wc, convert.wc)
{
  # Obtain whether vector or matrix is input for 'convert.wc'
  if(is.vector(convert.wc)){

    # Error if not same length
    if(length(convert.wc) != length(target.wc))
    {stop("'convert.wc' must be same length as 'target.wc'")}

    # Number of comparisons
    n <- 1

  }else{

    # Make sure 'convert.wc' is matrix
    convert.wc <- as.matrix(convert.wc)

    # Number of comparisons
    n <- ncol(convert.wc)

  }

  # Initialize conversion matrix
  convert.mat <- matrix(NA, nrow = length(target.wc), ncol = n)
  ## Get node names
  if(!is.null(names(target.wc))){
    row.names(convert.mat) <- names(target.wc)
  }

  # Check target wc for NAs
  if(any(is.na(target.wc))){
    na.wc <- which(is.na(target.wc))
    target.wc <- target.wc[-na.wc]
  }

  # Identify target membership within bootstrapped memberships
  for(i in 1:n){

    # New membership vector
    new.vec <- convert.wc[,i]

    # Remove NAs from new membership vector
    if(exists("na.wc")){
      new.vec <- new.vec[-na.wc]
    }

    # Unique new membership
    new.uniq <- unique(new.vec)

    # Check if dimensionality solution was reached
    if(length(unique(na.omit(new.vec))) == length(new.vec)){

      # Submit NAs
      final.vec <- rep(NA, length = length(target.wc))
      names(final.vec) <- names(target.wc)

    }else if(all(target.wc == new.vec, na.rm = TRUE)){
      # Check if all dimensions are the same
      final.vec <- new.vec
      names(final.vec) <- names(target.wc)
    }else if(length(na.omit(unique(target.wc))) > length(unique(na.omit(new.vec)))){
      # Converge based on maximum number of dimensions

      # Initialize rand and length vector
      rand <- vector("numeric", length = length(unique(na.omit(new.vec))))
      names(rand) <- na.omit(new.uniq)
      len <- rand

      for(j in new.uniq){

        # Target nodes
        target <- which(new.vec==j)

        # Lengths of target
        len[paste(j)] <- length(target)

        # Compute rand index
        if(length(target.wc[target]) == length(unique(target.wc[target]))){
          rand[paste(j)] <- ifelse(
            length(unique(target.wc[target])) == 1, 1, 0
          )
        }else if(length(new.vec[target]) == length(unique(new.vec[target]))){
          rand[paste(j)] <- ifelse(
            length(unique(new.vec[target])) == 1, 1, 0
          )
        }else{
          rand[paste(j)] <- igraph::compare(new.vec[target],target.wc[target],method="rand")
        }

      }

      # Remove NAs
      rand <- na.omit(rand)
      rand <- rand[!is.NA(names(rand))]
      len <- na.omit(ifelse(len == 0, NA, len))

      # Order rand by highest rand index and then number of items
      rand.ord <- rand[order(rand, len, decreasing = TRUE)]

      # Initialize final vector
      final.vec <- rep(NA, length = length(target.wc))
      names(final.vec) <- names(target.wc)

      # Insert new values into final vector
      for(j in as.numeric(names(rand.ord))){

        # Identify target
        new.target <- which(new.vec==j)

        # Identify mode
        target.mode <- mode(target.wc[new.target], final.vec)

        # Insert into final vector
        final.vec[new.target] <- rep(target.mode)
      }

    }else if(length(na.omit(unique(target.wc))) < length(unique(na.omit(new.vec)))){

      # Initialize rand and length vector
      rand <- vector("numeric", length = length(unique(na.omit(new.vec))))
      names(rand) <- na.omit(new.uniq)
      len <- rand

      for(j in new.uniq){

        # Target nodes
        target <- which(new.vec==j)

        # Lengths of target
        len[paste(j)] <- length(target)

        # Compute rand index
        if(length(target.wc[target]) == length(unique(target.wc[target]))){
          rand[paste(j)] <- ifelse(
            length(unique(target.wc[target])) == 1, 1, 0
          )
        }else if(length(new.vec[target]) == length(unique(new.vec[target]))){
          rand[paste(j)] <- ifelse(
            length(unique(new.vec[target])) == 1, 1, 0
          )
        }else{
          rand[paste(j)] <- igraph::compare(new.vec[target],target.wc[target],method="rand")
        }

      }

      # Remove NAs
      rand <- na.omit(rand)
      rand <- rand[!is.NA(names(rand))]
      len <- na.omit(ifelse(len == 0, NA, len))

      # Order rand by highest rand index and then number of items
      rand.ord <- rand[order(rand, len, decreasing = TRUE)]

      # Initialize final vector
      final.vec <- rep(NA, length = length(target.wc))
      names(final.vec) <- names(target.wc)

      # Insert new values into final vector
      for(j in as.numeric(names(rand.ord))){

        # Identify target
        new.target <- which(new.vec==j)

        # Identify mode
        target.mode <- mode(target.wc[new.target], final.vec)

        # Insert into final vector
        final.vec[new.target] <- rep(target.mode)
      }

      # Identify number of extra dimensions
      extra.dim <- unique(new.vec[which(is.na(final.vec))])

      # Initialize extra dimension length vector
      extra.len <- rep(NA, length = length(extra.dim))
      names(extra.len) <- extra.dim

      # Initialize count
      count <- 0

      # Order length of extra dimensions
      for(j in extra.dim){

        # Increase count
        count <- count + 1

        # Length of extra dimensions
        extra.len[count] <- length(which(new.vec==j))
      }

      el.ord <- extra.len[order(extra.len, decreasing = TRUE)]

      # Reset count
      count <- 0

      # Insert extra dimensions into final vector
      for(j in 1:length(el.ord)){

        # Increase count
        count <- count + 1

        # Target extra dimension
        target.ed <- as.numeric(names(el.ord)[j])

        # Insert dimensions into final vector
        final.vec[which(new.vec==target.ed)] <- (max(target.wc) + count)
      }

    }else{

      # Initialize rand and length vector
      rand <- vector("numeric", length = length(unique(na.omit(new.vec))))
      names(rand) <- na.omit(new.uniq)
      len <- rand

      for(j in new.uniq){

        # Target nodes
        target <- which(new.vec==j)

        # Lengths of target
        len[paste(j)] <- length(target)

        # Compute rand index
        if(length(target.wc[target]) == length(unique(target.wc[target]))){
          rand[paste(j)] <- ifelse(
            length(unique(target.wc[target])) == 1, 1, 0
          )
        }else if(length(new.vec[target]) == length(unique(new.vec[target]))){
          rand[paste(j)] <- ifelse(
            length(unique(new.vec[target])) == 1, 1, 0
          )
        }else{
          rand[paste(j)] <- igraph::compare(new.vec[target], target.wc[target], method="rand")
        }

      }

      # Remove NAs
      rand <- na.omit(rand)
      rand <- rand[!is.NA(names(rand))]
      len <- na.omit(ifelse(len == 0, NA, len))

      # Order rand by highest rand index and then number of items
      rand.ord <- rand[order(rand, len, decreasing = TRUE, na.last = FALSE)]

      # Initialize final vector
      final.vec <- rep(NA, length = length(target.wc))
      names(final.vec) <- names(target.wc)

      # Insert new values into final vector
      for(j in as.numeric(names(rand.ord)))
      {
        # Identify target
        new.target <- which(new.vec==j)

        # Identify mode
        target.mode <- mode(target.wc[new.target], final.vec)

        # Insert into final vector
        final.vec[new.target] <- rep(target.mode)
      }
    }

    # Insert final vector into final matrix
    if(exists("na.wc")){
      convert.mat[-na.wc,i] <- final.vec
    }else{
      convert.mat[,i] <- final.vec
    }
  }

  return(convert.mat)
}

#' @noRd
# Checks for missing dimensions
# For itemStability
# Updated 26.02.2021
missing.dimension.check <- function (proportion, membership, bootstrap)
{
  # Get names of proportion
  proportion.names <- as.numeric(colnames(proportion))

  # Check for missing dimensions
  missing.dimensions <- setdiff(membership, proportion.names)

  # If there are missing dimensions, then add them (all zeroes)
  if(length(missing.dimensions) != 0){

    for(i in seq_along(missing.dimensions)){

      # Append missing dimension
      proportion <- cbind(proportion, rep(0, nrow(proportion)))
      colnames(proportion)[ncol(proportion)] <- paste(missing.dimensions[i])

    }

  }

  # Check for NA in bootstrap membership
  if(any(is.na(bootstrap))){

    if(any(colnames(proportion) == "NA")){
      proportion[,"NA"] <- colMeans(apply(bootstrap, 1, is.na))
    }else{

      # Get NA proportions
      NA.proportions <- colMeans(apply(bootstrap, 1, is.na))

      # Append proportion
      proportion <- cbind(proportion, NA.proportions)
      colnames(proportion)[ncol(proportion)] <- "NA"

    }

  }

  return(proportion)

}

#' @noRd
# Plot configuration for itemStability
# For itemStability
# Updated 17.07.2021
itemStability.plot <- function (res, bootega.obj)
{
  # Obtain empirical membership
  empirical.membership <- res$membership$empirical

  # Obtain unique membership
  unique.membership <- res$membership$unique

  # Obtain empirical stability
  empirical.stability <- res$item.stability$empirical.dimensions

  # Organize plot
  organize.plot <- data.frame(Item = names(empirical.membership),
                              Replication = empirical.stability,
                              Community = factor(empirical.membership,
                                                 unique.membership[order(
                                                   unique.membership
                                                 )]))

  # Item stability plot
  IS.plot <- ggpubr::ggdotchart(organize.plot, x = "Item", y = "Replication",
                                group = "Community", color = "Community",
                                legend.title = "Empirical EGA Communities",
                                add = "segments",
                                rotate = TRUE,
                                dot.size = 6,
                                label = round(organize.plot$Replication, 2),
                                font.label = list(color = "black", size = 8,
                                                  vjust = 0.5),
                                ggtheme = ggpubr::theme_pubr()
  )

  # Adjust y-axis and legend title
  IS.plot <- IS.plot + ggplot2::ylim(c(0,1)) + ggplot2::theme(
    legend.title = ggplot2::element_text(face = "bold"),
    axis.title = ggplot2::element_text(face = "bold")
  )

  # Manually change alpha
  IS.plot$layers[[2]]$aes_params$alpha <- 0.7

  # Adjust item label sizes based on
  sizes <- seq(6,12,.25)
  ## Number of nodes
  nodes <- rev(seq(0, 200, length.out = length(sizes)))
  n.size <- min(which(length(empirical.membership) > nodes))
  ## Number of characters in item name
  chars <- rev(seq(0,100, length.out = length(sizes)))
  ### Maximum characters in item name
  max.chars <- max(unlist(lapply(row.names(organize.plot),nchar)))
  c.size <- min(which(max.chars > chars))
  # Text size
  text.size <- sizes[min(c(n.size,c.size))]

  # Change text size
  IS.plot <- IS.plot +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size=text.size))

  # Change color.palette (if necessary)
  if(!ggplot2::is.ggplot(bootega.obj$plot.typical.ega)){

    if(!is.null(bootega.obj$color.palette)){

      IS.plot <- suppressMessages(
        IS.plot + ggplot2::scale_color_manual(values = color_palette_EGA(bootega.obj$color.palette,
                                                                         empirical.membership,
                                                                         sorted = TRUE),
                                              breaks = sort(empirical.membership))
      )

    }else{

      IS.plot <- suppressMessages(
        IS.plot + ggplot2::scale_color_manual(values = color_palette_EGA("rainbow",
                                                                         empirical.membership,
                                                                         sorted = TRUE),
                                              breaks = sort(empirical.membership))
      )

    }

  }else{
    if(bootega.obj$color.palette != "Set1"){
      IS.plot <- suppressMessages(
        IS.plot + ggplot2::scale_color_manual(values = color_palette_EGA(bootega.obj$color.palette,
                                                                         empirical.membership,
                                                                         sorted = TRUE),
                                              breaks = sort(empirical.membership, na.last = TRUE)
        )
      )
    }
  }

  # Reverse ordering
  IS.plot <- IS.plot + ggplot2::scale_x_discrete(limits = rev(IS.plot$data$Item))

  # Insert plot into results
  res$plot <- IS.plot

  return(res)
}

#' @noRd
# Average network loadings for itemStability
# For itemStability
# Updated 14.04.2021
itemStability.loadings <- function(res, bootega.obj)
{
  # Get graphs
  graphs <- bootega.obj$bootGraphs

  # Get memberships
  memberships <- res$membership$bootstrap

  # Maximum number of dimensions
  max.dimensions <- max(memberships, na.rm = TRUE)

  # Get number of iterations
  iterations <- ncol(memberships)

  # Make list of indices
  index.list <- as.list(1:iterations)

  # Loop through (for loop is slightly faster than lapply)
  loadings <- list()

  for(i in 1:iterations){

    if(all(is.na(memberships[,i]))){
      loadings[[i]] <- NULL
    }else{
      loadings[[i]] <- net.loads(
        A = graphs[[i]],
        wc = memberships[,i]
      )$std
    }

  }

  # Remove NULL loadings
  if(any(unlist(lapply(loadings, is.null)))){
    loadings <- loadings[-which(unlist(lapply(loadings, is.null)))]
  }

  # Initialize final loadings array
  loadings.array <- array(NA,
                          dim = c(nrow(memberships), # number of items
                                  max(memberships, na.rm = TRUE), # number of dimensions
                                  iterations), # number of bootstrap replicates
                          dimnames = list(row.names(memberships),
                                          1:max.dimensions,
                                          NULL))

  # Loop through
  for(i in 1:length(loadings)){

    loadings.array[row.names(loadings[[i]]), # get available loadings
                   colnames(loadings[[i]]), # get available dimensions
                   i] <- as.matrix(loadings[[i]]) # insert loadings into array

  }

  # Obtain average loadings
  mean.loadings <- apply(loadings.array, 1:2, mean, na.rm = TRUE)

  return(mean.loadings)
}

#%%%%%%%%%%%%%%%%%%%
# boot.ergoInfo ----
#%%%%%%%%%%%%%%%%%%%

#' Expand grid with unique edges
#' @noRd
# Updated 07.07.2022
expand.grid.unique <- function(x, y, include.equals = FALSE)
{
  x <- unique(x)
  
  y <- unique(y)
  
  g <- function(i)
  {
    z <- setdiff(y, x[seq_len(i-include.equals)])
    
    if(length(z)) cbind(x[i], z, deparse.level=0)
  }
  
  do.call(rbind, lapply(seq_along(x), g))
}

#' Rewiring function
#' @noRd
# Updated 16.09.2022
rewire <- function(network, min = 0.20, max = 0.40, noise = 0.10)
{
  
  # Number of edges
  edges <- sum(ifelse(network != 0, 1, 0)) / 2
  
  # Add noise
  if(!is.null(noise)){
    
    # Lower triangle of network
    lower_network <- network[lower.tri(network)]
    
    # Only add to existing edges
    lower_network[lower_network != 0] <- lower_network[lower_network != 0] +
      runif(
        n = edges,
        min = -noise,
        max = noise
      )
    
    # Replace lower network
    network[lower.tri(network)] <- lower_network
    
    # Replace upper network
    network <- t(network)
    network[lower.tri(network)] <- lower_network
    
  }
  
  # Set random proportion
  proportion <- runif(1, min = min, max = max)
  
  # Obtain proportion of connections to change
  rewire_number <- floor(edges * proportion)
  
  # Obtain edge list
  lower_network <- network
  lower_network[upper.tri(lower_network)] <- 0
  edge_list <- which(lower_network != 0, arr.ind = TRUE)
  
  # Obtain edges to rewire
  rewire_list <- edge_list[sample(
    1:edges, rewire_number, replace = FALSE
  ),]
  
  # Ensure rewire list is matrix
  if(!is.matrix(rewire_list)){
    rewire_list <- matrix(rewire_list, ncol = 2)
  }
  colnames(rewire_list) <- c("row", "col")
  
  # Edge list
  edge_list <- expand.grid.unique(
    1:ncol(network), 1:ncol(network)
  )
  
  # Ensure edge list is matrix
  if(!is.matrix(edge_list)){
    edge_list <- matrix(edge_list, ncol = 2)
  }
  colnames(edge_list) <- c("row", "col")
  
  # Edges to replace
  replace_list <- edge_list[sample(
    1:nrow(edge_list), rewire_number, replace = FALSE
  ),]
  
  # Ensure replace list is matrix
  if(!is.matrix(replace_list)){
    replace_list <- matrix(replace_list, ncol = 2)
  }
  colnames(replace_list) <- c("row", "col")
  
  # Rewire edges
  for(i in 1:nrow(rewire_list)){
    
    # Obtain row and column
    row <- rewire_list[i, "row"]
    column <- rewire_list[i, "col"]
    
    # Obtain replace row and column
    replace_row <- replace_list[i, "row"]
    replace_column <- replace_list[i, "col"]
    
    # Obtain values
    original_value <- network[row, column]
    replace_value <- network[replace_row, replace_column]
    
    # Swap values
    network[replace_row, replace_column] <- original_value
    network[replace_column, replace_row] <- original_value
    network[row, column] <- replace_value
    network[column, row] <- replace_value
    
  }
  
  # Return the rewired network
  return(network)
  
}


#%%%%%%%%%%%%
# dynEGA ----
#%%%%%%%%%%%%

#' @noRd
# Updated 28.08.2022
# GLLA Estimation:
glla.multi <- function(data, n.embed = n.embed, tau = tau, delta = delta, order = order){
  order.deriv <- paste0("Ord",seq(from = 0, to = order))
  data.est <- vector("list")
  for(i in 1:ncol(data)){
    data.est[[i]] <- as.data.frame(EGAnet::glla(data[,i], n.embed = n.embed, tau = tau, delta = delta, order = order))
    data.est[[i]] <- as.data.frame(data.est[[i]])
  }
  data.est2 <- vector("list")
  for(i in 0:order+1){
    data.est2[[i]] <- sapply(data.est, "[[", i)
  }
  data.estimates <- data.frame(Reduce(cbind, data.est2))
  colnames(data.estimates) <- paste(colnames(data), rep(order.deriv, each = ncol(data)), sep = ".")
  return(data.estimates)
}

#%%%%%%%%%%%%%%%%%
# infoCluster ----
#%%%%%%%%%%%%%%%%%

#' @noRd
# Variation of information
# Updated 30.07.2022
vi <- function(wc1, wc2)
{
  
  # Obtain non-NA memberships
  nonNA <- !is.na(wc1) & !is.na(wc2)
  wc1 <- wc1[nonNA]
  wc2 <- wc2[nonNA]
  
  # Compute maximum VI
  ## Set max memberships
  max_wc1 <- rep(1, length(wc1))
  max_wc2 <- 1:length(wc2)
  max_wc2[length(max_wc2)] <- 1
  max_vi <- igraph::compare(
    max_wc1, max_wc2
  )
  
  # Obtain variation of information
  vi <- igraph::compare(
    wc1, wc2, method = "vi"
  )
  
  # Normalize VI by maximum value
  vi <- vi / max_vi
  
  # Return
  return(vi)
  
}

#' @noRd
# Root Mean Square Error (for matrices)
# Updated 30.07.2022
matrix_rmse <- function(matrix1, matrix2)
{
  # Check for symmetric
  if(
    isSymmetric(unname(as.matrix(matrix1))) &
    isSymmetric(unname(as.matrix(matrix2)))
  ){
    
    # Compute lower triangles
    matrix1 <- matrix1[lower.tri(matrix1)]
    matrix2 <- matrix2[lower.tri(matrix2)]
    
  }
  
  # Compute RMSE
  rmse <- sqrt(mean((matrix1 - matrix2)^2, na.rm = TRUE))
  
  # Return RMSE
  return(rmse)
  
}

#' @noRd
# Rescaled Laplacian matrix
# Updated 30.07.2022
rescaled_laplacian <- function(net)
{
  # Ensure diagonal is zero
  diag(net) <- 0
  
  # Make network absolute
  # net <- abs(net)
  # net <- ifelse(net != 0, 1, 0)
  
  # Laplacian matrix
  rescaled_L <- (diag(colSums(net)) - net) / sum(net)
  
  # Return
  return(rescaled_L)
  
}

#' @noRd
# Von Neumann Entropy
# Updated 06.07.2022
vn_entropy <- function(L_mat)
{
  
  # Eigenvalues
  eigenvalues <- eigen(L_mat)$values
  
  # Von Neumann entropy
  vn_entropy <- suppressWarnings(
    -sum(eigenvalues * log2(eigenvalues), na.rm = TRUE)
  )
  
  # Return
  return(vn_entropy)
  
}