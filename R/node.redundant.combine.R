#' Combines Redundant Nodes
#' @description Allows user to combine redundant nodes into
#' sum scores and latent variables to reduce the redundancy of
#' variables in their data
#'
#' @param node.redundant.obj A \code{\link[EGAnet]{node.redundant}} object
#'
#' @param type Character.
#' Method to use to combine redundant variables.
#'
#' \itemize{
#' \item{\code{"sum"}}
#' {Computes sum scores (i.e., means) of the variables}
#'
#' \item{\code{"latent"}}
#' {Computes latent variable scores using \code{[lavaan]{cfa}}}
#' }
#'
#' Defaults to \code{"latent"}
#'
#' @param estimator Character.
#' Estimator to use for latent variables.
#' Defaults to \code{"WLSMV"}.
#' See \code{[lavaan]{cfa}} for more options
#'
#' @param auto NOT RECOMMENDED. Boolean.
#' Should redundant nodes be automatically combined?
#' Defaults to \code{FALSE}.
#' If set to \code{TRUE}, then redundant nodes will combined using
#' the following heuristics:
#'
#' 1. Redundant nodes that form a 3-clique (i.e., a triangle)
#' with the target node are automatically redundant
#'
#' 2. If there are no 3-cliques, then the 2-clique with the
#' largest regularized partial correlation is selected
#'
#' @param ... Options to be passed onto \code{[lavaan]{cfa}}
#'
#' @return Returns a list:
#'
#' \item{data}{New data with redundant variables merged}
#'
#' \item{merged}{A matrix containing the variables that were
#' decided to be redundant with one another}
#'
#' @examples
#' # obtain SAPA items
#' items <- psychTools::spi[,-c(1:10)]
#'
#' \donttest{
#' # weighted topological overlap
#' redund <- node.redundant(items, method = "wTO", type = "adapt")
#'
#' # partial correlation
#' redund <- node.redundant(items, method = "pcor", type = "adapt")
#'
#' # check redundancies
#' key.ind <- match(colnames(items), as.character(psychTools::spi.dictionary$item_id))
#' key <- as.character(psychTools::spi.dictionary$item[key.ind])
#'
#' # change names in redundancy output to questionnaire item description
#' named.nr <- node.redundant.names(redund, key)
#' }
#'
#' if(interactive())
#' {combine <- node.redundant.combine(named.nr)}
#'
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#'
#' @export
# Redundant Nodes Combination Function
# Updated 08.03.2020
node.redundant.combine <- function (node.redundant.obj,
                                    type = c("sum", "latent"),
                                    estimator = "WLSMV",
                                    auto = FALSE, ...)
{
  # Check for node.redundant object class
  if(class(node.redundant.obj) != "node.redundant")
  {stop("A 'node.redundant' object must be used as input")}
  
  # Missing type
  if(missing(type))
  {type <- "latent"
  }else{type <- match.arg(type)}
  
  # Warning for sum scores
  if(type == "sum")
  {
    warning('type = "sum" assumes all variables are reverse scored for a positive manifold',
            immediate. = TRUE)
  }
  
  # Redundant list
  redund <- node.redundant.obj$redundant
  
  # New data
  new.data <- node.redundant.obj$data
  
  # Weights
  net <- as.matrix(node.redundant.obj$network)
  
  # Track merged items
  merged <- list()
  
  # Track changed names
  name.chn <- vector("character")
  
  # Initalize count
  count <- 0
  
  if("key" %in% names(node.redundant.obj))
  {
    key <- node.redundant.obj$key
    names(key) <- names(node.redundant.obj$key)
  }else{
    key <- colnames(node.redundant.obj$data)
    names(key) <- key
  }
  
  # Loop through named node redundant list
  while(length(redund) != 0)
  {
    # Tracking list
    track <- redund
    
    # Targeting redundancy
    target.item <- names(redund)[1]
    
    # Potential redundancies
    pot <- redund[[1]]
    
    if(length(pot) != 0)
    {
      # Configure into list
      pot <- list(pot)
      names(pot) <- target.item
      
      # Initialize while escape
      escape <- FALSE
      
      # Initialize count
      count2 <- 1
      
      while(!escape)
      {
        # Identify target redundancies
        target <- redund[na.omit(match(pot[[count2]], names(redund)))]
        
        if(length(target) == 0)
        {
          # Escape while loop
          escape <- TRUE
        }else{
          # Increase count
          count2 <- count2 + 1
          
          # Input extended redundancies
          pot[[count2]] <- target
        }
      }
      
      # Possible options
      poss <- unique(unname(unlist(pot)))
      
      # Organize plot of redundancy connections
      mat <- matrix(0, nrow = length(poss) + 1, ncol = length(poss) + 1)
      colnames(mat) <- c(paste("Target"), 1:length(poss))
      row.names(mat) <- colnames(mat)
      
      mat["Target",paste(1:length(unlist(pot[[1]])))] <- net[names(key[match(target.item, key)]),names(key[match(unlist(pot[[1]]),key)])]
      mat[paste(1:length(unlist(pot[[1]]))),"Target"] <- net[names(key[match(target.item, key)]),names(key[match(unlist(pot[[1]]),key)])]
      
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
            elem <- match(names(target.ext)[j], poss)
            
            # Get elements redundant with element
            red.elem <- match(single, poss)
            
            # Put into matrix
            mat[paste(elem),paste(red.elem)] <- net[names(key[match(poss[elem], key)]),names(key[match(poss[red.elem],key)])]
            mat[paste(red.elem),paste(elem)] <- net[names(key[match(poss[elem], key)]),names(key[match(poss[red.elem],key)])]
          }
        }
      }
      
      if(!auto)
      {
        # Print target and potential options
        cat(paste("Target item: '", target.item, "'", sep = ""))
        cat("\n\nPotential redundancies:\n\n")
        cat("0. Do not combine with any")
        cat(paste("\n", 1:length(poss), ". ", "'", poss, "'", sep = ""),"\n")
        
        # Plot
        qgraph::qgraph(mat, color = c(paste("red"),rep("white",length(poss))), layout = "circle",
                       title = "Regularized partial correlations")
        
        # Get input
        input <- readline(prompt = "Enter numbers of items redundant with the target item (separate by commas): ")
        
        # Input check function
        in.check <- function(input)
        {
          inp <- suppressWarnings(as.numeric(unlist(strsplit(unlist(strsplit(input, split = " ")), split = ","))))
          
          ret.val <- FALSE
          
          if(any(is.na(inp)))
          {ret.val <- TRUE}
          
          if(length(inp) == 0)
          {ret.val <- TRUE}
          
          if(length(setdiff(inp, 0:length(poss))) != 0)
          {ret.val <- TRUE}
          
          return(ret.val)
        }
        
        # Redo input check
        re.input <- in.check(input)
        
        while(re.input)
        {
          # Print message to try again
          message("Inappropriate input. Try again.")
          
          # Get input
          input <- readline(prompt = "Enter numbers of items redundant with the target item (separate by commas): ")
          
          # Redo input check
          re.input <- in.check(input)
        }
        
      }else{
        
        # Get cliques
        cliq <- igraph::max_cliques(NetworkToolbox::convert2igraph(mat))
        
        # Initialize count 3
        count3 <- 0
        
        # Remove cliques not containing one
        for(n in 1:length(cliq))
        {
          # Increase count 3
          count3 <- count3 + 1
          
          if(!any(cliq[[count3]] == 1))
          {
            # Remove clique
            cliq[[count3]] <- NULL
            
            # Decrease count 3
            count3 <- count3 - 1
          }else{
            cliq[[count3]] <- cliq[[count3]] - 1
          }
        }
        
        # Clique lengths
        cliq.len <- unlist(lapply(cliq, length))
        
        # Get maximum clique size
        max.cliq <- max(cliq.len)
        
        # Heuristics
        if(max.cliq >= 3)
        {input <- paste(setdiff(unique(unlist(cliq[which(cliq.len == max.cliq)])),0), collapse = ", ")
        }else{
          
          # Obtain cliques
          cliq <- cliq[which(cliq.len == max.cliq)]
          
          # Weights vector
          wei <- numeric(length(cliq))
          
          # Identify maximum weight
          for(n in 1:length(cliq))
          {wei[n] <- net[names(key[match(poss[cliq[[n]]], key)]),names(key[match(target.item,key)])]}
          
          # Obtain node for input
          input <- paste(setdiff(cliq[[which.max(wei)]],0))
        }
      }
      
      if(input != "0")
      {
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
        if(type == "sum")
        {new.data[,tar.idx] <- rowMeans(new.data[,c(tar.idx, idx)], na.rm = TRUE)
        }else if(type == "latent")
        {
          # Latent variable
          ## create model
          mod <- paste(paste("comb =~ ",sep=""), paste(colnames(new.data[,c(tar.idx, idx)]), collapse = " + "))
          
          ## fit model
          fit <- suppressWarnings(lavaan::cfa(mod, data = new.data, ...))
          
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
            message("Some items are reverse coded (negative correlations with latent variable were found). Correlations with latent variable:")
            print(corrs)
            
            input2 <- "o"
            
            while(input2 != "y" && input2 != "n")
            {input2 <- readline("Reverse code for positive labelling (y/n): ")}
            
            if(input2 == "y")
            {new.vec <- -new.vec}
          }
          
          # input new vector
          new.data[,tar.idx] <- new.vec
        }
        
        # Ask for new label
        if(!auto)
        {lab <- readline(prompt = "New label for item (no quotations): ")
        }else{lab <- paste("cvar",count,sep="")}
        name.chn[count] <- lab
        col.idx <- match(tar.idx, colnames(new.data))
        colnames(new.data)[col.idx] <- lab
        
        # Remove redundant variables from data
        rm.idx <- match(idx, colnames(new.data))
        new.data <- new.data[,-rm.idx]
        
        # Remove variables from potential future options
        opts <- redund[na.omit(match(comb, names(redund)))]
        
        if(length(opts) != 0)
        {redund[names(opts)] <- NULL}
        
        # Remove target item
        redund[[1]] <- NULL
        
        # Remove variables within future options
        rm.var <- which(lapply(lapply(redund, match, comb), function(x){any(!is.na(x))}) == TRUE)
        
        if(length(rm.var) != 0)
        {
          for(j in 1:length(rm.var))
          {
            # Target option
            target.opt <- redund[rm.var][[j]]
            
            # Remove target variable(s)
            target.var <- na.omit(match(comb, target.opt))
            
            redund[rm.var][[j]] <- target.opt[-target.var]
          }
        }
        
      }else{
        # Map target item to column names of new data
        item.name <- names(key)[match(target.item, key)]
        target.col <- match(item.name, colnames(new.data))
        colnames(new.data)[target.col] <- target.item
        redund[[1]] <- NULL
      }
      
    }else{
      # Map target item to column names of new data
      item.name <- names(key)[match(target.item, key)]
      target.col <- match(item.name, colnames(new.data))
      colnames(new.data)[target.col] <- target.item
      redund[[1]] <- NULL
    }
    
  }
  
  # Transform merged list to matrix
  if(length(merged) != 0)
  {
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
  if(any(colnames(new.data) %in% names(key)))
  {
    # Target names
    target.names <- which(colnames(new.data) %in% names(key))
    
    # new.data names
    new.data.names <- colnames(new.data)[target.names]
    
    # Insert into new data
    colnames(new.data)[target.names] <- key[new.data.names]
  }
  
  # Check if 'm.mat' exists
  if(!exists("m.mat"))
  {m.mat <- NULL}
  
  # Initialize results list
  res <- list()
  res$data <- new.data
  res$merged <- m.mat
  
  return(res)
  
}
#----
