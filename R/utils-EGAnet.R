#%%%%%%%%%%%%%%%%%%%%
# NETWORKTOOLBOX ----
#%%%%%%%%%%%%%%%%%%%%

# pwr.r.test
#' @noRd
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

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# MULTI-FUNCTION SUB-ROUTINES ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

#%%%%%%%%%%%%
# oldUVA ----
#%%%%%%%%%%%%

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
    
    tom <- wto(net)
    
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
    
    check_package("fitdistrplus")
    
    # Distributions, initialize AIC vector
    distr <- c("norm", "gamma")
    aic <- numeric(length(distr))
    names(aic) <- c("normal", "gamma")
    
    ## Obtain distribution
    for(i in 1:length(distr)){
      capture.output(
        aic[i] <- fitdistrplus::fitdist(pos.vals, distr[i], method="mle")$aic
      )
    }
    
    ## Obtain parameters
    g.dist <- suppressWarnings(
      fitdistrplus::fitdist(
        pos.vals, distr[which.min(aic)]
      )$estimate
    )
    
    # Estimate p-values
    pval <- switch(names(aic)[which.min(aic)],
                   
                   normal = 1 - unlist(lapply(pos.vals, # positive values
                                              pnorm, # probability in normal distribution
                                              mean = g.dist["mean"], #mean of normal
                                              sd = g.dist["sd"]) #standard deviation of normal
                   ),
                   
                   gamma = 1 - unlist(lapply(pos.vals, # positive values
                                             pgamma, # probability in gamma distribution
                                             shape = g.dist["shape"], # shape of gamma
                                             rate = g.dist["rate"]) # rate of gamma
                   ),
    )
    
    # Check if using adaptive alpha
    if(type == "adapt"){
      sig <- adapt.a("cor", alpha = sig, n = length(pos.vals), efxize = "medium")$adapt.a
    }
    
    # Get redundant pairings
    redund <- pos.vals[which(pval <= sig)]
    
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
  
  # Check for {gridExtra}
  check_package("gridExtra")
  
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
# Updated 06.05.2022
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