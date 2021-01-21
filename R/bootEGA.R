#' Dimension Stability Analysis of \code{\link[EGAnet]{EGA}}
#'
#' \code{bootEGA} Estimates the number of dimensions of \emph{n} bootstraps
#' using the empirical (partial) correlation matrix (parametric) or resampling from
#' the empirical dataset (non-parametric). It also estimates a typical
#' median network structure, which is formed by the median or mean pairwise (partial)
#' correlations over the \emph{n} bootstraps.
#'
#' @param data Matrix or data frame.
#' Includes the variables to be used in the \code{bootEGA} analysis
#'
#' @param uni Boolean.
#' Should unidimensionality be checked?
#' Defaults to \code{TRUE}.
#' Set to \code{FALSE} to check for multidimensionality only.
#' If \code{TRUE}, then the same number of variables as the original
#' data (i.e., from argument \code{data}) are generated from a factor
#' model with one factor and loadings of .70. These data are then
#' appended to the original data and dimensionality is checked.
#' If the number of dimensions is one or two, then the original
#' data are unidimensional; otherwise, the data are multidimensional
#' (see Golino, Shi, et al., 2020 for more details)
#'
#' @param iter Numeric integer.
#' Number of replica samples to generate from the bootstrap analysis.
#' At least \code{500} is recommended
#'
#' @param type Character.
#' A string indicating the type of bootstrap to use.
#'
#' Current options are:
#'
#' \itemize{
#'
#' \item{\strong{\code{"parametric"}}}
#' {Generates \code{n} new datasets (multivariate normal random distributions) based on the
#' original dataset, via the \code{\link[MASS]{mvrnorm}} function}
#'
#' \item{\strong{\code{"resampling"}}}
#' {Generates n random subsamples of the original data}
#'
#' }
#'
#' @param model Character.
#' A string indicating the method to use.
#'
#' Current options are:
#'
#' \itemize{
#'
#' \item{\strong{\code{glasso}}}
#' {Estimates the Gaussian graphical model using graphical LASSO with
#' extended Bayesian information criterion to select optimal regularization parameter.
#' This is the default method}
#'
#' \item{\strong{\code{TMFG}}}
#' {Estimates a Triangulated Maximally Filtered Graph}
#'
#' }
#'
#' @param model.args List.
#' A list of additional arguments for \code{\link[EGAnet]{EBICglasso.qgraph}}
#' or \code{\link[NetworkToolbox]{TMFG}}
#'
#' @param algorithm A string indicating the algorithm to use or a function from \code{\link{igraph}}
#' Current options are:
#'
#' \itemize{
#'
#' \item{\strong{\code{walktrap}}}
#' {Computes the Walktrap algorithm using \code{\link[igraph]{cluster_walktrap}}}
#'
#' \item{\strong{\code{louvain}}}
#' {Computes the Walktrap algorithm using \code{\link[igraph]{cluster_louvain}}}
#'
#' }
#'
#' @param algorithm.args List.
#' A list of additional arguments for \code{\link[igraph]{cluster_walktrap}}, \code{\link[igraph]{cluster_louvain}},
#' or some other community detection algorithm function (see examples)
#'
#' @param typicalStructure Boolean.
#' If \code{TRUE}, returns the typical network of partial correlations
#' (estimated via graphical lasso or via TMFG) and estimates its dimensions.
#' The "typical network" is the median of all pairwise correlations over the \emph{n} bootstraps.
#' Defaults to \code{TRUE}
#'
#' @param plot.typicalStructure Boolean.
#' If \code{TRUE}, returns a plot of the typical network (partial correlations),
#' which is the median of all pairwise correlations over the \emph{n} bootstraps,
#' and its estimated dimensions.
#' Defaults to \code{TRUE}
#'
#' @param plot.type Character.
#' Plot system to use.
#' Current options are \code{\link[qgraph]{qgraph}} and \code{\link{GGally}}.
#' Defaults to \code{"GGally"}.
#'
#' @param plot.args List.
#' A list of additional arguments for the network plot.
#' For \code{plot.type = "qgraph"}:
#'
#' \itemize{
#'
#' \item{\strong{\code{vsize}}}
#' {Size of the nodes. Defaults to 6.}
#'
#'}
#' For \code{plot.type = "GGally"} (see \code{\link[GGally]{ggnet2}} for
#' full list of arguments):
#'
#' \itemize{
#'
#' \item{\strong{\code{vsize}}}
#' {Size of the nodes. Defaults to 6.}
#'
#' \item{\strong{\code{label.size}}}
#' {Size of the labels. Defaults to 5.}
#'
#' \item{\strong{\code{alpha}}}
#' {The level of transparency of the nodes, which might be a single value or a vector of values. Defaults to 0.7.}
#'
#' \item{\strong{\code{edge.alpha}}}
#' {The level of transparency of the edges, which might be a single value or a vector of values. Defaults to 0.4.}
#' 
#'  \item{\strong{\code{legend.names}}}
#' {A vector with names for each dimension}
#' 
#' \item{\strong{\code{color.palette}}}
#' {The color palette for the nodes. For custom colors,
#' enter HEX codes for each dimension in a vector.
#' See \code{\link[EGAnet]{color_palette_EGA}} for 
#' more details and examples}
#' 
#' }
#'
#'
#' @param ncores Numeric.
#' Number of cores to use in computing results.
#' Defaults to \code{parallel::detectCores() / 2} or half of your
#' computer's processing power.
#' Set to \code{1} to not use parallel computing
#'
#' If you're unsure how many cores your computer has,
#' then use the following code: \code{parallel::detectCores()}
#'
#' @param ... Additional arguments.
#' Used for deprecated arguments from previous versions of \code{\link{EGA}}
#'
#' @return Returns a list containing:
#'
#' \item{iter}{Number of replica samples in bootstrap}
#'
#' \item{boot.ndim}{Number of dimensions identified in each replica sample}
#'
#' \item{boot.wc}{Item allocation for each replica sample}
#'
#' \item{bootGraphs}{Networks of each replica sample}
#'
#' \item{summary.table}{Summary table containing number of replica samples, median,
#' standard deviation, standard error, 95\% confidence intervals, and quantiles (lower = 2.5\% and upper = 97.5\%)}
#'
#' \item{frequency}{Proportion of times the number of dimensions was identified
#' (e.g., .85 of 1,000 = 850 times that specific number of dimensions was found)}
#'
#' \item{EGA}{Output of the original \code{\link[EGAnet]{EGA}} results}
#'
#' \item{typicalGraph}{A list containing:
#'
#' \itemize{
#'
#' \item{\strong{\code{graph}}}
#' {Network matrix of the median network structure}
#'
#' \item{\strong{\code{typical.dim.variables}}}
#' {An ordered matrix of item allocation}
#'
#' \item{\strong{\code{wc}}}
#' {Item allocation of the median network}
#'    }
#' }
#'
#' @author Hudson Golino <hfg9s at virginia.edu> and Alexander P. Christensen <alexpaulchristensen@gmail.com>
#'
#' @examples
#' # Load data
#' wmt <- wmt2[,7:24]
#'
#' \donttest{# bootEGA glasso example
#' ## plot.type = "qqraph" used for CRAN checks
#' ## plot.type = "GGally" is the default
#' boot.wmt <- bootEGA(data = wmt, iter = 500, plot.type = "qgraph",
#' type = "parametric", ncores = 2)
#'
#' # bootEGA TMFG example
#' boot.wmt <- bootEGA(data = wmt, iter = 500, model = "TMFG",
#' plot.type = "qgraph", type = "parametric", ncores = 2)
#'
#' # bootEGA Louvain example
#' boot.wmt <- bootEGA(data = wmt, iter = 500, algorithm = "louvain",
#' plot.type = "qgraph", type = "parametric", ncores = 2)
#'
#' # bootEGA Spinglass example
#' boot.wmt <- bootEGA(data = wmt, iter = 500, model = "TMFG", plot.type = "qgraph",
#' algorithm = igraph::cluster_spinglass, type = "parametric", ncores = 2)
#' }
#'
#' # Load data
#' intwl <- intelligenceBattery[,8:66]
#'
#' \donttest{# Another bootEGA example
#' boot.intwl <- bootEGA(data = intwl, iter = 500,
#' plot.type = "qgraph", type = "parametric", ncores = 2)
#' }
#'
#' @references
#' # Original implementation of bootEGA \cr
#' Christensen, A. P., & Golino, H. (2019).
#' Estimating the stability of the number of factors via Bootstrap Exploratory Graph Analysis: A tutorial.
#' \emph{PsyArXiv}.
#' doi:\href{https://doi.org/10.31234/osf.io/9deay}{10.31234/osf.io/9deay}
#'
#' # Structural consistency (see \code{\link[EGAnet]{dimStability}}) \cr
#' Christensen, A. P., Golino, H., & Silvia, P. J. (in press).
#' A psychometric network perspective on the validity and validation of personality trait questionnaires.
#' \emph{European Journal of Personality}.
#' doi: \href{https://doi.org/10.1002/per.2265}{10.1002/per.2265}
#'
#' @seealso \code{\link[EGAnet]{EGA}} to estimate the number of dimensions of an instrument using EGA
#' and \code{\link[EGAnet]{CFA}} to verify the fit of the structure suggested by EGA using confirmatory factor analysis.
#'
#' @importFrom stats cov median sd qt quantile
#'
#' @export
#'
# Bootstrap EGA
# Updated 20.01.2021
bootEGA <- function(data, uni = TRUE, iter, type = c("parametric", "resampling"),
                    model = c("glasso", "TMFG"), model.args = list(),
                    algorithm = c("walktrap", "louvain"), algorithm.args = list(),
                    typicalStructure = TRUE, plot.typicalStructure = TRUE,
                    plot.type = c("GGally", "qgraph"),
                    plot.args = list(), ncores, ...) {

  #### DEPRECATED ARGUMENTS ####

  # Get additional arguments
  add.args <- list(...)

  # Check if n has been input as an argument
  if("n" %in% names(add.args)){

    # Give deprecation warning
    warning(
      paste(
        "The 'n' argument has been deprecated in the bootEGA function.\n\nInstead use: iter = ", add.args$n,
        sep = ""
      )
    )

    # Handle the number of iterations appropriately
    iter <- add.args$n
  }

  # Check if steps has been input as an argument
  if("steps" %in% names(add.args)){

    # Give deprecation warning
    warning(
      paste(
        "The 'steps' argument has been deprecated in all EGA functions.\n\nInstead use: algorithm.args = list(steps = ", add.args$steps, ")",
        sep = ""
      )
    )

    # Handle the number of steps appropriately
    algorithm.args$steps <- add.args$steps
  }

  #### DEPRECATED ARGUMENTS ####

  #### MISSING ARGUMENTS HANDLING ####

  if(missing(model)){
    model <- "glasso"
  }else{
    model <- match.arg(model)
  }

  if(missing(algorithm)){
    algorithm <- "walktrap"
  }else if(!is.function(algorithm)){
    algorithm <- match.arg(algorithm)
  }

  if(missing(type)){
    type <- "parametric"
  }else{
    type <- match.arg(type)
  }

  if(missing(ncores)){
    ncores <- ceiling(parallel::detectCores() / 2)
  }

  if(missing(plot.type)){
    plot.type <- "GGally"
  }else{
    plot.type <- match.arg(plot.type)
  }

  ## Check for input plot arguments
  color.palette <- "polychrome"
  if(plot.type == "GGally"){
    
    if(length(plot.args) == 0){
      
      default.args <- formals(GGally::ggnet2)
      ega.default.args <- list(size = 6, alpha = 0.4, label.size = 5,
                               edge.alpha = 0.7, layout.exp = 0.2)
      default.args[names(ega.default.args)]  <- ega.default.args
      default.args <- default.args[-length(default.args)]
      
    }else{
      
      default.args <- formals(GGally::ggnet2)
      ega.default.args <- list(size = 6, alpha = 0.4, label.size = 5,
                               edge.alpha = 0.7, layout.exp = 0.2)
      default.args[names(ega.default.args)]  <- ega.default.args
      default.args <- default.args[-length(default.args)]
      
      
      if("vsize" %in% names(plot.args)){
        plot.args$size <- plot.args$vsize
        plot.args$vsize <- NULL
      }
      
      if("color.palette" %in% names(plot.args)){
        color.palette <- plot.args$color.palette
      }
      
      if(any(names(plot.args) %in% names(default.args))){
        target.args <- plot.args[which(names(plot.args) %in% names(default.args))]
        default.args[names(target.args)] <- target.args
      }
      
    }
    
    plot.args <- default.args
    
  }

  #### MISSING ARGUMENTS HANDLING ####

  #number of cases
  cases <- nrow(data)

  #set inverse covariance matrix for parametric approach
  if(type=="parametric"){  # Use a parametric approach

    if(model=="glasso"){

      g <- -suppressMessages(EGA.estimate(data = data, n = cases, model = model, model.args = model.args)$network)
      diag(g) <- 1

    }else if(model=="TMFG"){

      g <- -suppressMessages(NetworkToolbox::LoGo(data, normal = TRUE, partial = TRUE))
      diag(g) <- 1

    }
  }

  #initialize data list
  datalist <- list()

  #initialize count
  count <- 0

  #let user know data generation has started
  message("\nGenerating data...", appendLF = FALSE)

  repeat{

    #increase count
    count <- count + 1

    #generate data
    if(type == "parametric"){

      datalist[[count]] <- MASS::mvrnorm(cases, mu = rep(0, ncol(g)), Sigma = corpcor::pseudoinverse(g))

    }else if(type == "resampling"){

      datalist[[count]] <- data[sample(1:cases, replace=TRUE),]

    }

    #break out of repeat
    if(count == iter)
    {break}
  }

  #let user know data generation has ended
  message("done", appendLF = TRUE)

  #Parallel processing
  cl <- parallel::makeCluster(ncores)

  #Export variables
  parallel::clusterExport(cl = cl,
                          varlist = c("datalist", "uni", "cases",
                                      "model", "model.args",
                                      "algorithm", "algorithm.args"),
                          envir=environment())

  #let user know data generation has started
  message("Estimating EGA networks...\n", appendLF = FALSE)

  #Estimate networks
  boots <- pbapply::pblapply(
    X = datalist, cl = cl,
    FUN = EGA,
    uni = uni,
    model = model, model.args = model.args,
    algorithm = algorithm, algorith.args = algorithm.args,
    plot.EGA = FALSE
  )

  parallel::stopCluster(cl)

  #let user know results are being computed
  message("Computing results...\n")

  #get networks
  bootGraphs <- lapply(boots, function(x, col.names){
    net <- x$network
    colnames(net) <- col.names
    row.names(net) <- col.names
    return(net)
  }, col.names = colnames(data))

  #get community membership
  boot.wc <- lapply(boots, function(x, col.names){
    wc <- x$wc
    names(wc) <- col.names
    return(wc)
  }, col.names = colnames(data))

  #get dimensions
  boot.ndim <- matrix(NA, nrow = iter, ncol = 2)
  colnames(boot.ndim) <- c("Boot.Number", "N.Dim")

  boot.ndim[,1] <- seq_len(iter)
  boot.ndim[,2] <- unlist(
    lapply(boots, function(x){
      x$n.dim
    })
  )

  if (typicalStructure){

    typical.Structure <- switch(model,
                                glasso = apply(simplify2array(bootGraphs),1:2, median),
                                TMFG = apply(simplify2array(bootGraphs),1:2, mean)
                         )

    # Sub-routine to following EGA approach (handles undimensional structures)
    typical.wc <- typicalStructure.network(A = typical.Structure,
                                           model = model, model.args = model.args,
                                           n = cases, uni = uni, algorithm = algorithm,
                                           algorithm.args = algorithm.args)

    typical.ndim <- max(typical.wc, na.rm = TRUE)
    dim.variables <- data.frame(items = colnames(data), dimension = typical.wc)
  }
  if (plot.typicalStructure) {
    if(plot.type == "qgraph"){
      plot.typical.ega <- qgraph::qgraph(typical.Structure, layout = "spring",
                                         vsize = plot.args$vsize, groups = as.factor(typical.wc))
    }else if(plot.type == "GGally"){
        network1 <- network::network(typical.Structure,
                                     ignore.eval = FALSE,
                                     names.eval = "weights",
                                     directed = FALSE)

      network::set.vertex.attribute(network1, attrname= "Communities", value = typical.wc)
      network::set.vertex.attribute(network1, attrname= "Names", value = network::network.vertex.names(network1))
      network::set.edge.attribute(network1, "color", ifelse( network::get.edge.value(network1, "weights") > 0, "darkgreen", "red"))
      network::set.edge.value(network1,attrname="AbsWeights",value=abs(typical.Structure))
      network::set.edge.value(network1,attrname="ScaledWeights",
                              value=matrix(#scales::rescale(typical.Structure),
                                rescale.edges(typical.Structure, plot.args$size),
                                           nrow = nrow(typical.Structure),
                                           ncol = ncol(typical.Structure)))

      # Layout "Spring"
      graph1 <- igraph::as.igraph(qgraph::qgraph(typical.Structure, DoNotPlot = TRUE))
      edge.list <- igraph::as_edgelist(graph1)
      layout.spring <- qgraph::qgraph.layout.fruchtermanreingold(edgelist = edge.list,
                                                                 weights =
                                                                   abs(igraph::E(graph1)$weight/max(abs(igraph::E(graph1)$weight)))^2,
                                                                 vcount = ncol(typical.Structure))


      set.seed(1234)
      plot.args$net <- network1
      plot.args$node.color <- "Communities"
      plot.args$node.alpha <- plot.args$alpha
      plot.args$node.shape <- plot.args$shape
      plot.args$node.size <- plot.args$size
      plot.args$edge.color <- "color"
      plot.args$edge.size <- "ScaledWeights"
      plot.args$color.palette <- "Set1"
      
      lower <- abs(typical.Structure[lower.tri(typical.Structure)])
      non.zero <- sqrt(lower[lower != 0])
      
      plot.args$edge.alpha <- non.zero
      plot.args$mode <- layout.spring
      plot.args$label <- colnames(typical.Structure)
      plot.args$node.label <- plot.args$label
      if(plot.args$label.size == "max_size/2"){plot.args$label.size <- plot.args$size/2}
      if(plot.args$edge.label.size == "max_size/2"){plot.args$edge.label.size <- plot.args$size/2}
      
      plot.typical.ega <- do.call(GGally::ggnet2, plot.args) + ggplot2::theme(legend.title = ggplot2::element_blank())
      
      plot.typical.ega <- suppressMessages(
        do.call(GGally::ggnet2, plot.args) + 
          ggplot2::theme(legend.title = ggplot2::element_blank()) +
          ggplot2::scale_color_manual(values = color_palette_EGA(color.palette, typical.wc),
                                      breaks = sort(typical.wc)) +
          ggplot2::guides(
            color = ggplot2::guide_legend(override.aes = list(
              size = plot.args$size,
              alpha = plot.args$alpha
            ))
          )
      )
      
      plot(plot.typical.ega)
    }
    
    set.seed(NULL)


  }
  Median <- median(boot.ndim[, 2])
  se.boot <- sd(boot.ndim[, 2])
  ciMult <- qt(0.95/2 + 0.5, nrow(boot.ndim) - 1)
  ci <- se.boot * ciMult
  quant <- quantile(boot.ndim[,2], c(.025, .975), na.rm = TRUE)
  summary.table <- data.frame(n.Boots = iter, median.dim = Median,
                              SE.dim = se.boot, CI.dim = ci,
                              Lower.CI = Median - ci, Upper.CI = Median + ci,
                              Lower.Quantile = quant[1], Upper.Quantile = quant[2])
  row.names(summary.table) <- NULL

  #compute frequency
  dim.range <- range(boot.ndim[,2])
  lik <- matrix(0, nrow = diff(dim.range)+1, ncol = 2)
  colnames(lik) <- c("# of Factors", "Frequency")
  count <- 0

  for(i in seq(from=min(dim.range),to=max(dim.range),by=1))
  {
    count <- count + 1
    lik[count,1] <- i
    lik[count,2] <- length(which(boot.ndim[,2]==i))/iter
  }

  result <- list()
  result$iter <- iter
  result$type <- type
  result$boot.ndim <- boot.ndim
  result$boot.wc <- boot.wc
  result$bootGraphs <- bootGraphs
  result$summary.table <- summary.table
  result$frequency <- lik
  result$EGA <- suppressMessages(suppressWarnings(EGA(data = data, uni = uni,
                                                      model = model, model.args = model.args,
                                                      algorithm = algorithm, algorith.args = algorithm.args,
                                                      plot.EGA = FALSE)))

  # Typical structure
  if (typicalStructure) {
    typicalGraph <- list()
    typicalGraph$graph <- typical.Structure
    typicalGraph$typical.dim.variables <- dim.variables[order(dim.variables[,2]), ]
    typicalGraph$wc <- typical.wc
    result$typicalGraph <- typicalGraph
    if(plot.typicalStructure){
      result$plot.typical.ega <- plot.typical.ega
    }
  }
  
  # Add plot arguments (for itemStability)
  result$color.palette <- color.palette

  class(result) <- "bootEGA"

  # Message that unidimensional structures were not checked
  if(!uni){
    message("\nEGA did not check for unidimensionality. Set argument 'uni' to TRUE to check for unidimensionality")
  }

  return(result)
}
#----
