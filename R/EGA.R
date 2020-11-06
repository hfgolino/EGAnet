#' Applies the Exploratory Graph Analysis technique
#'
#' Estimates the number of dimensions of a given dataset or correlation matrix
#' using the graphical lasso (\code{\link{EBICglasso.qgraph}}) or the
#' Triangulated Maximally Filtered Graph (\code{\link[NetworkToolbox]{TMFG}})
#' network estimation methods.
#'
#' Two community detection algorithms, Walktrap (Pons & Latapy, 2006) and
#' Louvain (Blondel et al., 2008), are pre-programmed because of their
#' superior performance in simulation studies on psychological
#' data generated from factor models (Christensen & Golino; 2020; Golino et al., 2020).
#' Notably, any community detection algorithm from the \code{\link{igraph}}
#' can be used to estimate the number of communities (see examples).
#'
#' @param data Matrix or data frame.
#' Variables (down columns) or correlation matrix.
#' If the input is a correlation matrix,
#' then argument \code{n} (number of cases) is \strong{required}
#'
#' @param n Integer.
#' Sample size if \code{data} provided is a correlation matrix
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
#' @param model Character.
#' A string indicating the method to use.
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
#'
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
#' @param plot.EGA Boolean.
#' If \code{TRUE}, returns a plot of the network and its estimated dimensions.
#' Defaults to \code{TRUE}
#'
#' @param plot.type Character.
#' Plot system to use.
#' Current options are \code{\link[qgraph]{qgraph}} and \code{\link{GGally}}.
#' Defaults to \code{"GGally"}
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
#' For \code{plot.type = "GGally"}:
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
#' {The level of transparency of the nodes, which might be a single value or a vector of values. Defaults to 0.4.}
#'
#' \item{\strong{\code{edge.alpha}}}
#' {The level of transparency of the edges, which might be a single value or a vector of values. Defaults to 0.7.}
#' }
#'
#' @param ... Additional arguments.
#' Used for deprecated arguments from previous versions of \code{\link{EGA}}
#'
#' @author Hudson Golino <hfg9s at virginia.edu>, Alexander P. Christensen <alexpaulchristensen at gmail.com>, Maria Dolores Nieto <acinodam at gmail.com> and Luis E. Garrido <garrido.luiseduardo at gmail.com>
#'
#' @return Returns a list containing:
#'
#' \item{network}{A symmetric network estimated using either the
#' \code{\link{EBICglasso.qgraph}} or \code{\link[NetworkToolbox]{TMFG}}}
#'
#' \item{wc}{A vector representing the community (dimension) membership
#' of each node in the network. \code{NA} values mean that the node
#' was disconnected from the network}
#'
#' \item{n.dim}{A scalar of how many total dimensions were identified in the network}
#'
#' \item{cor.data}{The zero-order correlation matrix}
#'
#' \item{Methods}{Arguments for creating a Methods section (see \code{\link[EGAnet]{EGA.methods.section}})}
#'
#' @examples
#' \donttest{# Estimate EGA
#' ## plot.type = "qqraph" used for CRAN checks
#' ## plot.type = "GGally" is the default
#' ega.wmt <- EGA(data = wmt2[,7:24], plot.type = "qgraph")
#'
#' # Summary statistics
#' summary(ega.wmt)
#'
#' # Estimate EGAtmfg
#' ega.wmt <- EGA(data = wmt2[,7:24], model = "TMFG", plot.type = "qgraph")
#'
#' # Estimate EGA with Louvain algorithm
#' ega.wmt <- EGA(data = wmt2[,7:24], algorithm = "louvain", plot.type = "qgraph")
#'
#' # Estimate EGA with Spinglass algorithm
#' ega.wmt <- EGA(data = wmt2[,7:24],
#' algorithm = igraph::cluster_spinglass, plot.type = "qgraph")
#'
#' # Estimate EGA
#' ega.intel <- EGA(data = intelligenceBattery[,8:66], model = "glasso", plot.EGA = FALSE)
#'
#' # Summary statistics
#' summary(ega.intel)
#' }
#'
#'  \dontshow{# Fast for CRAN checks
#' # Pearson's correlation matrix
#' wmt <- cor(wmt2[,7:24])
#'
#' # Estimate EGA
#' ega.wmt <- EGA(data = wmt, n = nrow(wmt2), uni = FALSE, model = "glasso", plot.EGA = FALSE)
#'
#' }
#'
#' @seealso \code{\link{bootEGA}} to investigate the stability of EGA's estimation via bootstrap
#' and \code{\link{CFA}} to verify the fit of the structure suggested by EGA using confirmatory factor analysis.
#'
#' @references
#' # Louvain algorithm \cr
#' Blondel, V. D., Guillaume, J.-L., Lambiotte, R., & Lefebvre, E. (2008).
#' Fast unfolding of communities in large networks.
#' \emph{Journal of Statistical Mechanics: Theory and Experiment}, \emph{2008}, P10008.
#' doi: \href{https://doi.org/10.1088/1742-5468/2008/10/P10008}{10.1088/1742-5468/2008/10/P10008}
#'
#' # Compared all \emph{igraph} community detections algorithms, introduced Louvain algorithm, simulation with continuous and polytomous data \cr
#' Christensen, A. P., & Golino, H. (under review).
#' Estimating factors with psychometric networks: A Monte Carlo simulation comparing community detection algorithms.
#' \emph{PsyArXiv}.
#' doi: \href{https://doi.org/10.31234/osf.io/hz89e}{10.31234/osf.io/hz89e}
#'
#' # Original simulation and implementation of EGA \cr
#' Golino, H. F., & Epskamp, S. (2017).
#' Exploratory graph analysis: A new approach for estimating the number of dimensions in psychological research.
#' \emph{PloS one}, \emph{12(6)}, e0174035..
#' doi: \href{https://doi.org/10.1371/journal.pone.0174035}{journal.pone.0174035}
#'
#' Golino, H. F., & Demetriou, A. (2017).
#' Estimating the dimensionality of intelligence like data using Exploratory Graph Analysis.
#' \emph{Intelligence}, \emph{62}, 54-70.
#' doi: \href{https://doi.org/10.1016/j.intell.2017.02.007}{j.intell.2017.02.007}
#'
#' # Current implementation of EGA, introduced unidimensional checks, continuous and dichotomous data \cr
#' Golino, H., Shi, D., Christensen, A. P., Garrido, L. E., Nieto, M. D., Sadana, R., & Thiyagarajan, J. A. (2020).
#' Investigating the performance of Exploratory Graph Analysis and traditional techniques to identify the number of latent factors: A simulation and tutorial.
#' \emph{Psychological Methods}, \emph{25}, 292-320.
#' doi: \href{https://doi.org/10.1037/met0000255}{10.1037/met0000255}
#'
#' # Walktrap algorithm \cr
#' Pons, P., & Latapy, M. (2006).
#' Computing communities in large networks using random walks.
#' \emph{Journal of Graph Algorithms and Applications}, \emph{10}, 191-218.
#' doi: \href{https://doi.org/10.7155/jgaa.00185}{10.7155/jgaa.00185}
#'
#' @importFrom stats cor rnorm runif na.omit
#'
#' @export
#'
# Updated 30.10.2020
## EGA Function to detect unidimensionality:
EGA <- function (data, n = NULL, uni = TRUE,
                 model = c("glasso", "TMFG"), model.args = list(),
                 algorithm = c("walktrap", "louvain"), algorithm.args = list(),
                 plot.EGA = TRUE, plot.type = c("GGally", "qgraph"), plot.args = list(),...) {

  # Get additional arguments
  add.args <- list(...)

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

  ## Check for input plot arguments
  if(missing(plot.args)){
    plot.args <-list(vsize = 6, alpha = 0.4, label.size = 5, edge.alpha = 0.7)}

  else{
    plot.args <- plot.args
    plots.arg1 <- list(vsize = 6, label.size = 5, alpha = 0.4, edge.alpha = 0.7)
    plot.args.use <- plot.args

    if(any(names(plots.arg1) %in% names(plot.args.use))){

      plot.replace.args <- plots.arg1[na.omit(match(names(plot.args.use), names(plots.arg1)))]

      plot.args <- c(plot.args.use,plots.arg1[names(plots.arg1) %in% names(plot.args.use)==FALSE])}
    }


  #### ARGUMENTS HANDLING ####

  if(missing(model)){
    model <- "glasso"
  }else{model <- match.arg(model)}

  if(missing(algorithm)){
    algorithm <- "walktrap"
  }else if(!is.function(algorithm)){
    algorithm <- tolower(match.arg(algorithm))
  }

  if(missing(plot.type)){
    plot.type <- "GGally"
  }else{plot.type <- match.arg(plot.type)}

  #### ARGUMENTS HANDLING ####

  # Check for correlation matrix or data
  if(nrow(data) == ncol(data)){ ## Correlation matrix

    # Check for number of cases
    if(missing(n)){
      stop("There is no input for argument 'n'. Number of cases must be input when the matrix is square.")
    }

    # Make cor.data == data
    cor.data <- as.data.frame(data)

    # Multidimensional result
    ## Ensures proper partial correlations
    multi.res <- EGA.estimate(data = data, n = n,
                              model = model, model.args = model.args,
                              algorithm = algorithm, algorithm.args = algorithm.args)

    # Unidimensional result
    if(uni){

      # Set one factor for simulated data
      nfact <- 1
      nvar <- ncol(cor.data)
      if(nvar > 12){
        nvar <- 12
      }

      # Generate data
      uni.data <- MASS::mvrnorm(n = n, mu = rep(0, nvar), Sigma = cor.data)

      # Simulate data from unidimensional factor model
      sim.data <- sim.func(data = uni.data, nvar = nvar, nfact = nfact, load = .70)

      # Estimate unidimensional EGA
      uni.res <- suppressMessages(EGA.estimate(data = sim.data, n = n,
                                               model = model, model.args = model.args,
                                               algorithm = algorithm, algorithm.args = algorithm.args))

      # Set up results
      if(uni.res$n.dim <= nfact + 1){ ## If unidimensional

        n.dim <- uni.res$n.dim
        cor.data <- cor.data
        estimated.network <- multi.res$network
        wc <- uni.res$wc[-c(1:(nvar*nfact))]
        if(model == "glasso"){
          gamma <- uni.res$gamma
          lambda <- uni.res$lambda
        }

      }else{ ## If not

        n.dim <- multi.res$n.dim
        cor.data <- multi.res$cor.data
        estimated.network <- multi.res$network
        wc <- multi.res$wc
        if(model == "glasso"){
          gamma <- multi.res$gamma
          lambda <- multi.res$lambda
        }

      }

    }else{ ## Multidimensional check only

      n.dim <- multi.res$n.dim
      cor.data <- multi.res$cor.data
      estimated.network <- multi.res$network
      wc <- multi.res$wc
      if(model == "glasso"){
        gamma <- multi.res$gamma
        lambda <- multi.res$lambda
      }

    }

  }else{ ## Data

    # Convert to data frame
    data <- as.data.frame(data)

    # Get number of cases
    n <- nrow(data)

    # Check for unidimensional structure
    if(uni){

      # Set one factor for simulated data
      nfact <- 1
      nvar <- ncol(data)
      if(nvar > 12){
        nvar <- 12
      }

      ## Simulate data from unidimensional factor model
      data.sim <- sim.func(data = data, nvar = nvar, nfact = nfact, load = .70)

      ## Compute correlation matrix
      cor.data <- qgraph::cor_auto(data.sim)

      # Unidimensional result
      uni.res <- EGA.estimate(data = cor.data, n = n,
                              model = model, model.args = model.args,
                              algorithm = algorithm, algorithm.args = algorithm.args)

      ## Remove simulated data for multidimensional result
      cor.data <- cor.data[-c(1:nvar),-c(1:nvar)]

      # Multidimensional result
      multi.res <- suppressMessages(EGA.estimate(cor.data, n = n,
                                                 model = model, model.args = model.args,
                                                 algorithm = algorithm, algorithm.args = algorithm.args))

      if(uni.res$n.dim <= nfact + 1){

        n.dim <- uni.res$n.dim
        cor.data <- cor.data
        estimated.network <- multi.res$network
        wc <- uni.res$wc[-c(1:(nvar*nfact))]
        if(model == "glasso"){
          gamma <- uni.res$gamma
          lambda <- uni.res$lambda
        }

      }else{

        n.dim <- multi.res$n.dim
        cor.data <- cor.data
        estimated.network <- multi.res$network
        wc <- multi.res$wc

        if(model == "glasso"){
          gamma <- multi.res$gamma
          lambda <- multi.res$lambda
        }

      }


    }else{ ## Multidimensional check only

      ## Compute correlation matrix
      cor.data <- qgraph::cor_auto(data)

      # Multidimensional result
      multi.res <- suppressMessages(EGA.estimate(cor.data, n = n,
                                                 model = model, model.args = model.args,
                                                 algorithm = algorithm, algorithm.args = algorithm.args))

      n.dim <- multi.res$n.dim
      cor.data <- cor.data
      estimated.network <- multi.res$network
      wc <- multi.res$wc

      if(model == "glasso"){
        gamma <- multi.res$gamma
        lambda <- multi.res$lambda
      }

    }

  }

  a <- list()
  # Returning only communities that have at least two items:
  if(length(unique(wc))>1){
    indices <- seq_along(wc)
    indices2 <- indices[wc %in% wc[duplicated(wc)]]
    wc[indices[-indices2]] <- NA
    a$n.dim <- length(unique(na.omit(wc)))
  }else{
    a$n.dim <- length(unique(wc))
  }

  a$correlation <- cor.data
  a$network <- estimated.network
  a$wc <- wc
  # check if data has column names
  if(is.null(colnames(data)))
  {
    dim.variables <- data.frame(items = paste("V", 1:ncol(data), sep = ""), dimension = a$wc)
  }else{dim.variables <- data.frame(items = colnames(data), dimension = a$wc)}
  dim.variables <- dim.variables[order(dim.variables[, 2]),]
  a$dim.variables <- dim.variables
  if (plot.EGA == TRUE) {
    if (plot.type == "qgraph"){
      if(a$n.dim < 2){
        plot.ega <- qgraph::qgraph(a$network, layout = "spring",
                                   vsize = plot.args$vsize, groups = as.factor(a$wc), label.prop = 1, legend = FALSE)
      }else{
        plot.ega <- qgraph::qgraph(a$network, layout = "spring",
                                   vsize = plot.args$vsize, groups = as.factor(a$wc), label.prop = 1, legend = TRUE)
      }
    }else if(plot.type == "GGally"){
      if(a$n.dim <= 2){
        # weighted  network
        network1 <- network::network(a$network,
                                     ignore.eval = FALSE,
                                     names.eval = "weights",
                                     directed = FALSE)
        network::set.vertex.attribute(network1, attrname= "Communities", value = a$wc)
        network::set.vertex.attribute(network1, attrname= "Names", value = network::network.vertex.names(network1))
        network::set.edge.attribute(network1, "color", ifelse(network::get.edge.value(network1, "weights") > 0, "darkgreen", "red"))
        network::set.edge.value(network1,attrname="AbsWeights",value=abs(a$network))
        network::set.edge.value(network1,attrname="ScaledWeights",
                                value=matrix(scales::rescale(as.vector(a$network),
                                                             to = c(.001, 1.75)),
                                             nrow = nrow(a$network),
                                             ncol = ncol(a$network)))

        # Layout "Spring"
        graph1 <- NetworkToolbox::convert2igraph(a$network)
        edge.list <- igraph::as_edgelist(graph1)
        layout.spring <- qgraph::qgraph.layout.fruchtermanreingold(edgelist = edge.list,
                                                                   weights =
                                                                     abs(igraph::E(graph1)$weight/max(abs(igraph::E(graph1)$weight)))^2,
                                                                   vcount = ncol(a$network))

        set.seed(1234)
        plot.ega <- GGally::ggnet2(network1, edge.size = "ScaledWeights", palette = "Set1",
                                   color = "Communities", edge.color = c("color"),
                                   alpha = plot.args$alpha, #0.7,
                                   size = plot.args$vsize, #12,
                                   edge.alpha = plot.args$edge.alpha, #0.4,
                                   mode =  layout.spring,
                                   label.size = plot.args$label.size, #5
                                   label = colnames(a$network)) +
          ggplot2::theme(legend.title = ggplot2::element_blank())

        plot(plot.ega)

      }else{
        # weighted  network
        network1 <- network::network(a$network,
                                     ignore.eval = FALSE,
                                     names.eval = "weights",
                                     directed = FALSE)

        network::set.vertex.attribute(network1, attrname= "Communities", value = a$wc)
        network::set.vertex.attribute(network1, attrname= "Names", value = network::network.vertex.names(network1))
        network::set.edge.attribute(network1, "color", ifelse(network::get.edge.value(network1, "weights") > 0, "darkgreen", "red"))
        network::set.edge.value(network1,attrname="AbsWeights",value = abs(a$network))
        network::set.edge.value(network1,attrname="ScaledWeights",
                                value=matrix(scales::rescale(as.vector(a$network),
                                                             to = c(.001, 1.75)),
                                             nrow = nrow(a$network),
                                             ncol = ncol(a$network)))

        # Layout "Spring"

        graph1 <- NetworkToolbox::convert2igraph(a$network)
        edge.list <- igraph::as_edgelist(graph1)
        layout.spring <- qgraph::qgraph.layout.fruchtermanreingold(edgelist = edge.list,
                                                                   weights =
                                                                     abs(igraph::E(graph1)$weight/max(abs(igraph::E(graph1)$weight)))^2,
                                                                   vcount = ncol(a$network))

        set.seed(1234)
        plot.ega <- GGally::ggnet2(network1, edge.size = "ScaledWeights", palette = "Set1",
                                  color = "Communities", edge.color = c("color"),
                                  alpha = plot.args$alpha, #0.7,
                                  size = plot.args$vsize, #12,
                                  edge.alpha = plot.args$edge.alpha, #0.4,
                                  mode =  layout.spring,
                                  label.size = plot.args$label.size, #5
                                  label = colnames(a$network)) +
          ggplot2::theme(legend.title = ggplot2::element_blank(), legend.position = "none")

        plot(plot.ega)
      }
    }
  }else{plot.ega <- qgraph::qgraph(a$network, DoNotPlot = TRUE)}

  # check for variable labels in qgraph
  if(plot.type == "qgraph"){
    if(is.null(names(plot.ega$graphAttributes$Nodes$labels)))
    {names(plot.ega$graphAttributes$Nodes$labels) <- paste(1:ncol(data))}

    row.names(a$dim.variables) <- plot.ega$graphAttributes$Nodes$labels[match(a$dim.variables$items, names(plot.ega$graphAttributes$Nodes$labels))]
  }

  a$EGA.type <- ifelse(a$n.dim <= 2, "Unidimensional EGA", "Traditional EGA")
  a$Plot.EGA <- plot.ega

  # Get arguments
  args <- list()

  ## Get model and algorithm arguments
  args$model <- model
  args$algorithm <- algorithm

  ## Check if glasso was used
  if(model == "glasso")
  {
    args$gamma <- gamma
    args$lambda <- lambda
  }

  ## Check if walktrap was used
  if(!is.function(algorithm)){

    if(algorithm == "walktrap"){

      if("steps" %in% names(algorithm.args)){
        args$steps <- algorithm.args$steps
      }else{args$steps <- 4}

    }

  }

  a$Methods <- args

  class(a) <- "EGA"

  # Message that unidimensional structures were not checked
  if(!uni){
    message("\nEGA did not check for unidimensionality. Set argument 'uni' to TRUE to check for unidimensionality")
  }

  # Return estimates:
  return(a)
}
#----
