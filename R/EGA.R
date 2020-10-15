#' Applies the Exploratory Graph Analysis technique
#'
#' Estimates the number of dimensions of a given dataset/instrument
#' using graphical lasso (\code{\link{EBICglasso.qgraph}}) or the
#' Triangulated Maximally Filtered Graph (\code{\link[NetworkToolbox]{TMFG}})
#' method and the walktrap community detection algorithm (\code{\link[igraph]{cluster_walktrap}}).
#' The glasso regularization parameter is set via EBIC.
#'
#' This algorithm includes checking for whether the data is unidimensional.
#'
#' @param data A dataframe with the variables to be used in the analysis or a correlation matrix.
#' If the data used is a correlation matrix, the argument \code{n} will need to be specified.
#'
#' @param n Integer.
#' Sample size, if the data provided is a correlation matrix
#'
#' @param plot.EGA Logical.
#' If TRUE, returns a plot of the network and its estimated dimensions.
#' Defaults to TRUE
#'
#' @param plot Character.
#' Plot system to use.
#' Current options are \code{\link[qgraph]{qgraph}} and \code{\link[GGally]{GGally}}.
#' Defaults to \code{"GGally"}.
#'
#' @param model A string indicating the method to use.
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
#' @param algorithm A string indicating the algorithm to use.
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
#' @param steps Number of steps to be used in \code{\link[igraph]{cluster_walktrap}} algorithm.
#' Defaults to 4.
#'
#' @param nvar Number of variables to use in the simulation part of the unidimensionality check. Defaults to 4.
#'
#' @param nfact Number of factors to be simulated (part of the unidimensionality check algorithm). Defaults to 1.
#'
#' @param load Factor loadings (used in the unidimensionality check algorithm). Defaults to 0.70.
#'
#' @param ... Additional arguments to be passed to \code{\link{EBICglasso.qgraph}}
#' or \code{\link[NetworkToolbox]{TMFG}}
#'
#' @author Hudson F. Golino <hfg9s at virginia.edu>, Alexander P. Christensen <alexpaulchristensen at gmail.com>, Maria Dolores Nieto <acinodam at gmail.com> and Luis E. Garrido <garrido.luiseduardo at gmail.com>
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
#'
#' \donttest{
#' \dontrun{
#' #estimate EGA
#' ega.wmt <- EGA(data = wmt2[,7:24], model = "glasso", plot.EGA = TRUE)
#' ega.wmt$Plot.EGA
#'
#' #estimate EGAtmfg
#' ega.wmt <- EGA(data = wmt2[,7:24], model = "TMFG", plot.EGA = TRUE)
#' ega.wmt$Plot.EGA
#'
#' #summary statistics
#' summary(ega.wmt)
#'
#' #plot
#' plot(ega.wmt)
#'
#' #estimate EGA
#' ega.intel <- EGA(data = intelligenceBattery[,8:66], model = "glasso", plot.EGA = TRUE)
#' ega.intel$Plot.EGA
#'
#' #summary statistics
#' summary(ega.intel)
#'
#' #plot
#' plot(ega.intel)
#' }
#' }
#'
#' @seealso \code{\link{bootEGA}} to investigate the stability of EGA's estimation via bootstrap
#' and \code{\link{CFA}} to verify the fit of the structure suggested by EGA using confirmatory factor analysis.
#'
#' @references
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
#' Golino, H., Shi, D., Christensen, A. P., Garrido, L. E., Nieto, M. D., Sadana, R., & Thiyagarajan, J. A. (in press).
#' Investigating the performance of Exploratory Graph Analysis and traditional techniques to identify the number of latent factors: A simulation and tutorial.
#' \emph{Psychological Methods}, \emph{25}, 292-320.
#' doi: \href{https://doi.org/10.1037/met0000255}{10.1037/met0000255}
#'
#' @importFrom stats cor rnorm runif na.omit
#'
#' @export
#'
# Updated 15.10.2020
## EGA Function to detect unidimensionality:
EGA <- function (data, model = c("glasso", "TMFG"),
                 algorithm = c("walktrap", "louvain"),
                 plot.EGA = TRUE, plot = c("GGally", "qgraph"), n = NULL,
                 steps = 4, nvar = 4, nfact = 1, load = 0.70, ...) {

  ##################################
  #### DATA SIMULATION FUNCTION ####
  ##################################

  sim.func <- function(data, nvar, nfact, load)
  {
    # Check for unidimensional structure
    ## Set up data simulation
    n <- nrow(data)
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
    Z = mvtnorm::rmvnorm(n, sigma = diag(J))                                  ## Obtain sample matrix of continuous variables
    X = Z%*%U
    colnames(X) <- paste0("X", 1:ncol(X))

    data.sim <- cbind(X, data)

    return(data.sim)
  }

  ##################################
  #### DATA SIMULATION FUNCTION ####
  ##################################

  #### MISSING ARGUMENTS HANDLING ####

  if(missing(model))
  {model <- "glasso"
  }else{model <- match.arg(model)}

  if(missing(algorithm))
  {algorithm <- "walktrap"
  }else{algorithm <- match.arg(algorithm)}

  if(missing(plot))
  {plot <- "GGally"
  }else{plot <- match.arg(plot)}

  #### MISSING ARGUMENTS HANDLING ####

  # Check for data or correlation matrix
  if(nrow(data) == ncol(data))
  {
    # Multidimensional correlation result
    multi.cor.res <- EGA.estimate(data = data, model = model, algorithm = algorithm, steps = steps, n = n, ...)

    # Unidimensional correlation result
    uni.data <- MASS::mvrnorm(n = n, mu = rep(0, ncol(data)), Sigma = multi.cor.res$cor.data)
    sim.data <- sim.func(data = uni.data, nvar = nvar, nfact = nfact, load = load)
    uni.cor.res <- suppressMessages(EGA.estimate(data = sim.data, model = model, algorithm = algorithm, steps = steps, n = n, ...))

    # Set up results
    if(uni.cor.res$n.dim <= nfact + 1)
    {
      n.dim <- uni.cor.res$n.dim
      cor.data <- multi.cor.res$cor.data
      estimated.network <- multi.cor.res$network
      wc <- uni.cor.res$wc[-c(1:(nvar*nfact))]
      if(model == "glasso")
      {
        gamma <- uni.res$gamma
        lambda <- uni.res$lambda
      }

    }else{
      n.dim <- multi.cor.res$n.dim
      cor.data <- multi.cor.res$cor.data
      estimated.network <- multi.cor.res$network
      wc <- multi.cor.res$wc
      if(model == "glasso")
      {
        gamma <- multi.res$gamma
        lambda <- multi.res$lambda
      }
    }

  }else{

    # Convert to data frame
    data <- as.data.frame(data)

    #-------------------------------------------------------------------------
    ## EGA WITH SIMULATED DATA + ORIGINAL DATA (UNIDIMENSIONALITY CHECK)
    #-------------------------------------------------------------------------

    n <- nrow(data)

    cor.data <- qgraph::cor_auto(data)

    data.sim <- sim.func(data = data, nvar = nvar, nfact = nfact, load = load)

    uni.res <- EGA.estimate(data.sim, model = model, algorithm = algorithm, steps = steps, n = n, ...)

    if(uni.res$n.dim <= nfact + 1)
    {
      n.dim <- uni.res$n.dim
      cor.data <- cor.data
      estimated.network <- suppressMessages(EGA.estimate(cor.data, model = model, algorithm = algorithm, steps = steps, n = n, ...)$network)
      wc <- uni.res$wc[-c(1:(nvar*nfact))]
      if(model == "glasso")
      {
        gamma <- uni.res$gamma
        lambda <- uni.res$lambda
      }

    }else{

      #-------------------------------------------------------------------------
      ## TRADITIONAL EGA (IF NUMBER OF FACTORS > 2)
      #-------------------------------------------------------------------------

      multi.res <- suppressMessages(EGA.estimate(cor.data, model = model, algorithm = algorithm, steps = steps, n = n, ...))

      n.dim <- multi.res$n.dim
      cor.data <- cor.data
      estimated.network <- multi.res$network
      wc <- multi.res$wc
      if(model == "glasso")
      {
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
    if (plot == "qgraph"){
    if(a$n.dim < 2){
      plot.ega <- qgraph::qgraph(a$network, layout = "spring",
                                 vsize = 6, groups = as.factor(a$wc), label.prop = 1, legend = FALSE)
    }else{
      plot.ega <- qgraph::qgraph(a$network, layout = "spring",
                                 vsize = 6, groups = as.factor(a$wc), label.prop = 1, legend = TRUE)
    }
      }else if(plot == "GGally"){
      if(a$n.dim <= 2){
        # weighted  network
        network1 <- network::network(a$network,
                                     ignore.eval = FALSE,
                                     names.eval = "weights",
                                     directed = FALSE)
        network::set.vertex.attribute(network1, attrname= "Communities", value = a$wc)
        network::set.vertex.attribute(network1, attrname= "Names", value = network::network.vertex.names(network1))
        network::set.edge.attribute(network1, "color", ifelse(get.edge.value(network1, "weights") > 0, "darkgreen", "red"))
        network::set.edge.value(network1,attrname="AbsWeights",value=abs(a$network))
        network::set.edge.value(network1,attrname="ScaledWeights",
                                value=matrix(scales::rescale(as.vector(a$network),
                                                             to = c(.001, 1.75)),
                                             nrow = nrow(a$network),
                                             ncol = ncol(a$network)))

        # Layout "Spring"
        graph1 <- igraph::as.igraph(qgraph::qgraph(a$network, DoNotPlot = TRUE))
        edge.list <- igraph::as_edgelist(graph1)
        layout.spring <- qgraph::qgraph.layout.fruchtermanreingold(edgelist = edge.list,
                                                                   weights =
                                                                     abs(E(graph1)$weight/max(abs(E(graph1)$weight)))^2,
                                                                   vcount = ncol(a$network))

        set.seed(1234)
        plot.ega <- GGally::ggnet2(network1, edge.size = "ScaledWeights", palette = "Set1",
                       color = "Communities", edge.color = c("color"),
                       alpha = 0.5, size = 6, edge.alpha = 0.5,
                       mode =  layout.spring,
                       label.size = 2.4,
                       label = colnames(a$network))+theme(legend.title = element_blank())
      }else{
        # weighted  network
        network1 <- network::network(a$network,
                                     ignore.eval = FALSE,
                                     names.eval = "weights",
                                     directed = FALSE)

        network::set.vertex.attribute(network1, attrname= "Communities", value = a$wc)
        network::set.vertex.attribute(network1, attrname= "Names", value = network::network.vertex.names(network1))
        network::set.edge.attribute(network1, "color", ifelse(get.edge.value(network1, "weights") > 0, "darkgreen", "red"))
        network::set.edge.value(network1,attrname="AbsWeights",value=abs(a$network))
        network::set.edge.value(network1,attrname="ScaledWeights",
                                value=matrix(scales::rescale(as.vector(a$network),
                                                             to = c(.001, 1.75)),
                                             nrow = nrow(a$network),
                                             ncol = ncol(a$network)))

        # Layout "Spring"
        graph1 <- igraph::as.igraph(qgraph::qgraph(a$network, DoNotPlot = TRUE))
        edge.list <- igraph::as_edgelist(graph1)
        layout.spring <- qgraph::qgraph.layout.fruchtermanreingold(edgelist = edge.list,
                                                                   weights =
                                                                     abs(E(graph1)$weight/max(abs(E(graph1)$weight)))^2,
                                                                   vcount = ncol(a$network))


        set.seed(1234)
        plot.ega <-GGally::ggnet2(network1, edge.size = "ScaledWeights", palette = "Set1",
                       color = "Communities", edge.color = c("color"),
                       alpha = 0.5, size = 6, edge.alpha = 0.5,
                       mode =  layout.spring,
                       label.size = 2.4,
                       label = colnames(a$network))+theme(legend.title = element_blank(), legend.position = "none")
      }
    }
  }else{plot.ega <- qgraph::qgraph(a$network, DoNotPlot = TRUE)}

  # check for variable labels in qgraph
  if(plot == "qgraph"){
    if(is.null(names(plot.ega$graphAttributes$Nodes$labels)))
    {names(plot.ega$graphAttributes$Nodes$labels) <- paste(1:ncol(data))}

    row.names(a$dim.variables) <- plot.ega$graphAttributes$Nodes$labels[match(a$dim.variables$items, names(plot.ega$graphAttributes$Nodes$labels))]
  }

  a$EGA.type <- ifelse(a$n.dim <= 2, "Unidimensional EGA", "Traditional EGA")
  a$Plot.EGA <- plot.ega
  a$Plot.EGA
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
  if(algorithm == "walktrap")
  {args$steps <- steps}

  a$Methods <- args

  class(a) <- "EGA"

  # Return estimates:
  return(a)
}
#----
