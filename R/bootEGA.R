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
#' @param n Numeric integer.
#' Number of replica samples to generate from the bootstrap analysis.
#' At least \code{500} is recommended
#'
#' @param model Character.
#' A string indicating the method to use.
#' Defaults to \code{"glasso"}.
#'
#' Current options are:
#'
#' \itemize{
#'
#' \item{\strong{\code{"glasso"}}}
#' {Estimates the Gaussian graphical model using graphical LASSO with
#' extended Bayesian information criterion to select optimal regularization parameter.
#' See \code{\link[EGAnet]{EBICglasso.qgraph}}}
#'
#' \item{\strong{\code{"TMFG"}}}
#' {Estimates a Triangulated Maximally Filtered Graph.
#' See \code{\link[NetworkToolbox]{TMFG}}}
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
#' @param type Character.
#' A string indicating the type of bootstrap to use.
#'
#' Current options are:
#'
#' \itemize{
#'
#' \item{\strong{\code{"parametric"}}}
#' {Generates \code{n} new datasets (multivariate normal random distributions) based on the
#' original dataset, via the \code{\link[mvtnorm]{Mvnorm}} function of the mvtnorm package}
#'
#' \item{\strong{\code{"resampling"}}}
#' {Generates n random subsamples of the original data}
#'
#' }
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
#' #' @param plot.type Character.
#' Plot system to use.
#' Current options are \code{\link[qgraph]{qgraph}} and \code{\link[GGally]{GGally}}.
#' Defaults to \code{"GGally"}.
#'
#' @param ncores Numeric.
#' Number of cores to use in computing results.
#' Defaults to \code{parallel::detectCores() / 2} or half of your
#' computer's processing power.
#' Set to \code{1} to not use parallel computing.
#' Recommended to use maximum number of cores minus one
#'
#' If you're unsure how many cores your computer has,
#' then use the following code: \code{parallel::detectCores()}
#'
#' @param ... Additional arguments to be passed to \code{\link{EBICglasso.qgraph}}
#' or \code{\link[NetworkToolbox]{TMFG}}
#'
#' @return Returns a list containing:
#'
#' \item{n}{Number of replica samples in bootstrap}
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
#'
#'     }
#' }
#'
#' @author Hudson F. Golino <hfg9s at virginia.edu> and Alexander P. Christensen <alexpaulchristensen@gmail.com>
#'
#' @examples
#'
#' # Load data
#' wmt <- wmt2[,7:24]
#'
#' \dontrun{
#'
#' # bootEGA glasso example
#' boot.wmt <- bootEGA(data = wmt, n = 100, typicalStructure = TRUE,
#' plot.typicalStructure = TRUE, model = "glasso", type = "parametric", ncores = 4)
#' }
#'
#' # Load data
#' intwl <- intelligenceBattery[,8:66]
#'
#' \dontrun{
#' # bootEGA TMFG example
#' boot.intwl <- bootEGA(data = intelligenceBattery[,8:66], n = 100, typicalStructure = TRUE,
#' plot.typicalStructure = TRUE, model = "TMFG", type = "parametric", ncores = 4)
#'}
#'
#' @references
#' Christensen, A. P., & Golino, H. F. (2019).
#' Estimating the stability of the number of factors via Bootstrap Exploratory Graph Analysis: A tutorial.
#' \emph{PsyArXiv}.
#' doi:\href{https://doi.org/10.31234/osf.io/9deay}{10.31234/osf.io/9deay}
#'
#' @seealso \code{\link[EGAnet]{EGA}} to estimate the number of dimensions of an instrument using EGA
#' and \code{\link[EGAnet]{CFA}} to verify the fit of the structure suggested by EGA using confirmatory factor analysis.
#'
#' @importFrom stats cov median sd qt quantile
#'
#' @export
#'
# Bootstrap EGA
# Updated 10.15.2020
bootEGA <- function(data, n,
                    model = c("glasso", "TMFG"), algorithm = c("walktrap", "louvain"),
                    type = c("parametric", "resampling"),
                    typicalStructure = TRUE, plot.typicalStructure = TRUE,
                    plot.type = c("GGally", "qgraph"), ncores, ...) {

  #### MISSING ARGUMENTS HANDLING ####

  if(missing(model))
  {model <- "glasso"
  }else{model <- match.arg(model)}

  if(missing(algorithm))
  {algorithm <- "walktrap"
  }else{algorithm <- match.arg(algorithm)}

  if(missing(type))
  {type <- "parametric"
  }else{type <- match.arg(type)}

  if(missing(ncores))
  {ncores <- ceiling(parallel::detectCores() / 2)
  }else{ncores}

  if(missing(plot.type))
  {plot.type <- "GGally"
  }else{plot.type <- match.arg(plot.type)}

  #### MISSING ARGUMENTS HANDLING ####

  #number of cases
  cases <- nrow(data)

  #set inverse covariance matrix for parametric approach
  if(type=="parametric")  # Use a parametric approach:
  {
    if(model=="glasso")
    {
      g <- -EBICglasso.qgraph(qgraph::cor_auto(data), n = cases, lambda.min.ratio = 0.1, returnAllResults = FALSE, ...)
      diag(g) <- 1
    }else if(model=="TMFG")
    {
      g <- -NetworkToolbox::LoGo(data, normal = TRUE, partial=TRUE, ...)
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
    if(type == "parametric")
    {datalist[[count]] <- mvtnorm::rmvnorm(cases, sigma = corpcor::pseudoinverse(g))
    }else if(type == "resampling")
    {datalist[[count]] <- data[sample(1:cases, replace=TRUE),]}

    #break out of repeat
    if(count == n)
    {break}
  }

  #let user know data generation has ended
  message("done", appendLF = TRUE)

  #initialize correlation matrix list
  corlist <- list()

  #let user know data generation has started
  message("\nComputing correlation matrices...\n", appendLF = FALSE)

  #Parallel processing
  cl <- parallel::makeCluster(ncores)

  #Export variables
  parallel::clusterExport(cl = cl,
                          varlist = c("datalist", "corlist", "cases", ...),
                          envir=environment())

  #Compute correlation matrices
  corlist <- pbapply::pblapply(X = datalist, cl = cl,
                               FUN = qgraph::cor_auto)

  #let user know data generation has started
  message("Estimating networks...\n", appendLF = FALSE)

  #Estimate networks
  if(model == "glasso")
  {
    boots <- pbapply::pblapply(X = corlist, cl = cl,
                               FUN = EBICglasso.qgraph,
                               n = cases,
                               lambda.min.ratio = 0.1,
                               returnAllResults = FALSE,
                               ...)
  }else if(model == "TMFG")
  {
    boots <- pbapply::pblapply(X = corlist, cl = cl,
                               FUN = NetworkToolbox::TMFG,
                               normal = TRUE,
                               ...)

    for(i in 1:n)
    {boots[[i]] <- boots[[i]]$A}
  }

  parallel::stopCluster(cl)

  #let user know results are being computed
  message("Computing results...", appendLF = FALSE)

  bootGraphs <- vector("list", n)
  for (i in 1:n) {
    bootGraphs[[i]] <- boots[[i]]
    colnames(bootGraphs[[i]]) <- colnames(data)
    rownames(bootGraphs[[i]]) <- colnames(data)
  }
  boot.igraph <- vector("list", n)
  for (l in 1:n) {
    boot.igraph[[l]] <- NetworkToolbox::convert2igraph(abs(bootGraphs[[l]]))
  }
  boot.wc <- vector("list", n)
  for (m in 1:n) {

    boot.wc[[m]] <- switch(algorithm,
                           walktrap = igraph::cluster_walktrap(boot.igraph[[m]]),
                           louvain = igraph::cluster_louvain(boot.igraph[[m]])
                           )
  }
  boot.ndim <- matrix(NA, nrow = n, ncol = 2)
  for (m in 1:n) {
    boot.ndim[m, 2] <- max(boot.wc[[m]]$membership)
  }

  colnames(boot.ndim) <- c("Boot.Number", "N.Dim")

  boot.ndim[, 1] <- seq_len(n)
  if (typicalStructure == TRUE) {
    if(model=="glasso")
    {typical.Structure <- apply(simplify2array(bootGraphs),1:2, median)
    }else if(model=="TMFG")
    {typical.Structure <- apply(simplify2array(bootGraphs),1:2, mean)}
    typical.igraph <- NetworkToolbox::convert2igraph(abs(typical.Structure))

    typical.wc <- switch(algorithm,
                         walktrap = igraph::cluster_walktrap(typical.igraph),
                         louvain = igraph::cluster_louvain(typical.igraph)
                         )

    typical.ndim <- max(typical.wc$membership)
    dim.variables <- data.frame(items = colnames(data), dimension = typical.wc$membership)
  }
  if (plot.typicalStructure == TRUE) {
    if(plot.type == "qgraph"){
      plot.typical.ega <- qgraph::qgraph(typical.Structure, layout = "spring",
                                         vsize = 6, groups = as.factor(typical.wc$membership))
    }else if(plot.type == "GGally"){
        network1 <- network::network(typical.Structure,
                                     ignore.eval = FALSE,
                                     names.eval = "weights",
                                     directed = FALSE)

      network::set.vertex.attribute(network1, attrname= "Communities", value = typical.wc$membership)
      network::set.vertex.attribute(network1, attrname= "Names", value = network::network.vertex.names(network1))
      network::set.edge.attribute(network1, "color", ifelse( network::get.edge.value(network1, "weights") > 0, "darkgreen", "red"))
      network::set.edge.value(network1,attrname="AbsWeights",value=abs(typical.Structure))
      network::set.edge.value(network1,attrname="ScaledWeights",
                              value=matrix(scales::rescale(as.vector(typical.Structure),
                                                           to = c(.001, 1.75)),
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
      plot.typical.ega <-GGally::ggnet2(network1, edge.size = "ScaledWeights", palette = "Set1",
                                color = "Communities", edge.color = c("color"),
                                alpha = 0.5, size = 6, edge.alpha = 0.5,
                                mode =  layout.spring,
                                label.size = 2.4,
                                label = colnames(typical.Structure))+ggplot2::theme(legend.title = ggplot2::element_blank())
      plot(plot.typical.ega)
      }


  }
  Median <- median(boot.ndim[, 2])
  se.boot <- sd(boot.ndim[, 2])
  ciMult <- qt(0.95/2 + 0.5, nrow(boot.ndim) - 1)
  ci <- se.boot * ciMult
  quant <- quantile(boot.ndim[,2], c(.025, .975), na.rm = TRUE)
  summary.table <- data.frame(n.Boots = n, median.dim = Median,
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
    lik[count,2] <- length(which(boot.ndim[,2]==i))/n
  }

  #let user know results have been computed
  message("done", appendLF = TRUE)

  result <- list()
  result$n <- n
  result$boot.ndim <- boot.ndim
  result$boot.wc <- boot.wc
  result$bootGraphs <- bootGraphs
  result$summary.table <- summary.table
  result$frequency <- lik
  result$EGA <- suppressMessages(suppressWarnings(EGA(data = data, model = model, plot.EGA = FALSE)))

  # Typical structure
  if (typicalStructure == TRUE) {
    typicalGraph <- list()
    typicalGraph$graph <- typical.Structure
    typicalGraph$typical.dim.variables <- dim.variables[order(dim.variables[,2]), ]
    typicalGraph$wc <- typical.wc$membership
    result$typicalGraph <- typicalGraph
    result$plot.typical.ega <- plot.typical.ega
  }

  class(result) <- "bootEGA"
  return(result)
}
#----
