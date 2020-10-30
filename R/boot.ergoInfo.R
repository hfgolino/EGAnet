#' Bootstrap Test for the Ergodicity Information Index
#'
#' @description Computes a parametric Bootstrap Test for the Ergodicity Information Index, comparing the
#' empirical Ergodicity Information index to values obtained in data generated using \code{N} parametric bootstraps of the correlation matrix estimated using the
#' \code{\link[EGAnet]{dynEGA}} function, for the population structure. The p-values in the bootstrap test can be calculated as \code{(sum(EII>=boot.EII)+1)/(iter+1)} and as
#' \code{(sum(EII<=boot.EII)+1)/(iter+1)}, where EII is the empirical Ergodicity Information Index, boot.EII is the values of the Ergodicity Information Index obtained
#' in the bootstraped samples, and \code{iter} is the number of random samples generated in the simulation. The two-sided p-value is computed as two times the lowest p-value. In the bootstrap Test for the Ergodicity Information Index,
#' the null hypothesis is that the empirical value of EII is equal to the values of EII obtained in multiple individuals with the same structure as the population structure estimated
#' via \code{\link[EGAnet]{dynEGA}}.
#' Small values of p indicate that is very unlikely to obtain an EII as large as the one obtained in the empirical sample if the null hypothesis is true (i.e. all individuals have the same structure as the population structure), thus there is convincing evidence that the empirical Ergodicity Information Index is
#' different than it could be expected if all individuals had a similar latent structure.
#'
#' @param dynEGA.pop A dynEGA or a dynEGA.pop.ind object.
#'
#' @param iter Numeric integer.
#' Number of random samples to generate in the Monte-Carlo simulation.
#' At least \code{500} is recommended
#'
#' @param EII Numeric.
#' Empirical Ergodicity Information Index obtained via the \code{\link[EGAnet]{ergoInfo}} function.
#'
#' @param use Character.
#' A string indicating what network element will be used to compute the algorithm complexity in the \code{\link[EGAnet]{ergoInfo}} function,
#' the list of edges or the weights of the network.
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
#' @param corr Type of correlation matrix to compute. The default uses \code{\link[qgraph]{cor_auto}}.
#' Current options are:
#'
#' \itemize{
#'
#' \item{\strong{\code{cor_auto}}}
#' {Computes the correlation matrix using the \code{\link[qgraph]{cor_auto}} function from
#' \code{\link[qgraph]{qgraph}}}.
#'
#' \item{\strong{\code{pearson}}}
#' {Computes Pearson's correlation coefficient using the pairwise complete observations via
#' the \code{\link[stats]{cor}}} function.
#'
#' \item{\strong{\code{spearman}}}
#' {Computes Spearman's correlation coefficient using the pairwise complete observations via
#' the \code{\link[stats]{cor}}} function.
#' }
#'
#' @param model Character.
#' A string indicating the method to use. Defaults to \code{glasso}.
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
#' @param comparison String.
#' The population structure to be used as a reference in the Ergodicity Information Index should reflect the simulated data for all individuals (population),
#' or just a single individual (random.individual)?
#'
#'#' Current options are:
#'
#' \itemize{
#'
#' \item{\strong{\code{population.structure}}}
#' {Computes the ergodicity information index using the "population network" of the dynEGA object as
#' the reference (population) network.}
#'
#' \item{\strong{\code{boot.pop.structure}}}
#' {Computes the ergodicity information index having the population network estimated as the network structure of
#' the collection of bootstrap samples as the reference (population) network.}
#' }
#'
#' @param ... Additional arguments.
#' Used for deprecated arguments from previous versions of \code{\link{EGA}}
#'
#' @examples
#'
#' \dontrun{
#' \donttest{
#' dyn1 <- dynEGA.ind.pop(data = sim.dynEGA[,-c(22)], n.embed = 5, tau = 1,
#'                       delta = 1, id = 21, use.derivatives = 1,
#'                     model = "glasso", ncores = 2, corr = "pearson")
#'
#' eii1 <- ergoInfo(data = dyn1)$EII
#'
#' testing.ergoinfo <- boot.ergoInfo(dynEGA.pop = dyn1, iter = 10,EII = eii1,
#' model = "glasso", ncores = 2, corr = "pearson")
#' }}
#'
#' @return Returns a list containing:
#'
#' \item{boot.ergoInfo}{The values of the Ergodicity Information Index obtained in the Monte-Carlo Simulation}
#'
#' \item{p.value.twosided}{The p-value of the Monte-Carlo test for the Ergodicity Information Index.
#' The null hypothesis is that the empirical Ergodicity Information index is equal to the expected value of the EII if the all individuals
#' had similar latent structures.}
#'
#' \item{effect}{Indicates whether the empirical EII is greater or less then the bootstraped values EII, under the condition
#' that all individuals have a structure that is similar to the population structure.}
#'
#' \item{plot.dist}{Histogram of the bootstrapped ergodicity information index}
#'
#' @author Hudson Golino <hfg9s at virginia.edu>
#'
#' @export
# Bootstrap Test for the Ergodicity Information Index
# Updated 30.10.2020


boot.ergoInfo <- function(dynEGA.pop,
                            iter,
                            EII,
                            use,
                            model, model.args = list(),
                            algorithm = c("walktrap", "louvain"),
                            algorithm.args = list(),
                            corr, ncores,
                            comparison, ...
){


  #### MISSING ARGUMENTS HANDLING ####

  if(missing(dynEGA.pop))
  {  # Warning
    warning(
        "The 'dynEGA.pop' argument is missing. Please, provide the name of the object
        created using the 'dynEGA' or the 'dynEGA.pop.ind' functions."
    )
  }


  if(missing(use))
  {use <- "edge.list"
  }else{use}


  if(missing(model))
  {model <- "glasso"
  }else{model}

  if(missing(corr))
  {corr <- "cor_auto"
  }else{corr}

  if(missing(algorithm))
  {algorithm <- "walktrap"
  }else{algorithm}

  if(missing(ncores))
  {ncores <- ceiling(parallel::detectCores() / 2)
  }else{ncores}

  if(missing(comparison))
  {comparison <- "population.structure"
  }else{comparison}


  # Initialize Data list
  data.sim <- vector("list", length = iter)

  #let user know data generation has started
  message("\nGenerating the Data...\n", appendLF = FALSE)

  N <- nrow(dynEGA.pop$Derivatives$EstimatesDF)
  unique.ids <- unique(dplyr::last(dynEGA.pop$Derivatives$EstimatesDF))
  time.points <- floor(N/length(unique.ids))


 if(class(dynEGA.pop)=="dynEGA"){
   for(i in 1:iter){
     for(j in 1:time.points){
       data.sim[[i]][[j]] <- MASS::mvrnorm(n = time.points, mu = rep(0, ncol(dynEGA.pop$dynEGA$cor.dat)), Sigma = dynEGA.pop$dynEGA$cor.dat)
     }
   }
 } else if(class(dynEGA.pop)=="dynEGA.ind.pop"){
   for(i in 1:iter){
     for(j in 1:time.points){
     data.sim[[i]][[j]] <- MASS::mvrnorm(n = time.points, mu = rep(0, ncol(dynEGA.pop$dynEGA.pop$cor.data)), Sigma = as.matrix(Matrix::nearPD(dynEGA.pop$dynEGA.pop$cor.data)$mat))

   }
   }
 }

  data.sim.df <- vector("list", length = iter)
  for(i in 1:iter){
    data.sim.df[[i]] <- purrr::map_df(data.sim[[i]], ~as.data.frame(.))
    #data.sim.df[[i]]$ID <- rep(1:length(unique.ids), each = time.points)
  }

  #let user know data generation has ended
  message("done", appendLF = TRUE)

  #initialize correlation matrix list
  boots <- vector("list", length = iter)

  #Parallel processing
  cl <- parallel::makeCluster(ncores)

  #Export variables
  parallel::clusterExport(cl = cl,
                          varlist = c("data.sim.df", "N",
                                      "model", "model.args",
                                      "algorithm", "algorithm.args"),
                          envir=environment())

  #let user know data generation has started
  message("Estimating the Network Structures...\n", appendLF = FALSE)

  #Estimate networks
  boots <- pbapply::pblapply(
    X = data.sim.df, cl = cl,
    FUN = EGA,
    uni = FALSE,
    model = model, model.args = model.args,
    algorithm = algorithm, algorith.args = algorithm.args,
    plot.EGA = FALSE
  )

  parallel::stopCluster(cl)

   data.pop1 <- vector("list")
  if(class(dynEGA.pop)=="dynEGA"){
    data.pop1 <- dynEGA.pop$dynEGA}
  else{data.pop1 <- dynEGA.pop$dynEGA.pop}

   if(comparison == "population.structure"){
     data.pop1 <- data.pop1
   } else if(comparison == "boot.pop.structure"){
     data.pop1 <- data.pop1
     data <- purrr::map_df(data.sim.df, ~as.data.frame(.))
     ega.pop <- suppressMessages(suppressWarnings(EGAnet::EGA(data, plot.EGA = FALSE)))
     data.pop1$network <- ega.pop$network
     data.pop1$n.dim <- ega.pop$n.dim
   }



  ergoInfo2 <- function(data.ind, use = use, data.pop){


    ## Complexity - Individual Networks:
    # Number of edges:

    igraph.Network <- igraph::as.igraph(qgraph::qgraph(data.ind$network, DoNotPlot=TRUE))

    # Get the adjacency matrix:
    adj.net <- as.matrix(igraph::get.adjacency(igraph.Network,type="both"))

    # Transforming the matrices into a grid (ncol and nrow = ID 1)

    grid <- expand.grid(x = 1:ncol(data.ind$network), y = 1:ncol(data.ind$network))


    # Computing the Prime-Weight Transformation

    encode<- apply(grid, 1,
                                    function(x){return(
                                      ifelse(data.ind$network[x[1], x[2]]==0,1,(adj.net[x[1], x[2]]*2)^(data.ind$network[x[1], x[2]]))
                                    )})

    mat.encode <- matrix(encode, ncol = ncol(data.ind$network), nrow = nrow(data.ind$network))
    mat.encode <- ifelse(mat.encode==1, 0, mat.encode)
    mat.encode.igraph <- igraph::as.igraph(qgraph::qgraph(mat.encode, layout = "spring", DoNotPlot = TRUE))
    edge.list.mat.encode <- igraph::get.edgelist(mat.encode.igraph)

    if(use == "edge.list"){
      bits <- vector("list")
      compression <- vector("list")
      kcomp <- vector("list")

      for(i in 1:1000){
        bits[[i]] <- toString(edge.list.mat.encode[sample(1:nrow(edge.list.mat.encode), size = nrow(edge.list.mat.encode), replace = TRUE)])
        compression[[i]] <- memCompress(bits[[i]], "gzip")
        kcomp[[i]] <- length(compression[[i]])
      }
    }else{
      edge.list.mat.encode2 <- edge.list.mat.encode
      edge.list.mat.encode2 <- cbind(edge.list.mat.encode2, NA)
      for(i in 1:nrow(edge.list.mat.encode2)){
        edge.list.mat.encode2[i,3] <- mat.encode[edge.list.mat.encode2[i,1],edge.list.mat.encode2[i,2]]
      }

      bits <- vector("list")
      compression <- vector("list")
      kcomp <- vector("list")

      for(i in 1:1000){
        bits[[i]] <- toString(edge.list.mat.encode2[sample(1:nrow(edge.list.mat.encode2), size = nrow(edge.list.mat.encode2), replace = TRUE),3])
        compression[[i]] <- memCompress(bits[[i]], "gzip")
        kcomp[[i]] <- length(compression[[i]])
      }
    }


    ## Complexity - Population Network (reference):
    # Number of edges:

    igraph.Network.pop <- igraph::as.igraph(qgraph::qgraph(data.pop$network, DoNotPlot=TRUE))

    # Get the adjacency matrix:
    adj.net.pop <- as.matrix(igraph::get.adjacency(igraph.Network.pop,type="both"))

    # Transforming the matrices into a grid (ncol and nrow = ID 1)

    grid.pop <- expand.grid(x = 1:ncol(data.pop$network), y = 1:ncol(data.pop$network))


    # Computing the Prime-Weight Transformation

    encode.pop<- apply(grid, 1,
                   function(x){return(
                     ifelse(data.pop$network[x[1], x[2]]==0,1,(adj.net.pop[x[1], x[2]]*2)^(data.pop$network[x[1], x[2]]))
                   )})


    mat.encode.pop <- matrix(encode.pop, ncol = ncol(data.pop$network), nrow = nrow(data.pop$network))
    mat.encode.pop <- ifelse(mat.encode.pop==1, 0, mat.encode.pop)
    mat.encode.igraph.pop <- igraph::as.igraph(qgraph::qgraph(mat.encode.pop, layout = "spring", DoNotPlot = TRUE))
    edge.list.mat.encode.pop <- igraph::get.edgelist(mat.encode.igraph.pop)

    if(use == "edge.list"){
      edge.list.mat.encode.pop2 <- edge.list.mat.encode.pop
      edge.list.mat.encode.pop2 <- cbind(edge.list.mat.encode.pop2, NA)
      for(i in 1:nrow(edge.list.mat.encode.pop2)){
        edge.list.mat.encode.pop2[i,3] <- mat.encode[edge.list.mat.encode.pop2[i,1],edge.list.mat.encode.pop2[i,2]]
      }

      bits.pop <- vector("list")
      compression.pop <- vector("list")
      kcomp.pop <- vector("list")

      for(i in 1:1000){
        bits.pop[[i]] <- toString(edge.list.mat.encode.pop2[sample(1:nrow(edge.list.mat.encode.pop2), size = nrow(edge.list.mat.encode.pop2), replace = TRUE),3])
        compression.pop[[i]] <- memCompress(bits.pop[[i]], "gzip")
        kcomp.pop[[i]] <- length(compression.pop[[i]])
      }
    }else{

      bits.pop <- vector("list")
      compression.pop <- vector("list")
      kcomp.pop <- vector("list")


      for(i in 1:1000){
        bits.pop[[i]] <- toString(edge.list.mat.encode.pop[sample(1:nrow(edge.list.mat.encode.pop), size = nrow(edge.list.mat.encode.pop), replace = TRUE)])
        compression.pop[[i]] <- memCompress(bits.pop[[i]], "gzip")
        kcomp.pop[[i]] <- length(compression.pop[[i]])
      }
    }
    # Kolmogorov Complexity:
    results <- list()
    results$PrimeWeight.pop <- encode.pop
    results$KComp <- mean(unlist(kcomp))
    results$KComp.pop <- mean(unlist(kcomp.pop))
    ergo.info.index<- sqrt(data.pop$n.dim)^((results$KComp.pop/results$KComp)/log(sum(!results$PrimeWeight.pop==0)))
    results$EII <- ergo.info.index
    return(results)
  }



  #Parallel processing
  cl <- parallel::makeCluster(ncores)

  #Export variables
  parallel::clusterExport(cl = cl,
                          varlist = c("boots", "use",
                                      "data.pop1"),
                          envir=environment())

  #let user know data generation has started
  message("Estimating the Ergodicity Information Index\n", appendLF = FALSE)


  complexity.estimates <- pbapply::pblapply(X = boots, cl = cl,
                                            FUN = ergoInfo2, use = use, data.pop = data.pop1)
  parallel::stopCluster(cl)


  #let user know results are being computed
  message("Computing results...\n")


  complexity.estimates2 <- vector("list")
  for(i in 1:length(complexity.estimates)){
    complexity.estimates2[[i]] <- complexity.estimates[[i]]$EII
  }

  complexity.estimates2 <- unlist(complexity.estimates2)

  ## Compute the P-value of the bootstrap test:
  p.greater <- (sum(EII>=complexity.estimates2)+1)/(iter+1)
  p.lower <- (sum(EII<=complexity.estimates2)+1)/(iter+1)
  p.values <- c(p.greater, p.lower)
  two.sided <- 2*min(p.values)

  # Plot:
  complexity.df <- data.frame(EII = complexity.estimates2, ID = 1:length(complexity.estimates2))
  plot.bootErgoInfo <- suppressWarnings(suppressMessages(ggpubr::gghistogram(complexity.df, x = "EII",
                                                                             add = "mean",
                                                                             fill = "#00AFBB",
                                                                             color = "black",
                                                                             rug = TRUE)))

  ## Return Results:
  results <- vector("list")
  results$boot.ergoInfo <- complexity.estimates2
  results$p.value.twosided <- two.sided
  results$effect <- ifelse(p.greater<p.lower, "Greater", "Less")
  results$plot.dist <- plot.bootErgoInfo
  return(results)
}
#----
