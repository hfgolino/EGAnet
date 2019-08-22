#' Apply the Exploratory Graph Analysis technique
#'
#' Estimates the number of dimensions of a given dataset/instrument
#' using graphical lasso (\code{\link{EBICglasso.qgraph}}) or the
#' Triangulated Maximally Filtered Graph (\code{\link[NetworkToolbox]{TMFG}})
#' method and the walktrap community detection algorithm (\code{\link[igraph]{cluster_walktrap}}).
#' The glasso regularization parameter is set via EBIC.
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
#' @examples
#'
#' \donttest{
#' #estimate EGA
#' ega.wmt <- EGA(data = wmt2[,7:24], model = "glasso", plot.EGA = TRUE)
#'
#'
#' #estimate EGAtmfg
#' ega.wmt <- EGA(data = wmt2[,7:24], model = "TMFG", plot.EGA = TRUE)
#'
#' #summary statistics
#' summary(ega.wmt)
#'
#' #plot
#' plot(ega.wmt)
#'
#' #estimate EGA
#' ega.intel <- EGA(data = intelligenceBattery[,8:66], model = "glasso", plot.EGA = TRUE)
#'
#' #summary statistics
#' summary(ega.intel)
#'
#' #plot
#' plot(ega.intel)
#' }
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
#' Golino, H., Shi, D., Garrido, L. E., Christensen, A. P., Nieto, M. D., Sadana, R., & Thiyagarajan, J. A. (2018).
#' Investigating the performance of Exploratory Graph Analysis and traditional techniques to identify the number of latent factors: A simulation and tutorial.
#' \emph{PsyArXiv}.
#' doi: \href{https://psyarxiv.com/gzcre/}{10.31234/osf.io/gzcre}
#'
#' @importFrom iterators iter nextElem
#' @importFrom stats cor rnorm runif na.omit
#'
#' @export
#'
# EGA default function - 05/30/2019
## EGA Function to detect unidimensionality:

EGA <- function (data, model = c("glasso", "TMFG"), plot.EGA = TRUE, n = NULL,
                 steps = 4, nvar = 4, nfact = 1, load = 0.70, ...) {

    #-------------------------------------------------------------------------
    ## CHECKING DATA STRUCTURE
    #-------------------------------------------------------------------------
    if(nrow(data)!=ncol(data)){
        data <- as.data.frame(data)
      n <- nrow(data)
    }

    if (nrow(data) == ncol(data)){
        cor.data <- data
        if (missing(model)) {
            model = "glasso"
        }
        else {
            model = match.arg(model)
        }
        if (model == "glasso") {
            gamma.values <- c(0.50, 0.25, 0)
            gvals <- iterators::iter(gamma.values)
            repeat {
                estimated.network <- EBICglasso.qgraph(data = cor.data,
                                                       n = n,
                                                       lambda.min.ratio = 0.1,
                                                       returnAllResults = FALSE,
                                                       gamma = iterators::nextElem(gvals), ...)
                if (any(NetworkToolbox::strength(estimated.network)>0)) {
                    break
                }
            }
        }
        else if (model == "TMFG") {
            estimated.network <- NetworkToolbox::TMFG(cor.data, ...)$A
        }
        graph <- NetworkToolbox::convert2igraph(abs(estimated.network))
        wc <- igraph::walktrap.community(graph, steps = steps)
        names(wc$membership) <- colnames(data)
        n.dim <- max(wc$membership)
        wc <- wc$membership
    } else{


        #-------------------------------------------------------------------------
        ## DATA GENERATION
        #-------------------------------------------------------------------------

        nvar <- nvar
        nfac <- nfact
        n <- nrow(data)
        load <- load
        corf <- 0
        J       = nvar*nfac
        sdcross = 0

        ## GENERATE SAMPLE DATA MATRIX
        check.eig <- TRUE
        check.com <- TRUE

        while(check.eig == TRUE|check.com == TRUE){

            SATF = matrix(0, J, nfac)

            for(j in 1:nfac){
                SATF[(j*nvar-nvar+1):(j*nvar),j]<-runif(nvar, load-.10, load+.10)
                if(nfac>1){
                    CROSS.L = apply(as.matrix(SATF[(j*nvar-nvar+1+2):(j*nvar),-c(j)]), ## Generate cross-loadings
                                    2,
                                    function(x) rnorm((nvar-2), 0, sdcross))
                    SATF[(j*nvar-nvar+1+2):(j*nvar),-c(j)] = CROSS.L
                }
            }
            SATF # Population factor loading matrix with cross-loadings and marker items

            FCOR      = matrix(corf, nfac, nfac); diag(FCOR)<-1 ## Factor correlation matrix
            R         = SATF%*%FCOR%*%t(SATF)                          ## Rr
            check.com = any(diag(R) > .90)                                  ## Check communalities values
            diag(R)   = 1                                                                    ## Insert ones in the diagonal of Rr
            R                                                                                       ## Rp
            check.eig = any(eigen(R)$values <= 0)                      ## Check eigenvalues
        }

        U = chol(R)                                                                       ## Cholesky decomposition of Rp
        Z = mvtnorm::rmvnorm(n, sigma = diag(J))                                  ## Obtain sample matrix of continuous variables
        X = Z%*%U
        X <-as.data.frame(X)
        colnames(X) <- paste0("V", 1:ncol(X))

        data.sim <- data.frame(X, data)

        #-------------------------------------------------------------------------
        ## EGA WITH SIMULATED DATA + ORIGINAL DATA
        #-------------------------------------------------------------------------

        if (missing(model)) {
            model = "glasso"
        }
        else {
            model = match.arg(model)
        }

        cor.data.sim <- qgraph::cor_auto(data.sim)

        if (model == "glasso") {
            gamma.values <- c(0.50, 0.25, 0)
            gvals <- iterators::iter(gamma.values)
            repeat {
                estimated.network.sim <- EBICglasso.qgraph(data = cor.data.sim,
                                                           n = n,
                                                           lambda.min.ratio = 0.1,
                                                           returnAllResults = FALSE,
                                                           gamma = iterators::nextElem(gvals), ...);
                if (any(NetworkToolbox::strength(estimated.network.sim)>0)) {
                    break
                }
            }
        }
        else if (model == "TMFG") {
            estimated.network.sim <- NetworkToolbox::TMFG(cor.data.sim, ...)$A
        }
        graph.sim <- NetworkToolbox::convert2igraph(abs(estimated.network.sim))
        wc.sim <- igraph::walktrap.community(graph.sim, steps = steps)
        names(wc.sim$membership) <- colnames(data.sim)
        n.dim.sim <- max(wc.sim$membership)
        if(n.dim.sim <= 2){
            n.dim <- n.dim.sim
            cor.data <- cor.data.sim[-c(1:nvar),-c(1:nvar)]
            estimated.network <- estimated.network.sim[-c(1:nvar),-c(1:nvar)]
            wc <- wc.sim$membership[-c(1:nvar)]
        }
        #-------------------------------------------------------------------------
        ## TRADITIONAL EGA (IF NUMBER OF FACTORS > 2)
        #-------------------------------------------------------------------------

        else{

            cor.data <- cor.data.sim[-c(1:nvar),-c(1:nvar)]

            if (missing(model)) {
                model = "glasso"
            }
            else {
                model = match.arg(model)
            }
            if (model == "glasso") {
                gamma.values <- c(0.50, 0.25, 0)
                gvals <- iterators::iter(gamma.values)
                repeat {
                    estimated.network <- EBICglasso.qgraph(data = cor.data,
                                                           n = n,
                                                           lambda.min.ratio = 0.1,
                                                           returnAllResults = FALSE,
                                                           gamma = iterators::nextElem(gvals), ...)
                    if (any(NetworkToolbox::strength(estimated.network)>0)) {
                        break
                    }
                }
            }
            else if (model == "TMFG") {
                estimated.network <- NetworkToolbox::TMFG(cor.data, ...)$A
            }
            graph <- NetworkToolbox::convert2igraph(abs(estimated.network))
            wc <- igraph::walktrap.community(graph, steps = steps)
            names(wc$membership) <- colnames(data)
            n.dim <- max(wc$membership)
            wc <- wc$membership
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
    dim.variables <- data.frame(items = colnames(data), dimension = a$wc)
    dim.variables <- dim.variables[order(dim.variables[, 2]),                                ]
    a$dim.variables <- dim.variables
    if (plot.EGA == TRUE) {
        if(a$n.dim <= 2){
            plot.ega <- qgraph::qgraph(a$network, layout = "spring",
                                       vsize = 6, groups = as.factor(a$wc), label.prop = 1, legend = FALSE)
        }else{
            plot.ega <- qgraph::qgraph(a$network, layout = "spring",
                                       vsize = 6, groups = as.factor(a$wc), label.prop = 1, legend = TRUE)
        }
    }

  a$EGA.type <- ifelse(a$n.dim <= 2, "Unidimensional EGA", "Traditional EGA")
  class(a) <- "EGA"
  return(a)
}
#----
