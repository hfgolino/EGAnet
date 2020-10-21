#' Automated Methods Section for \code{\link[EGAnet]{EGA}} Objects
#' 
#' @description This function accepts an \code{\link[EGAnet]{EGA}} object
#' and generates a Methods section for your analysis. The output is
#' an HTML containing the descriptions of the methods and parameters
#' along with a Reference section for appropriate citation.
#' 
#' @param EGA.object An \code{\link[EGAnet]{EGA}} object
#'
#' @return Automated HTML Methods section in your default browser
#' 
#' @examples
#' EGA.methods.section(ega.wmt)
#'
#' @export
#'
# Methods Section----
# Updated 09.10.2020
EGA.methods.section <- function(EGA.object)
{
  
  # Check for EGA object
  if(!class(EGA.object) %in% c("EGA", "bootEGA", "EGA.fit", "dynEGA"))
  {stop("Object input into 'EGA.object' is not an EGA object.")}
  
  #! Still need to work on bootEGA, EGA.fit, and dynEGA !#
  
  # Get class
  METHOD <- class(EGA.object)
  
  # Input arguments
  INPUT <- EGA.object$Methods
  
  # Year
  year <- unlist(strsplit(as.character(Sys.Date()), split = "\\-"))[1]
  
  # Version
  version <- packageVersion("EGAnet")
  
  # Version Year
  version.year <- packageDescription("EGAnet")
  version.year <- unlist(strsplit(as.character(version.year$Date), split = "\\-"))[1]
  
  # References
  refs <- list()
  
  refs$golino1 <- paste("Golino, H., & Christensen, A. P. (", version.year, "). ",
                        "EGAnet: Exploratory Graph Analysis -- A framework for estimating the number of dimensions in multivariate data using network psychometrics. ",
                        "Retrieved from https://cran.r-project.org/package=EGAnet",
                        sep = "")
  
  # Markdown YAML
  YAML <- c("---",
            "title: EGA Methods Section",
            "output: html_document",
            "---"
  )
  
  # For EGA
  model <- INPUT$model
  algorithm <- INPUT$algorithm
  
  # Set up text
  ## Introduction
  intro.header <- "# Exploratory Graph Analysis"
  intro.text <- paste("Exploratory graph analysis (EGA) is a recently developed method to estimate ",
                      "the number of dimensions in multivariate data using undirected network models ",
                      "(Golino & Epskamp, 2017; Golino et al., 2020). EGA first applies a network ",
                      "estimation method followed by a community detection algorithm for weighted ",
                      "networks (Fortunato, 2010). EGA has been shown to be as accurate or more accurate ",
                      "than more traditional factor analytic methods such as parallel analysis ",
                      "(Christensen & Golino, 2020; Golino et al., 2020). EGA was applied using the *EGAnet* package ",
                      "(version ", version, "; Golino & Christensen, ", version.year, ") in R (R Core Team, ", year, ").",
                      sep = "")
  
  refs$christensenB2020 <- paste("Christensen, A. P., & Golino, H. (2020).",
                                 "Estimating factors with psychometric networks: A Monte Carlo simulation comparing community detection algorithms.",
                                 "<em>PsyArXiv</em>.",
                                 "https://doi.org/10.31234/osf.io/hz89e")
  
  refs$golinoA2017 <- paste("Golino, H., & Epskamp, S. (2017).",
                            "Exploratory Graph Analysis: A new approach for estimating the number of dimensions in psychological research.",
                            "<em>PloS ONE</em>, <em>12</em>, e0174035.",
                            "https://doi.org/10.1371/journal.pone.0174035")
  
  refs$golinoB2020 <- paste("Golino, H., Shi, D., Christensen, A. P., Garrido, L. E., Nieto, M. D., Sadana, R., ... Martinez-Molina, A. (2020).",
                            "Investigating the performance of Exploratory Graph Analysis and traditional techniques to identify the number of latent factors: A simulation and tutorial.",
                            "<em>Psychological Methods</em>, <em>25</em>, 292--320.",
                            "https://doi.org/10.1037/met0000255")
  
  refs$fortunato2009 <- paste("Fortunato, S. (2010).",
                              "Community detectionin graphs.",
                              "<em>Physics Reports</em>, <em>3--5</em>, 75--174.",
                              "https://doi.org/10.1037/met0000255")
  
  ## Description of network estimation method
  model.header <- "## Network Estimation Method"
  
  if(tolower(model) == "glasso")
  {
    
    lambda <- INPUT$lambda
    gamma <- INPUT$gamma
    
    model.text <- paste("This study applied the graphical least absolute shrinkage and selection operator ",
                        "(GLASSO; Friedman, Haste, & Tibshirani, 2008, 2014), which estimates a Gaussian ",
                        "Graphical Model (GGM; Lauritzen, 1996) where nodes (circles) represent variables ",
                        "and edges (lines) represent the conditional dependence (or partial correlations) ",
                        "between nodes given all other nodes in the network. The least absolute shrinkage ",
                        "and selection operator (LASSO; Tibshirani, 1996) of the GLASSO is a regularization ",
                        "technique that reduces parameter estimates with some estimates becoming exactly zero. ",
                        "\n\n",
                        "The LASSO uses a parameter called lambda ($\\lambda$), which controls the sparsity of the network. ",
                        "Lower values of $\\lambda$ remove fewer edges, increasing the possibility of including ",
                        "spurious correlations, while larger values of $\\lambda$ remove more edges, increasing ",
                        "the possibility of removing relevant edges. When $\\lambda$ = 0, then the estimates are ",
                        "equal to the ordinary least squares solution for the partial correlation matrix. ",
                        "In this study, the ratio of the minimum and maximum $\\lambda$ was set to ", lambda, ".",
                        "\n\n",
                        "The popular approach in the network psychometrics literature is to compute models ",
                        "across several values of $\\lambda$ (usually 100) and to select the model that minimizes ",
                        "the extended Bayesian information criterion (EBIC; Chen & Chen, 2008; Epskamp & Fried, 2018). ",
                        "The EBIC model selection uses a hyperparameter gamma ($\\gamma$) to control how much it prefers simpler models ",
                        "(i.e., models with fewer edges; Foygel & Drton, 2010). Larger $\\gamma$ values lead to simpler models, ",
                        "while smaller $\\gamma$ values lead to denser models. When $\\gamma$ = 0, the EBIC is equal to the Bayesian ",
                        "information criterion. In this study, $\\gamma$ was set to ", gamma, ". In network psychometrics literature, ",
                        "this approach has been termed *EBICglasso* and is applied using the *qgraph* package (Epskamp et al., 2012) ",
                        "in R.",
                        sep = ""
    )
    
    refs$friedman2008 <- paste("Friedman, J., Hastie, T., & Tibshirani, R. (2008).",
                               "Sparse inverse covariance estimation with the graphical lasso.",
                               "<em>Biostatistics</em>, <em>9</em>, 432--441.",
                               "https://doi.org/10.1093/biostatistics/kxm045")
    
    refs$friedman2014 <- paste("Friedman, J., Hastie, T., & Tibshirani, R. (2014).",
                               "<em>glasso: Graphical lasso - estimation of Gaussian graphical models.</em>",
                               "Retrived from https://CRAN.R-project.org/package=glasso")
    
    refs$lauritzen1996 <- paste("Lauritzen, S. L. (1996).",
                                "<em>Graphical models.</em>",
                                "Oxford, UK: Clarendon Press.")
    
    refs$tibshirani1996 <- paste("Tibshirani, R. (1996).",
                                 "Regression shrinkage and selection via the lasso.",
                                 "<em>Journal of the Royal Statistical Society. Series B (Methodological)</em>, 267--288.",
                                 "https://doi.org/10.1111/j.2517-6161.1996.tb02080.x")
    
    refs$chen2008 <- paste("Chen, J., & Chen, Z. (2008).",
                           "Extended bayesian information criteria for model selection with large model spaces.",
                           "<em>Biometrika</em>, <em>95</em>, 759--771.",
                           "https://doi.org/10.1093/biomet/asn034")
    
    refs$epskampA2018 <- paste("Epskamp, S., & Fried, E. I. (2018).",
                               "A tutorial on regularized partial correlation networks.",
                               "<em>Psychological Methods</em>, <em>23</em>, 617--634.",
                               "https://doi.org/10.1037/met0000167")
    
    refs$foygel2010 <- paste("Foygel, R., & Drton, M. (2010).",
                             "Extended Bayesian information criteria for Gaussian graphical models.",
                             "In J. D. Lafferty, C. K. I. Williams, J. Shawe-Taylor, R. S., Zemel, & A. Culotta (Eds.),",
                             "<em>Advances in neural information processing systems</em> (pp. 604--612).",
                             "Retrieved from http://papers.nips.cc/paper/4087-extended-bayesianinformation-criteria-for-gaussian-graphical-models")
    
    refs$epskampB2012 <- paste("Epskamp, S., Cramer, A. O. J., Waldorp, L. J., Schmittmann, V. D., & Borsboom, D. (2012).",
                               "qgraph: Network visualizations of relationships in psychometric data.",
                               "<em>Journal of Statistical Software</em>, <em>48</em>, 1--18.",
                               "https://doi.org/10.18637/jss.v048.i04")
    
  }else if(model == "TMFG")
  {
    model.text <- paste("This study applied the Triangulated Maximally Filtered Graph (TMFG; Christensen et al., 2018; Massara, Di Matteo, & Aste, 2016), ",
                        "which applies a structural constraint on the zero-order correlation matrix. This constraint ",
                        "restrains the network to retain a certain number of edges (3*n*--6, where *n* is the number ",
                        "of nodes). This network is comprised of three- and four-node cliques (i.e., sets of connected ",
                        "nodes; a triangle and tetrahedron, respectively).",
                        "\n\n",
                        "Network estimation starts with a tetrahedron that is comprised of the four nodes ",
                        "that have the high sum of correlations that are greater than the average correlation in ",
                        "the correlation matrix. Next, the algorithm identifies a node that is not connected in the ",
                        "network and maximizes its sum of correlations to three nodes already in the network. This ",
                        "node is then connected to those three nodes. This process continues iteratively until ",
                        "every node is connected in the network.",
                        "\n\n",
                        "The resulting network is a *planar* network or a network that *could* be drawn ",
                        "such that no edges are crossing (Tumminello, Aste, Di Matteo, & Mantegna, 2005). One ",
                        'property of these networks is that they form a "nested hierarchy" such that its constituent ',
                        "elements (3-node cliques) form sub-networks in the overall network (Song, Di Matteo, & Aste, 2012). ",
                        "The TMFG method was applied using the *NetworkToolbox* package (Christensen, 2018) in R.",
                        sep = ""
    )
    
    refs$christensenC2018 <- paste("Christensen, A. P., Kenett, Y. N., Aste, T., Silvia, P. J., & Kwapil, T. R. (2018).",
                                   "Network structure of the Wisconsin Schizotypy Scales-Short Forms: Examining psychometric network filtering approaches.",
                                   "<em>Behavior Research Methods</em>, <em>50</em>, 2531--2550.",
                                   "https://doi.org/10.3758/s13428-018-1032-9")
    
    refs$massara2016 <- paste("Massara, G. P., Di Matteo, T., & Aste, T. (2016).",
                              "Network filtering for big data: Triangulated maximally filtered graph.",
                              "<em>Journal of Complex Networks</em>, <em>5</em>, 161--178.",
                              "https://doi.org/10.1093/comnet/cnw015")
    
    refs$tumminello2005 <- paste("Tumminello, M., Aste, T., Di Matteo, T., & Mantegna, R. N. (2005).",
                                 "A tool for filtering information in complex systems.",
                                 "<em>Proceedings of the National Academy of Sciences</em>, <em>102</em>, 10421--10426.",
                                 "https://doi.org/10.1073/pnas.0500298102")
    
    refs$song2012 <- paste("Song, W.-M., Di Matteo, T., & Aste, T. (2012).",
                           "Hierarchical information clustering by means of topologically embedded graphs.",
                           "<em>PLoS ONE</em>, <em>7</em>, e31929.",
                           "https://doi.org/10.1371/journal.pone.0031929")
    
    refs$christensenA2018 <- paste("Christensen, A. P. (2018).",
                                   "NetworkToolbox: Methods and measures for brain, cognitive, and psychometric network analysis in R.",
                                   "<em>The R Journal</em>, <em>10</em>, 422--439.",
                                   "https://doi.org/10.32614/RJ-2018-065")
  }
  
  ## Description of community detection algorithm
  algorithm.header <- "## Community Detection Algorithm"
  
  if(tolower(algorithm) == "walktrap")
  {
    steps <- INPUT$steps
    
    algorithm.text <- paste("The Walktrap algorithm (Pons & Latapy, 2006) is a commonly applied community detection algorithm in ",
                            "the psychometric network literature (Golino & Epskamp, 2017; Golino et al., 2020). The algorithm begins ",
                            "by computing a transition matrix where each element represents the probability of one node traversing to ",
                            "another (based on node strength or the sum of the connections to each node). Random walks are then initiated ",
                            "for a certain number of steps (e.g., ", steps, ") using the transition matrix for probable destinations. Using ",
                            "Ward's agglomerative clustering approach (Ward, 1963), each node starts as its own cluster and merges ",
                            "with adjacent clusters (based on squared distances between each cluster) in a way that minimizes the sum of ",
                            "squared distances between other clusters. Modularity (Newman, 2006) is then used to determine the optimal ",
                            "partition of clusters (i.e., communities). The Walktrap algorithm was implemented using the *igraph* ",
                            "package (Csardi & Nepusz, 2006) in R.",
                            sep = ""
    )
    
    refs$pons2006 <- paste("Pons, P., & Latapy, M. (2006).",
                           "Computing communities in large networks using random walks.",
                           "<em>Journal of Graph Algorithms and Applications</em>, <em>10</em>, 191--218.",
                           "https://doi.org/10.7155/jgaa.00185")
    
    refs$ward1963 <- paste("Ward, J. H. (1963).",
                           "Hierarchical clustering to optimise an objective function.",
                           "<em>Journal of the American Statistical Association</em>, <em>58</em>, 238--244.")
    
    refs$newman2006 <- paste("Newman, M. E. J. (2006).",
                             "Modularity and community structure in networks.",
                             "<em>Proceedings of the National Academy of Sciences</em>, <em>103</em>, 8577--8582.",
                             "https://doi.org/10.1073/pnas.0601602103")
    
    refs$csardi2006 <- paste("Csardi, G., & Nepusz, T. (2006).",
                             "The igraph software package for complex network research.",
                             "<em>InterJournal, Complex Systems</em>, <em>1695</em>, 1--9.",
                             "Retrieved from https://pdfs.semanticscholar.org/1d27/44b83519657f5f2610698a8ddd177ced4f5c.pdf")
    
    # EGA fit description
    #if(METHOD == "EGA.fit")
    #{
    #  fit.header <- "## Optimizing Fit"
    #  
    #  fit.text <- paste("",
    #                    sep = "")
    #}
    
  }else if(algorithm == "louvain")
  {
    algorithm.text <- paste("The Louvain algorithm (also referred to as Multi-level; Blondel, Guillaume, Lambiotte, & Lefebvre, 2008)",
                            "is one of the most commonly applied in network science (Gates, Henry, Steinley, & Fair, 2016). The algorithm",
                            "begins by randomly sorting nodes into communities with their neighbors and then uses",
                            "modularity (Newman, 2006) to iteratively optimize its communities partitions by exchanging nodes between communities",
                            "and evaluating the change in modularity until it no longer improves.",
                            "Then, the algorithm collapses the communities into latent nodes and identifies edge weights with other observed and latent",
                            "nodes, which provides a multi-level structure (Gates et al., 2016). In this study, the",
                            "algorithm was not used to identify hierarchical community structures in the network.",
                            "The Louvain algorithm was implemented using the *igraph* package (Csardi & Nepusz, 2006) in R.",
                            "It's also important to note that the algorithm implemented in *igraph* is deterministic;",
                            "however, other implementations are not (Gates et al., 2016)."
    )
    
    refs$blondel2008 <- paste("Blondel, V. D., Guillaume, J.-L., Lambiotte, R., & Lefebvre, E. (2008).",
                              "Fast unfolding of communities in large networks.",
                              "<em>Journal of Statistical Mechanics: Theory and Experiment</em>, <em>2008</em>, P10008.",
                              "https://doi.org/10.1088/1742-5468/2008/10/P10008")
    
    refs$gates2016 <- paste("Gates, K. M., Henry, T., Steinley, D., & Fair, D. A. (2016).",
                            "A Monte Carlo evaluation of weighted community detection algorithms.",
                            "<em>Frontiers in Neuroinformatics</em>, <em>10</em>, 45.",
                            "https://doi.org/10.3389/fninf.2016.00045")
    
    refs$newman2006 <- paste("Newman, M. E. J. (2006).",
                             "Modularity and community structure in networks.",
                             "<em>Proceedings of the National Academy of Sciences</em>, <em>103</em>, 8577--8582.",
                             "https://doi.org/10.1073/pnas.0601602103")
  }
  
  
  ## References
  references.header <- "## References"
  
  ## Order alphabetically
  references.text <- refs[order(names(refs))]
  
  references.text <- paste(references.text, collapse = "\n\n")
  
  markobj <- paste(YAML,
                   intro.header, intro.text,
                   model.header, model.text,
                   algorithm.header, algorithm.text,
                   references.header, references.text,
                   sep = "\n")
  
  tempHTML <- paste(tempdir(), "EGA_method.html", sep = "\\")
  tempHTML <- gsub("\\\\", "/", tempHTML)
  
  markdown::markdownToHTML(text = knitr::knit(text = markobj), output = tempHTML)
  
  browseURL(tempHTML)
  
}
