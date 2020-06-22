#' Loadings Comparison Test
#'
#' An algorithm to identify whether data were generated from a
#' random, factor, or network model using factor and network loadings.
#' The algorithm uses heuristics based on theory and simulation.
#'
#' @param data Matrix or data frame.
#' A dataframe with the variables to be used in the test or a correlation matrix.
#' If the data used is a correlation matrix, the argument \code{n} will need to be specified
#'
#' @param n Integer.
#' Sample size (if the data provided is a correlation matrix)
#' 
#' @param iter Integer.
#' Number of replicate samples to be drawn from a multivariate
#' normal distribution (uses \code{mvtnorm::mvrnorm}).
#' Defaults to \code{100}
#'
#' @author Hudson F. Golino <hfg9s at virginia.edu> and Alexander P. Christensen <alexpaulchristensen at gmail.com>
#'
#' @return Returns a character string with the model suggested
#' from the three heuristics (see Details)
#'
#' @examples
#' \donttest{# Compute LCT
#' ## Network model
#' LCT(data = wmt2[,7:24])
#' 
#' ## Factor model
#' LCT(data = depression[,48:68])}
#' 
#' @references
#' Christensen, A. P., & Golino, H. (2020).
#' Statistical equivalency of factor and network loadings.
#' \emph{PsyArXiv}.
#' doi:\href{https://doi.org/10.31234/osf.io/xakez}{10.31234/osf.io/xakez}
#' 
#' @importFrom utils setTxtProgressBar txtProgressBar
#'
#' @export
#'
# Loadings Comparison Test----
# Updated 05.06.2020
LCT <- function (data, n, iter = 100)
{
  # Convert data to matrix
  data <- as.matrix(data)
  
  # Number of cases
  if(nrow(data) == ncol(data))
  {
    if(missing(n))
    {stop("Argument 'n' must be supplied for an m x m matrix")}
    
    cases <- n
  }else{cases <- nrow(data)}
  
  # Initialize network loading matrix
  nl <- matrix(0, nrow = iter, ncol = 6)
  colnames(nl) <- c("Below", "Small", "Moderate", "Large", "Dominant", "Cross")
  fl <- nl
  
  # Initialize count
  count <- 1
  
  # Initialize progress bar
  pb <- txtProgressBar(max = iter, style = 3)
  
  repeat{
    
    # Good sample?
    good <- FALSE
    
    while(!good)
    {
      # Generate data
      if(nrow(data) != ncol(data))
      {dat <- mvtnorm::rmvnorm(cases, sigma = cov(data, use = "pairwise.complete.obs"))
      }else{dat <- mvtnorm::rmvnorm(cases, sigma = data)}
      
      # Get correlation matrix
      cor.mat <- cor(dat)
      colnames(cor.mat) <- paste("V", 1:ncol(cor.mat), sep = "")
      
      # Estimate network
      net <- suppressWarnings(suppressMessages(EGA.estimate(cor.mat, n = cases)))
      
      if(length(unique(net$wc)) == 1 | length(net$wc) == length(unique(net$wc)))
      {good <- FALSE
      }else{
        # Try to estimate network loadings
        n.loads <- try(abs(as.matrix(net.loads(net$network, net$wc)$std)), silent = TRUE)
        
        if(any(class(n.loads) == "try-error"))
        {good <- FALSE
        }else if(ncol(n.loads) == 1)
        {good <- FALSE
        }else{
          
          # Reorder network loadings
          n.loads <- n.loads[names(net$wc),]
          
          # Get network loading proportions
          n.below <- mean(n.loads < 0.15, na.rm = TRUE)
          n.low <- mean(n.loads >= 0.15, na.rm = TRUE)
          n.mod <- mean(n.loads >= 0.25, na.rm = TRUE)
          n.high <- mean(n.loads >= 0.35, na.rm = TRUE)
          
          # Initialize dominate loadings
          n.dom <- numeric(ncol(data))
          n.loads2 <- n.loads
          
          for(i in 1:ncol(n.loads))
          {
            n.dom[which(net$wc == i)] <- n.loads[which(net$wc == i), i]
            n.loads2[which(net$wc == i), i] <- 0
          }
          
          # Get dominant and cross-loading proportions
          n.dom <- mean(n.dom >= 0.15)
          n.cross <- mean(ifelse(n.loads2 == 0, NA, n.loads2) >= 0.15, na.rm = TRUE)
          
          nl[count,] <- c(n.below, n.low, n.mod, n.high, n.dom, n.cross)
          
          # Get factor loading proportions
          f.loads <- suppressWarnings(abs(as.matrix(psych::fa(cor.mat, nfactors = ncol(n.loads), n.obs = cases)$loadings[,1:ncol(n.loads)])))
          f.below <- mean(f.loads < 0.40, na.rm = TRUE)
          f.low <- mean(f.loads >= 0.40, na.rm = TRUE)
          f.mod <- mean(f.loads >= 0.55, na.rm = TRUE)
          f.high <- mean(f.loads >= 0.70, na.rm = TRUE)
          
          # Organize loadings
          org <- numeric(ncol(data))
          
          for(i in 1:ncol(data))
          {org[i] <- which.max(f.loads[i,])}
          
          # Initialize dominate loadings
          f.dom <- numeric(ncol(data))
          f.loads2 <- f.loads
          
          for(i in 1:max(org))
          {
            f.dom[which(org == i)] <- f.loads[which(org == i), i]
            f.loads2[which(org == i), i] <- 0
          }
          
          # Get dominant and cross-loading proportions
          f.dom <- mean(f.dom >= 0.40)
          f.cross <- mean(ifelse(f.loads2 == 0, NA, f.loads2) >= 0.40, na.rm = TRUE)
          
          fl[count,] <- c(f.below, f.low, f.mod, f.high, f.dom, f.cross)
          
          # Increase count
          count <- count + 1
          
          # Update progress
          setTxtProgressBar(pb, count)
          
          # Good data!
          good <- TRUE
        }
      }
    }
    
    # Break out of repeat
    if(count == (iter+1))
    {break}
  }
  
  # Close progress bar
  close(pb)
  
  # Convert to data frames
  nl <- as.data.frame(nl)
  fl <- as.data.frame(fl)
  
  # Random vs. non-random model
  if(mean(fl$Large == fl$Moderate, na.rm = TRUE) >= .20)
  {return("Random")}
  
  # Ratios
  small.ratio <- mean(nl$Small / fl$Small, na.rm = TRUE)
  fl.ratio <- log(mean(fl$Dominant, na.rm = TRUE) / mean(fl$Cross, na.rm = TRUE))
  
  # Network vs. factor model
  if(small.ratio > 1.5)
  {return("Network")
  }else{
    # For sparse network models
    if(fl.ratio > 5 || is.infinite(fl.ratio))
    {return("Factor")
    }else{return("Network")}
  }
}
#----
