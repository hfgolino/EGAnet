#' Loadings Comparison Test
#'
#' An algorithm to identify whether data were generated from a
#' random, factor, or network model using factor and network loadings.
#' The algorithm uses heuristics based on theory and simulation. These
#' heuristics were then submitted to several deep learning neural networks
#' with 240,000 samples per model with varying parameters.
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
#' @return Returns a list containing:
#' 
#' \item{empirical}{Prediction of model based on empirical dataset only}
#' 
#' \item{bootstrap}{Prediction of model based on means of the loadings across
#' the bootstrap replicate samples}
#' 
#' \item{bootstrapProportions}{Proportions of models suggested across bootstraps}
#' 
#' \item{omnibus}{An omnibus prediction based on a consensus of empirical,
#' bootstrap, and bootstrap proportions prediction. A consensus corresponds to
#' any combination of two predictions returning the same prediction}
#'
#' @examples
#' \donttest{# Compute LCT
#' ## Network model
#' LCT(data = wmt2[,7:24])$omnibus
#' 
#' ## Factor model
#' LCT(data = NetworkToolbox::neoOpen)$omnibus}
#' 
#' @references
#' Christensen, A. P., & Golino, H. (2020).
#' On the equivalency of factor and network loadings.
#' \emph{PsyArXiv}.
#' doi:\href{https://doi.org/10.31234/osf.io/xakez}{10.31234/osf.io/xakez}
#' 
#' @importFrom utils setTxtProgressBar txtProgressBar
#'
#' @export
#'
# Loadings Comparison Test----
# Updated 05.07.2020
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
  nl <- matrix(0, nrow = iter, ncol = 5)
  colnames(nl) <- c("Small", "Moderate", "Large", "Dominant", "Cross")
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
      if(nrow(data) != ncol(data)) {
        
        if(count == 1) {
          dat <- data
        } else {
          dat <- mvtnorm::rmvnorm(cases, sigma = cov(data, use = "pairwise.complete.obs"))
        }
      
        cor.mat <- cor(dat, use = "pairwise.complete.obs")
        
      }else{
        
        if(count == 1) {
          cor.mat <- data
        } else {
          
          dat <- mvtnorm::rmvnorm(cases, sigma = data)
          
          cor.mat <- cor(dat, use = "pairwise.complete.obs")
        }
        
      }
      
      # Make sure there are column names
      colnames(cor.mat) <- paste("V", 1:ncol(cor.mat), sep = "")
      
      # Estimate network
      net <- try(suppressWarnings(suppressMessages(EGA.estimate(cor.mat, n = cases))), silent = TRUE)
      
      if(any(class(net) == "try-error"))
      {good <- FALSE
      }else{
        
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
            n.cross <- ifelse(is.na(n.cross), 0, n.cross)
            
            nl[count,] <- c(n.low, n.mod, n.high, n.dom, n.cross)
            
            # Get factor loading proportions
            f.loads <- suppressWarnings(abs(as.matrix(psych::fa(cor.mat, nfactors = ncol(n.loads), n.obs = cases)$loadings[,1:ncol(n.loads)])))
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
            f.cross <- ifelse(is.na(f.cross), 0, f.cross)
            
            fl[count,] <- c(f.low, f.mod, f.high, f.dom, f.cross)
            
            # Increase count
            count <- count + 1
            
            # Update progress
            setTxtProgressBar(pb, count)
            
            # Good data!
            good <- TRUE
          }
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
  loads.mat <- as.matrix(cbind(nl, fl))
  dimnames(loads.mat) <- NULL
  loads.mat <- ifelse(is.na(loads.mat), 0, loads.mat)
  
  # Predictions
  predictions <- list()
  
  # Without bootstrap
  wo.boot <- paste(dnn.predict(loads.mat[1,]))
  
  wo.boot <- switch(wo.boot,
                    "1" = "Random",
                    "2" = "Factor",
                    "3" = "Network"
             )
                      
  predictions$empirical <- wo.boot
  
  # Bootstrap prediction
  boot <- paste(dnn.predict(colMeans(loads.mat, na.rm = TRUE)))
  
  boot <- switch(boot,
                 "1" = "Random",
                 "2" = "Factor",
                 "3" = "Network"
  )
  
  predictions$bootstrap <- boot
  
  # Bootstrap proportions
  boot.prop <- apply(na.omit(loads.mat), 1, dnn.predict)
  
  boot.prop <- colMeans(prop.table(as.matrix(boot.prop)))
  
  prop <- vector("numeric", length = 3)
  names(prop) <- c("Random", "Factor", "Network")
  
  prop[1:length(boot.prop)] <- boot.prop
  
  predictions$bootstrapProportions <- round(prop, 3)
  
  # Omnibus prediction
  omni.prop <- c(wo.boot, boot, names(prop)[which.max(prop)])
  omni.table <- table(omni.prop)
  
  if(any(omni.table > 1))
  {omni.pred <- names(omni.table)[which.max(omni.table)]
  }else{omni.pred <- "No consensus prediction. Check proportion and bootstrap predictions."}
  
  predictions$omnibus <- omni.pred
  
  return(predictions)
  
}
#----
