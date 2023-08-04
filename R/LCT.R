#' Loadings Comparison Test
#'
#' An algorithm to identify whether data were generated from a
#' factor or network model using factor and network loadings.
#' The algorithm uses heuristics based on theory and simulation. These
#' heuristics were then submitted to several deep learning neural networks
#' with 240,000 samples per model with varying parameters.
#'
#' @param data Matrix or data frame.
#' A data frame with the variables to be used in the test or a correlation matrix.
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
#' @param dynamic Boolean.
#' Is the dataset a time series where rows are time points and
#' columns are variables?
#' Defaults to \code{FASLE}.
#' 
#' @param dynamic.args List.
#' Arguments to be used in \code{\link[EGAnet]{dynEGA}}.
#' Defaults:
#' 
#' \itemize{
#' 
#' \item{\code{n.embed}}
#' {Number of embeddings: 4}
#' 
#' \item{\code{tau}}
#' {Lag: 1}
#' 
#' \item{\code{delta}}
#' {Delta: 1}
#' 
#' \item{\code{use.derivatives}}
#' {Derivatives: 1}
#' 
#' }
#' 
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
#' \item{proportion}{Proportions of models suggested across bootstraps}
#'
#' @examples
#' \donttest{# Compute LCT
#' ## Factor model
#' LCT(data = psychTools::bfi[,1:25])}
#' 
#' \dontrun{
#' # Dynamic LCT
#' LCT(sim.dynEGA[sim.dynEGA$ID == 1,1:24], dynamic = TRUE)}
#' 
#' 
#' @references
#' Christensen, A. P., & Golino, H. (2021).
#' Factor or network model? Predictions from neural networks.
#' \emph{Journal of Behavioral Data Science}, \emph{1}(1), 85-126.
#' 
#' @importFrom utils setTxtProgressBar txtProgressBar
#'
#' @noRd
#'
# Loadings Comparison Test----
# Updated 07.07.2022
LCT <- function (data, n, iter = 100,
                 dynamic = FALSE,
                 dynamic.args = list(
                   n.embed = 4, tau = 1, delta = 1,
                   use.derivatives = 1
                   )
                 )
{
  
  # Send not refactored warning (for now)
  not_refactored("LCT")
  
  # Convert data to matrix
  data <- as.matrix(data)
  
  # Number of cases
  if(nrow(data) == ncol(data)){
    if(missing(n))
    {stop("Argument 'n' must be supplied for an m x m matrix")}
    
    cases <- n
  }else{cases <- nrow(data)}
  
  # Initialize network loading matrix
  nl <- matrix(0, nrow = iter, ncol = 5)
  fl <- nl
  
  # Calculate total computations
  total_computations <- iter
  
  # Count computations
  count_computations <- 0
  
  # Initialize runtime updates
  runtime_update <- seq(0, total_computations, floor(total_computations / 100))
  runtime_update <- c(runtime_update, total_computations)
  
  # Obtain start time
  if(count_computations == 0){
    start_time <- Sys.time()
  }
  
  repeat{
    
    # Good sample?
    good <- FALSE
    
    while(!good){
      
      # Turn off pblapply progress bar
      if(isTRUE(dynamic)){
        opb <- pbapply::pboptions(type = "none")
      }
      
      # Generate data
      if(nrow(data) != ncol(data)) {
        
        if(count_computations == 1) {
          dat <- data
        } else {
          dat <- MASS_mvrnorm(cases, mu = rep(0, ncol(data)), Sigma = cov(data, use = "pairwise.complete.obs"))
          
        }
        
        # Static or dynamic
        if(isTRUE(dynamic)){
          
          # Organize time series
          dat <- dyn.org(data, dat)
          
          # Organize for dynamic EGA
          dat <- cbind(dat, rep(1, nrow(dat)), rep(1, nrow(dat)))
          colnames(dat)[(ncol(dat)-1):ncol(dat)] <- c("ID", "Group")
          dat <- as.data.frame(dat)
          
        }else{
          
          # Compute correlation
          cor.mat <- suppressMessages(
            auto.correlate(dat)
          )
          
        }
        
      }else{
        
        if(isTRUE(dynamic)){
          
          stop("Dynamic LCT requires the raw data. A correlation matrix cannot be used as input")
          
        }else{
          
          if(count_computations == 1) {
            cor.mat <- data
          } else {
            
            dat <- MASS_mvrnorm(cases, mu = rep(0, ncol(data)), Sigma = data)
            
            cor.mat <-  suppressMessages(
              auto.correlate(dat)
            )
          }
          
        }
        
      }
      
      # Make sure there are column names
      if(!isTRUE(dynamic)){
        if(is.null(colnames(cor.mat)))
        {colnames(cor.mat) <- paste("V", 1:ncol(cor.mat), sep = "")}
      }
      
      # Estimate network
      if(isTRUE(dynamic)){
        
        net <- try(suppressWarnings(suppressMessages(
          dynEGA(dat, n.embed = dynamic.args$n.embed,
                 tau = dynamic.args$tau, delta = dynamic.args$delta,
                 use.derivatives = dynamic.args$use.derivatives,
                 id = ncol(dat) - 1, group = ncol(dat),
                 model = "glasso", algorithm = "walktrap",
                 corr = "pearson", ncores = 2)$dynEGA
        )), silent = TRUE)
        
        cor.mat <- net$correlation
        
      }else{
        net <- try(suppressWarnings(suppressMessages(EGA(cor.mat, n = cases, plot.EGA = FALSE))), silent = TRUE)
      }
      
      if(any(class(net) == "try-error"))
      {good <- FALSE
      }else{
        
        if(length(net$wc) == length(unique(net$wc)))
        {good <- FALSE
        }else{
          
          # Remove variables missing dimension
          rm.NA <- which(is.na(net$wc))
          if(length(rm.NA) != 0){
            net$wc <- net$wc[-rm.NA]
            net$network <- net$network[-rm.NA, -rm.NA]
          }
          
          # Try to estimate network loadings
          n.loads <- try(abs(as.matrix(net.loads(net$network, net$wc)$std)), silent = TRUE)
          
          if(any(class(n.loads) == "try-error"))
          {good <- FALSE
          }else{
            
            # Reorder network loadings
            n.loads <- as.matrix(n.loads[match(names(net$wc), row.names(n.loads)),])
            
            # Get network loading proportions
            n.low <- mean(n.loads >= 0.15, na.rm = TRUE)
            n.mod <- mean(n.loads >= 0.25, na.rm = TRUE)
            n.high <- mean(n.loads >= 0.35, na.rm = TRUE)
            
            if(ncol(n.loads) != 1)
            {
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
              
              
            }else{
              n.dom <- NA
              n.cross <- NA
            }
            
            nl[count_computations,] <- c(n.low, n.mod, n.high, n.dom, n.cross)
            
            # Get factor loading proportions
            if(length(rm.NA) != 0){
              cor.mat <- cor.mat[-rm.NA, -rm.NA]
            }
            
            f.loads <- suppressWarnings(abs(as.matrix(psych::fa(cor.mat, nfactors = ncol(n.loads), n.obs = cases)$loadings[,1:ncol(n.loads)])))
            f.loads <- as.matrix(f.loads[match(names(net$wc), row.names(f.loads)),])
            f.low <- mean(f.loads >= 0.40, na.rm = TRUE)
            f.mod <- mean(f.loads >= 0.55, na.rm = TRUE)
            f.high <- mean(f.loads >= 0.70, na.rm = TRUE)
            
            # Organize loadings
            org <- numeric(ncol(cor.mat))
            
            for(i in 1:ncol(cor.mat))
            {org[i] <- which.max(f.loads[i,])}
            
            if(ncol(f.loads) != 1)
            {
              # Initialize dominate loadings
              f.dom <- numeric(ncol(cor.mat))
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
            }else{
              f.dom <- NA
              f.cross <- NA
            }
            
            fl[count_computations,] <- c(f.low, f.mod, f.high, f.dom, f.cross)
            
            # Increase count
            count_computations <- count_computations + 1
            
            # Update progress
            if(count_computations < 5){
              
              # Update progress
              custom_progress(
                i = count_computations,
                max = total_computations,
                start_time = "calculating"
              )
              
            }else if(count_computations %in% runtime_update){
              
              # Update progress
              custom_progress(
                i = count_computations,
                max = total_computations,
                start_time = start_time
              )
              
            }
            
            # Good data!
            good <- TRUE
          }
        }
        
      }
      
    }
    
    # Break out of repeat
    if(count_computations == (iter+1))
    {break}
  }
  
  # Convert to data frames
  loads.mat <- as.matrix(cbind(nl, fl))
  dimnames(loads.mat) <- NULL
  loads.mat <- ifelse(is.na(loads.mat), 0, loads.mat)
  
  # Predictions
  predictions <- list()
  
  # Without bootstrap
  wo.boot <- paste(dnn.predict(loads.mat[1,]))
  
  wo.boot <- switch(wo.boot,
                 "1" = "Factor",
                 "2" = "Network"
  )
  
  predictions$empirical <- wo.boot
  
  # Bootstrap prediction
  boot <- paste(dnn.predict(colMeans(loads.mat, na.rm = TRUE)))
  
  boot <- switch(boot,
                 "1" = "Factor",
                 "2" = "Network"
  )
  
  predictions$bootstrap <- boot
  
  # Bootstrap proportions
  boot.prop <- apply(na.omit(loads.mat), 1, dnn.predict)
  
  boot.prop <- colMeans(proportion.table(as.matrix(boot.prop)))
  
  prop <- vector("numeric", length = 2)
  names(prop) <- c("Factor", "Network")
  
  prop[1:length(boot.prop)] <- boot.prop
  
  predictions$proportion <- round(prop, 3)
  
  # Reset pboptions
  if(isTRUE(dynamic)){
    on.exit(pbapply::pboptions(opb))
  }
  
  return(predictions)
  
}
