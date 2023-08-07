#' @noRd
# PRNG tests based on {RcppZiggurat}
# Results should be non-significant
# (with likely [false] positives around 0.05)
# See https://cran.r-project.org/web/packages/RcppZiggurat/vignettes/RcppZiggurat.pdf
# Updated 31.07.2023
rng_test <- function(
    FUN,
    cases = 10000, samples = 100, 
    seed = NULL, type = c("uniform", "normal")
)
{
  
  # Get total draws
  total_draws <- cases * samples
  
  # Initialize starting and ending indices
  start <- seq.int(1, total_draws, samples)
  end <- seq.int(samples, total_draws, samples)
  
  # Get arguments needed for function
  FUN_ARGS <- obtain_arguments(
    FUN, list(n = total_draws, seed = seed)
  )
  
  # Branch based on type
  if(type == "uniform"){

    # Get all values
    all_values <- do.call(
      what = FUN,
      args = FUN_ARGS
    )
    
    # Get sum of sample draws
    sample_sums <- nvapply(
      seq_along(start), function(i){
        sum(all_values[start[i]:end[i]])
      }
    )
    
    # Get probabilities
    probabilities <- pnorm(
      sample_sums,
      mean = samples / 2,
      sd = sqrt(samples / 12)
    )
    
  }else if(type == "normal"){
    
    # Get all values
    all_values <- do.call(
      what = FUN,
      args = FUN_ARGS
    )
    
    # Get sum of sample draws
    sample_sums <- nvapply(
      seq_along(start), function(i){
        sum(all_values[start[i]:end[i]])
      }
    )
    
    # Get probabilities
    probabilities <- pnorm(
      sample_sums,
      mean = 0,
      sd = sqrt(samples)
    )
    
  }
  
  # Return tests
  return(
    list(
      ks = ks.test(
        probabilities,
        "punif", 0, 1,
        exact = TRUE
      ),
      wilcox = wilcox.test(
        probabilities, mu = 0.5
      )
    )
  )
  
}





