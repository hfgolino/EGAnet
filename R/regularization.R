#%%%%%%%%%%%%%%%%%%%%%%#
#### Regularization ####
#%%%%%%%%%%%%%%%%%%%%%%#

#%%%%%%%%%%%%%%%%
## Penalties ----
#%%%%%%%%%%%%%%%%

#' @noRd
# Updated 25.07.2025
atan_penalty <- function(x, lambda, gamma = 0.01, ...)
{
  return(lambda * (gamma + 2 / pi) * atan(abs(x) / gamma))
}

#' @noRd
# Updated 22.11.2025
bridge_penalty <- function(x, lambda, gamma = 1, ...)
{
  return(lambda * abs(x)^gamma)
}

#' @noRd
# Updated 13.01.2026
cauchy_penalty <- function(x, lambda, gamma = 0.01, ...)
{
  return(lambda * (1 / pi) * atan(abs(x) / gamma) + 0.5)
}

#' @noRd
# Updated 10.01.2026
exp_penalty <- function(x, lambda, gamma = 0.01, ...)
{

  # Pre-compute components
  x <- abs(x)

  # Return penalty
  return(lambda * (1 - exp(-(x / gamma))))

}

#' @noRd
# Updated 09.02.2026
gumbel_penalty <- function(x, lambda, gamma = 0.01, ...)
{

  # Pre-compute
  exp_1 <- exp(-1)

  return((lambda / (1 - exp_1)) * (exp(-exp(-abs(x) / gamma)) - exp_1))
  # theoretically, the `- exp(-1)` is necessary for the
  # penalty to converge at zero for the sparsity condition
  # in practice, this addition does not change the derivative,
  # which is used in the LLA
  # `1 - exp(-1)` is to scale lambda
}

#' @noRd
# Updated 25.07.2025
l1_penalty <- function(x, lambda, ...)
{
  return(lambda * abs(x))
}

#' @noRd
# Updated 25.07.2025
l2_penalty <- function(x, lambda, ...)
{
  return(lambda * x^2)
}

#' @noRd
# Updated 25.07.2025
mcp_penalty <- function(x, lambda, gamma = 3, ...)
{

  # Pre-compute components
  x <- abs(x)

  # Return penalty
  return(
    swiftelse(
      x <= (gamma * lambda),
      lambda * x - x^2 / (2 * gamma),
      gamma * lambda^2 / 2
    )
  )

}

#' @noRd
# Updated 25.07.2025
scad_penalty <- function(x, lambda, gamma = 3.7, ...)
{

  # Pre-compute components
  x <- abs(x)
  gamma_lambda <- gamma * lambda
  lambda_sq <- lambda^2

  # Return penalty
  return(
    (x <= lambda) * lambda * x +  # region 1
    ((x > lambda) & (x <= gamma_lambda)) * # region 2
    (-(x^2 - 2 * gamma * lambda * x + lambda_sq) / (2 * (gamma - 1))) +
    (x > gamma_lambda) * ((gamma + 1) * lambda_sq / 2)  # Region 3
  )

}

#' @noRd
# Updated 05.02.2026
weibull_penalty <- function(x, lambda, gamma = 0.01, shape, ...)
{

  # Pre-compute components
  x <- abs(x)

  # Return penalty
  return(lambda * (1 - exp(-(x / gamma)^shape)))

}

#%%%%%%%%%%%%%%%%%%
## Derivatives ----
#%%%%%%%%%%%%%%%%%%

#' @noRd
# Updated 27.02.2026
atan_derivative <- function(x, lambda, gamma = 0.01, ...)
{
  return(lambda * (gamma * (gamma + 2 / pi)) / (gamma^2 + x^2))
}

#' @noRd
# Updated 27.02.2026
bridge_derivative <- function(x, lambda, gamma = 1, eps = 1e-08, ...)
{

  # Check for gamma equal to zero
  if(gamma == 0){
    return(rep(lambda, length(x)))
  }

  # Otherwise, obtain absolute x
  abs_x <- pmax(abs(x), eps)

  # Return derivative
  return(lambda * gamma * x * abs_x^(gamma - 2))

}

#' @noRd
# Updated 27.02.2026
cauchy_derivative <- function(x, lambda, gamma = 0.01, ...)
{

  # Return derivative
  return(lambda * (1 / pi) * (gamma / (x^2 + gamma^2)))

}

#' @noRd
# Updated 27.02.2026
exp_derivative <- function(x, lambda, gamma = 0.01, ...)
{

  # Return penalty
  return(lambda * (1 / gamma) * exp(-(abs(x) / gamma)))

}

#' @noRd
# Updated 27.02.2026
gumbel_derivative <- function(x, lambda, gamma = 0.01, ...)
{

  # Pre-compute values
  gamma_x <- abs(x) / gamma

  # Return derivative
  return((lambda / (1 - exp(-1))) * (1 / gamma) * exp(-gamma_x - exp(-gamma_x)))

}

#' @noRd
# Updated 27.02.2026
l1_derivative <- function(x, lambda, ...)
{
  return(lambda)
}

#' @noRd
# Updated 27.02.2026
l2_derivative <- function(x, lambda, ...)
{
  return(2 * lambda * abs(x))
}

#' @noRd
# Updated 27.02.2026
mcp_derivative <- function(x, lambda, gamma = 3, ...)
{

  # Obtain absolute
  abs_x <- abs(x)

  # Return derivative
  return((abs_x <= (gamma * lambda)) * (lambda - abs_x / gamma))

}

#' @noRd
# Updated 27.02.2026
scad_derivative <- function(x, lambda, gamma = 3.7, ...)
{

  # Pre-compute components
  abs_x <- abs(x)
  gamma_lambda <- gamma * lambda

  # Return derivative
  return(
    (abs_x <= lambda) * lambda +  # region 1
    ((abs_x > lambda) & (abs_x <= gamma_lambda)) * # region 2
    (gamma_lambda - abs_x) / (gamma - 1) +
    (abs_x > gamma_lambda) * 0  # region 3
  )

}

#' @noRd
# Updated 27.02.2026
weibull_derivative <- function(x, lambda, gamma = 0.01, shape, ...)
{
  # Pre-compute components
  abs_x <- abs(x)
  x_gamma <- abs_x / gamma

  # Compute derivative
  derivative <- (shape / gamma) * x_gamma^(shape - 1) * exp(-x_gamma^shape)

  # Return derivative
  return(lambda * swiftelse(is.infinite(derivative), 100, derivative))

}

#%%%%%%%%%%%%%%%%%%%%%%%%%
## Proximal Operators ----
#%%%%%%%%%%%%%%%%%%%%%%%%%

#' @noRd
# Updated 13.01.2026
atan_proximal <- function(x, lambda, gamma = 0.01, ...)
{
  return(l1_proximal(x, atan_derivative(x, lambda, gamma)))
}

#' @noRd
# Updated 22.11.2025
bridge_proximal <- function(x, lambda, gamma = 1, eps = 1e-08, ...)
{

  # Check for gamma equal to zero
  if(gamma == 0){
    values <- rep(lambda, length(x))
  }else{

    # Compute absolute values and obtain values
    abs_x <- pmax(abs(x), eps)
    values <- lambda * gamma * abs_x^(gamma - 1)

  }

  # Return with l1 proximal
  return(l1_proximal(x, values))

}

#' @noRd
# Updated 13.01.2026
cauchy_proximal <- function(x, lambda, gamma = 0.01, ...)
{
  return(l1_proximal(x, cauchy_derivative(x, lambda, gamma)))
}

#' @noRd
# Updated 05.02.2026
exp_proximal <- function(x, lambda, gamma = 0.01, ...)
{
  return(l1_proximal(x, exp_derivative(x, lambda, gamma)))
}

#' @noRd
# Updated 05.02.2026
gumbel_proximal <- function(x, lambda, gamma = 0.01, ...)
{
  return(l1_proximal(x, gumbel_derivative(x, lambda, gamma)))
}

#' @noRd
# Updated 25.07.2025
l1_proximal <- function(x, lambda, ...)
{
  return(sign(x) * pmax(abs(x) - lambda, 0))
}

#' @noRd
# Updated 25.07.2025
l2_proximal <- function(x, lambda, ...)
{
  return(x / (1 + 2 * lambda))
}

#' @noRd
# Updated 25.07.2025
mcp_proximal <- function(x, lambda, gamma = 3, ...)
{

  return(
    swiftelse(
      abs(x) <= (gamma * lambda),
      l1_proximal(x, lambda) / (1 - 1 / gamma),
      x
    )
  )

}

#' @noRd
# Updated 25.07.2025
scad_proximal <- function(x, lambda, gamma = 3.7, ...)
{

  # Pre-compute components
  abs_x <- abs(x)
  gamma_lambda <- gamma * lambda
  gamma_one <- gamma - 1

  # Return values
  return(
    swiftelse(
      abs_x <= 2 * lambda,
      l1_proximal(x, lambda),
      swiftelse(
        abs_x <= gamma_lambda,
        (gamma_one / (gamma - 2)) * l1_proximal(x, gamma_lambda / gamma_one),
        x
      )
    )
  )

}

#' @noRd
# Updated 22.11.2025
weibull_proximal <- function(x, lambda, gamma = 0.01, scale, ...)
{
  return(l1_proximal(x, weibull_derivative(x, lambda, gamma, scale)))
}