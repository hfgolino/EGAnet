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

  # (pre-computed 1 / pi)
  return(lambda * 0.3183099 * atan(abs(x) / gamma) + 0.5)

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
# Updated 04.02.2026
gumbel_penalty <- function(x, lambda, gamma = 0.01, ...)
{
  return(lambda * exp(-exp(-abs(x) / gamma)))
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
# Updated 13.01.2026
weibull_penalty <- function(x, lambda, gamma, shape, ...)
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
# Updated 25.07.2025
atan_derivative <- function(x, lambda, gamma = 0.01, ...)
{
  return(lambda * sign(x) * (gamma * (gamma + 2 / pi)) / (gamma^2 + abs(x)^2))
}

#' @noRd
# Updated 22.11.2025
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
# Updated 13.01.2026
cauchy_derivative <- function(x, lambda, gamma = 0.01, ...)
{

  # Return derivative (pre-computed 1 / pi)
  return(lambda * sign(x) * 0.3183099 * (gamma / (x^2 + gamma^2)))

}

#' @noRd
# Updated 10.01.2026
exp_derivative <- function(x, lambda, gamma = 0.01, ...)
{

  # Pre-compute components
  x <- abs(x)

  # Return penalty
  return(lambda * (1 / gamma) * exp(-(x / gamma)))

}

#' @noRd
# Updated 04.02.2026
gumbel_derivative <- function(x, lambda, gamma = 0.01, ...)
{

  # Pre-compute values
  gamma_x <- abs(x) / gamma

  # Return derivative
  return(lambda * sign(x) * (1 / gamma) * exp(-gamma_x - exp(-gamma_x)))

}

#' @noRd
# Updated 25.07.2025
l1_derivative <- function(x, lambda, ...)
{
  return(lambda * sign(x))
}

#' @noRd
# Updated 25.07.2025
l2_derivative <- function(x, lambda, ...)
{
  return(2 * lambda * x)
}

#' @noRd
# Updated 25.07.2025
mcp_derivative <- function(x, lambda, gamma = 3, ...)
{

  # Obtain absolute
  abs_x <- abs(x)

  # Return derivative
  return((abs_x <= (gamma * lambda)) * (lambda - abs_x / gamma) * sign(x))

}

#' @noRd
# Updated 08.08.2025
scad_derivative <- function(x, lambda, gamma = 3.7, ...)
{

  # Pre-compute components
  abs_x <- abs(x)
  sign_x <- x / abs_x
  sign_x <- swiftelse(is.na(sign_x), 0, sign_x)
  gamma_lambda <- gamma * lambda

  # Return derivative
  return(
    (abs_x <= lambda) * lambda * sign_x +  # region 1
      ((abs_x > lambda) & (abs_x <= gamma_lambda)) * # region 2
      (gamma_lambda - abs_x) * sign_x / (gamma - 1) +
      (abs_x > gamma_lambda) * 0  # region 3
  )

}

#' @noRd
# Updated 13.01.2026
weibull_derivative <- function(x, lambda, gamma, shape, ...)
{
  # Pre-compute components
  abs_x <- abs(x)
  x_gamma <- abs_x / gamma

  # Return derivative
  return(lambda * sign(x) * (shape / gamma) * x_gamma^(shape - 1) * exp(-x_gamma^shape))

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
weibull_proximal <- function(x, lambda, gamma, scale, ...)
{
  return(l1_proximal(x, weibull_derivative(x, lambda, gamma, scale)))
}