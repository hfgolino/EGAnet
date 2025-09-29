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
# Updated 27.07.2025
dbridge_penalty <- function(x, lambda, gamma = 5, ...)
{
  return(lambda * abs(x)^(2^(-gamma)))
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
# Updated 12.01.2025
lgp_penalty <- function(x, lambda, gamma = 5, ...)
{
  return(lambda * abs(x)^(lambda / gamma))
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
# Updated 12.01.2025
pop_penalty <- function(x, lambda, gamma = 4, ...)
{

  # Return lambdas
  return(lambda * (1 - (1 / (abs(x) + 1))^gamma))

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
# Updated 28.01.2025
spot_penalty <- function(x, lambda, gamma = 3, ...)
{
  return(2 * lambda / (1 + exp(-abs(x) * 2^gamma)) - lambda)
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
# Updated 27.07.2025
dbridge_derivative <- function(x, lambda, gamma = 5, ...)
{

  # Bridge value
  bridge_value <- 2^(-gamma)

  # Return derivative
  return(swiftelse(gamma == 0, lambda, lambda * bridge_value * abs(x)^(bridge_value - 1)))

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
# Updated 12.01.2025
lgp_derivative <- function(x, lambda, gamma = 5, ...)
{
  return(lambda^2 * abs(x)^((lambda - gamma) / gamma) / gamma)
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
# Updated 12.01.2025
pop_derivative <- function(x, lambda, gamma = 4, ...)
{

  # Return lambdas
  return((lambda * gamma) / (abs(x) + 1)^(gamma + 1))

}

#' @noRd
# Updated 08.08.2025
scad_derivative <- function(x, lambda, gamma = 3.7, ...)
{

  # Pre-compute components
  abs_x <- abs(x)
  sign_x <- x / abs_x
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
# Updated 12.01.2025
spot_derivative <- function(x, lambda, gamma = 3, ...)
{

  # Obtain exponent
  exponent <- exp(-abs(x) * 2^gamma)

  # Return lambdas
  return(lambda * 2^(gamma + 1) * exponent / (exponent + 1)^2)

}

#%%%%%%%%%%%%%%%%%%%%%%%%%
## Proximal Operators ----
#%%%%%%%%%%%%%%%%%%%%%%%%%

#' @noRd
# Updated 27.07.2025
dbridge_proximal <- function(x, lambda, gamma = 5, ...)
{
  return(l1_proximal(x, dbridge_derivative(x, lambda, gamma)))
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
# Updated 27.07.2025
lgp_proximal <- function(x, lambda, gamma = 5, ...)
{
  return(l1_proximal(x, lgp_derivative(x, lambda, gamma)))
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
# Updated 27.07.2025
pop_proximal <- function(x, lambda, gamma = 4, ...)
{
  return(l1_proximal(x, slam_derivative(x, lambda, gamma)))
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
# Updated 27.07.2025
spot_proximal <- function(x, lambda, gamma = 3, ...)
{
  return(l1_proximal(x, spot_derivative(x, lambda, gamma)))
}