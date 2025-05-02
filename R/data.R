#' Simulated Dataset for TMLE Examples
#'
#' A dataset containing simulated variables for demonstrating Targeted Maximum Likelihood
#' Estimation (TMLE) methods. The data includes a binary outcome, binary treatment,
#' various covariates, and a known true treatment effect (0.25) to evaluate
#' the performance of different TMLE methods.
#'
#' @format A data frame with 5000 rows (1000 observations × 5 simulations) and 8 variables:
#' \describe{
#'   \item{Y}{Binary outcome variable (0/1)}
#'   \item{A}{Binary treatment variable (0/1)}
#'   \item{X1}{Continuous covariate (standard normal distribution)}
#'   \item{X2}{Continuous covariate (standard normal distribution)}
#'   \item{X3}{Binary covariate (0/1 with 0.5 probability)}
#'   \item{X4}{Categorical covariate (factor with 3 levels)}
#'   \item{sim_id}{Simulation ID (1-5) for grouping observations}
#'   \item{true_effect}{True average treatment effect (0.25)}
#' }
#' @source Simulated data with data-generating process that ensures a true average 
#'   treatment effect of 0.25. The treatment assignment depends on covariates, 
#'   creating confounding that TMLE methods need to adjust for.
"tmle_example_data"

#' Create Example Dataset for TMLE Analysis
#'
#' Creates a simulated dataset with a binary outcome, binary treatment, and
#' covariates for demonstrating TMLE methods. The data-generating process
#' ensures a known true average treatment effect of 0.25, with confounding
#' introduced through covariate-dependent treatment assignment.
#'
#' @param n Number of observations per simulation
#' @param n_sims Number of simulation datasets to generate
#'
#' @return A data frame containing the simulated data with variables:
#'   \itemize{
#'     \item Y: Binary outcome variable (0/1)
#'     \item A: Binary treatment variable (0/1)
#'     \item X1-X4: Covariates (continuous and categorical)
#'     \item sim_id: Simulation ID
#'     \item true_effect: True treatment effect (0.25)
#'   }
#' @examples
#' # Generate a small example dataset
#' example_data <- create_example_data(n = 100, n_sims = 2)
#' head(example_data)
#' 
#' # Check the distribution of outcome by treatment
#' table(example_data$A, example_data$Y)
#' @export
create_example_data <- function(n = 1000, n_sims = 5) {
  set.seed(123)
  
  # Create data frame to store all simulations
  all_data <- data.frame()
  
  for (sim in 1:n_sims) {
    # Generate covariates
    X1 <- rnorm(n)
    X2 <- rnorm(n)
    X3 <- rbinom(n, 1, 0.5)
    X4 <- factor(sample(1:3, n, replace = TRUE))
    
    # Generate treatment based on covariates
    p_A <- plogis(-0.5 + 0.2*X1 - 0.1*X2 + 0.3*X3 + 0.2*(as.numeric(X4) - 2))
    A <- rbinom(n, 1, p_A)
    
    # Generate outcome based on treatment and covariates
    true_effect <- 0.25  # True ATE
    p_Y <- plogis(-1 + true_effect*A + 0.3*X1 + 0.2*X2 - 0.1*X3 + 0.15*(as.numeric(X4) - 2) + 0.1*A*X1)
    Y <- rbinom(n, 1, p_Y)
    
    # Create data frame for this simulation
    sim_data <- data.frame(
      Y = Y,
      A = A,
      X1 = X1,
      X2 = X2,
      X3 = X3,
      X4 = X4,
      sim_id = sim,
      true_effect = true_effect
    )
    
    # Append to all data
    all_data <- rbind(all_data, sim_data)
  }
  
  return(all_data)
}

# Run this in console! If you are using this for the first time!
# # Generate and save the example data (only do this once)
# tmle_example_data <- create_example_data()
# usethis::use_data(tmle_example_data, overwrite = TRUE)
# 
# # Update documentation
# devtools::document()