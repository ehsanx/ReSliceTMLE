#' Prepare data for TMLE analysis
#'
#' This function prepares a dataset for TMLE analysis by handling missing values,
#' converting variable types, and extracting outcome, treatment, and covariates.
#'
#' @param data A data frame containing the analysis variables
#' @param outcome_var Name of the outcome variable
#' @param treatment_var Name of the treatment variable
#' @param covariates Names of covariates to include in the analysis
#' @param binary_outcome Logical indicating if outcome is binary
#' @param missing_handling Strategy for handling missing values ("complete_case", "impute_mean", or "impute_model")
#'
#' @return A list containing prepared Y, A, and W components
#' @export
prepare_tmle_data <- function(data,
                              outcome_var,
                              treatment_var,
                              covariates,
                              binary_outcome = TRUE,
                              missing_handling = "complete_case") {

  # Check inputs
  if (!is.data.frame(data)) {
    stop("data must be a data frame")
  }

  if (!outcome_var %in% names(data)) {
    stop(paste("Outcome variable", outcome_var, "not found in data"))
  }

  if (!treatment_var %in% names(data)) {
    stop(paste("Treatment variable", treatment_var, "not found in data"))
  }

  missing_covars <- setdiff(covariates, names(data))
  if (length(missing_covars) > 0) {
    stop(paste("The following covariates are not found in data:",
               paste(missing_covars, collapse = ", ")))
  }

  # Extract needed variables
  needed_vars <- c(outcome_var, treatment_var, covariates)
  analysis_data <- data[, needed_vars, drop = FALSE]

  # Handle missing values
  if (missing_handling == "complete_case") {
    # Drop rows with any missing values
    complete_rows <- complete.cases(analysis_data)
    if (sum(!complete_rows) > 0) {
      message(paste("Removing", sum(!complete_rows), "rows with missing values"))
      analysis_data <- analysis_data[complete_rows, ]
    }
  } else if (missing_handling == "impute_mean") {
    # Impute missing values with column means
    for (col in names(analysis_data)) {
      if (is.numeric(analysis_data[[col]])) {
        missing <- is.na(analysis_data[[col]])
        if (any(missing)) {
          mean_val <- mean(analysis_data[[col]], na.rm = TRUE)
          analysis_data[missing, col] <- mean_val
          message(paste("Imputed", sum(missing), "missing values in", col, "with mean"))
        }
      } else if (is.factor(analysis_data[[col]]) || is.character(analysis_data[[col]])) {
        missing <- is.na(analysis_data[[col]])
        if (any(missing)) {
          mode_val <- names(sort(table(analysis_data[[col]]), decreasing = TRUE))[1]
          analysis_data[missing, col] <- mode_val
          message(paste("Imputed", sum(missing), "missing values in", col, "with mode"))
        }
      }
    }
  } else if (missing_handling == "impute_model") {
    # TODO: This would be implemented with mice or another imputation package
    warning("Model-based imputation not yet implemented!!")
    complete_rows <- complete.cases(analysis_data)
    if (sum(!complete_rows) > 0) {
      message(paste("Falling back to complete case analysis, removing",
                    sum(!complete_rows), "rows with missing values"))
      analysis_data <- analysis_data[complete_rows, ]
    }
  } else {
    stop("Invalid missing_handling method. Use 'complete_case', 'impute_mean', or 'impute_model'")
  }

  # Extract Y, A, and W
  Y <- analysis_data[[outcome_var]]
  A <- analysis_data[[treatment_var]]
  W <- analysis_data[, covariates, drop = FALSE]

  # Convert to numeric if needed
  if (!is.numeric(Y)) {
    if (binary_outcome) {
      Y <- as.numeric(as.factor(Y)) - 1
      if (!all(Y %in% c(0, 1))) {
        stop("Binary outcome must have exactly 2 levels")
      }
    } else {
      Y <- as.numeric(Y)
    }
  }

  if (!is.numeric(A)) {
    A <- as.numeric(as.factor(A)) - 1
    if (!all(A %in% c(0, 1))) {
      stop("Treatment must have exactly 2 levels")
    }
  }

  # Convert factor variables in W to dummy variables
  # W_processed <- W
  # for (col in names(W)) {
  #   if (is.factor(W[[col]]) || is.character(W[[col]])) {
  #     # Convert to dummy variables
  #     dummies <- stats::model.matrix(~ 0 + .data[[col]], data = W)
  #     colnames(dummies) <- paste0(col, "_", levels(as.factor(W[[col]])))
  #
  #     # Remove original column
  #     W_processed[[col]] <- NULL
  #
  #     # Add dummy columns
  #     for (dummy_col in colnames(dummies)) {
  #       W_processed[[dummy_col]] <- dummies[, dummy_col]
  #     }
  #   }
  # }

  return(list(
    Y = Y,
    A = A,
    W = W,
    original_data = analysis_data
  ))
}

#' Format TMLE results for easy reading and export
#'
#' @param results List of TMLE results from a TMLE function
#' @param confint Logical indicating whether to include confidence intervals
#' @param alpha Significance level for confidence intervals
#'
#' @return A formatted data frame with results
#' @export
format_results <- function(results, confint = TRUE, alpha = 0.05) {
  if (is.null(results)) {
    return(NULL)
  }

  # Initialize an empty tibble to store formatted results
  formatted <- dplyr::tibble()

  # Check if results contains a data frame of multiple results
  if ("results" %in% names(results) && is.data.frame(results$results) && nrow(results$results) > 0) {
    res_df <- results$results

    # For each row (result) in the data frame
    for (i in 1:nrow(res_df)) {
      row <- res_df[i, ]

      # Create formatted results for this row
      row_formatted <- dplyr::tibble(
        Estimand = c("E[Y(1)]", "E[Y(0)]", "ATE"),
        Estimate = c(row$r1, row$r0, row$rd)
      )

      if (confint) {
        z_score <- stats::qnorm(1 - alpha/2)

        # Use variance-adjusted estimates if available
        if ("mv1" %in% names(row) && !is.na(row$mv1[1])) {
          se <- sqrt(c(row$mv1, row$mv0, row$mvd))
        } else {
          se <- sqrt(c(row$v1, row$v0, row$vd))
        }

        row_formatted$SE <- se
        row_formatted$CI_Lower <- row_formatted$Estimate - z_score * se
        row_formatted$CI_Upper <- row_formatted$Estimate + z_score * se
        row_formatted$p_value <- 2 * stats::pnorm(-abs(row_formatted$Estimate / se))
      }

      # Add method if it exists
      if ("method" %in% names(row)) {
        row_formatted$Method <- rep(row$method, nrow(row_formatted))
      }

      # Add sim_id if it exists
      if ("sim_id" %in% names(row)) {
        row_formatted$sim_id <- rep(row$sim_id, nrow(row_formatted))
      }

      # Append to the formatted results
      formatted <- dplyr::bind_rows(formatted, row_formatted)
    }

    return(formatted)
  }
  # Handle original case of single result set
  else if ("results" %in% names(results)) {
    # For single TMLE runs
    res <- results$results

    # Create formatted results tibble
    formatted <- dplyr::tibble(
      Estimand = c("E[Y(1)]", "E[Y(0)]", "ATE"),
      Estimate = c(res$r1, res$r0, res$rd)
    )

    if (confint) {
      z_score <- stats::qnorm(1 - alpha/2)

      # Use variance-adjusted estimates if available
      if ("mv1" %in% names(res)) {
        se <- sqrt(c(res$mv1, res$mv0, res$mvd))
      } else {
        se <- sqrt(c(res$v1, res$v0, res$vd))
      }

      formatted$SE <- se
      formatted$CI_Lower <- formatted$Estimate - z_score * se
      formatted$CI_Upper <- formatted$Estimate + z_score * se
      formatted$p_value <- 2 * stats::pnorm(-abs(formatted$Estimate / se))
    }

    return(formatted)
  }
  else if ("aggregated_results" %in% names(results)) {
    # For multiple TMLE runs
    res <- results$aggregated_results

    # Create formatted results tibble
    formatted <- dplyr::tibble(
      Estimand = c("E[Y(1)]", "E[Y(0)]", "ATE"),
      Estimate = c(res$r1, res$r0, res$rd)
    )

    if (confint) {
      z_score <- stats::qnorm(1 - alpha/2)

      # Use variance-adjusted estimates if available
      if ("mv1" %in% names(res)) {
        se <- sqrt(c(res$mv1, res$mv0, res$mvd))
      } else {
        se <- sqrt(c(res$v1, res$v0, res$vd))
      }

      formatted$SE <- se
      formatted$CI_Lower <- formatted$Estimate - z_score * se
      formatted$CI_Upper <- formatted$Estimate + z_score * se
      formatted$p_value <- 2 * stats::pnorm(-abs(formatted$Estimate / se))
    }

    return(formatted)
  }
  else {
    stop("Unknown results format")
  }
}

#' Create SuperLearner library for TMLE
#'
#' @param include_glm Logical to include glm
#' @param include_gam Logical to include gam
#' @param include_rf Logical to include random forest
#' @param include_gbm Logical to include gradient boosting
#' @param include_nnet Logical to include neural network
#' @param include_mean Logical to include mean (intercept-only model)
#'
#' @return A character vector of SuperLearner library names
#' @export
create_sl_library <- function(
    include_glm = TRUE,
    include_gam = FALSE,
    include_rf = FALSE,
    include_gbm = FALSE,
    include_nnet = FALSE,
    include_mean = TRUE
) {
  library <- c()

  if (include_glm) library <- c(library, "SL.glm")
  if (include_gam) library <- c(library, "SL.gam")
  if (include_rf) library <- c(library, "SL.randomForest")
  if (include_gbm) library <- c(library, "SL.gbm")
  if (include_nnet) library <- c(library, "SL.nnet")
  if (include_mean) library <- c(library, "SL.mean")

  if (length(library) == 0) {
    warning("No SuperLearner algorithms selected, defaulting to SL.glm and SL.mean")
    library <- c("SL.glm", "SL.mean")
  }

  return(library)
}

#' Split Data for Simulation Studies
#'
#' @param data A data frame containing simulation data
#' @param sim_id_col Name of the column containing simulation IDs
#'
#' @return A list of data frames, one for each simulation ID
#' @export
split_simulation_data <- function(data, sim_id_col = "sim_id") {
  if (!sim_id_col %in% names(data)) {
    stop(paste("Simulation ID column", sim_id_col, "not found in data"))
  }

  # Get unique simulation IDs
  sim_ids <- unique(data[[sim_id_col]])

  # Split data by simulation ID
  split_data <- lapply(sim_ids, function(id) {
    data[data[[sim_id_col]] == id, ]
  })

  names(split_data) <- paste0("sim_", sim_ids)

  return(split_data)
}

#' Format TMLE results with simulation ID
#'
#' Given a list of TMLE results where each element corresponds to a different
#' simulation ID, this function returns a data frame with columns r1, r0, rd, v1, v0, vd,
#' and sim_id.
#'
#' @param results_list A list where each element contains results from a TMLE function
#' @param sim_ids Optional vector of simulation IDs. If NULL, uses numeric indices of the list elements.
#'
#' @return A tibble with columns r1, r0, rd, v1, v0, vd, and sim_id
#' @export
format_results_with_sim_id <- function(results_list, sim_ids = NULL) {
  # Create an empty data frame to store results
  combined_results <- dplyr::tibble()

  # If sim_ids is not provided, use the indices of results_list
  if (is.null(sim_ids)) {
    sim_ids <- 1:length(results_list)
  } else if (length(sim_ids) != length(results_list)) {
    warning("Length of sim_ids doesn't match length of results_list. Using indices instead.")
    sim_ids <- 1:length(results_list)
  }

  # Process each result in the list
  for (i in 1:length(results_list)) {
    result <- results_list[[i]]

    # Skip NULL results
    if (is.null(result)) {
      warning(paste("Skipping NULL result for simulation ID", sim_ids[i]))
      next
    }

    # Extract the relevant result component
    if ("results" %in% names(result)) {
      # For single TMLE runs
      res <- result$results
    } else if ("aggregated_results" %in% names(result)) {
      # For multiple TMLE runs
      res <- as.data.frame(t(result$aggregated_results))
    } else {
      warning(paste("Unknown result format for simulation ID", sim_ids[i]))
      next
    }

    # Add simulation ID
    res$sim_id <- sim_ids[i]

    # Append to combined results
    combined_results <- dplyr::bind_rows(combined_results, res)
  }

  return(combined_results)
}

#' Combine TMLE Results from Multiple Simulations
#'
#' @param results_list List of results from multiple simulations
#' @param true_values Optional true parameter values for coverage calculation
#'
#' @return A summary data frame of simulation results
#' @export
combine_simulation_results <- function(results_list, true_values = NULL) {
  # Check if we have any results
  if (length(results_list) == 0) {
    return(NULL)
  }

  # Extract results
  results <- dplyr::bind_rows(lapply(results_list, function(r) {
    if (is.null(r)) {
      return(NULL)
    }

    # For data frame results (usually from our run_tmle_analysis function)
    if (is.data.frame(r)) {
      return(r)
    }

    # For results with "results" component
    if ("results" %in% names(r)) {
      return(r$results)
    }

    # For results with "aggregated_results" component
    if ("aggregated_results" %in% names(r)) {
      return(as.data.frame(t(r$aggregated_results)))
    }

    # Return NULL if we can't extract results
    return(NULL)
  }), .id = "sim_id")

  # If no valid results, return NULL
  if (nrow(results) == 0) {
    return(NULL)
  }

  # Check if results has the required columns
  required_cols <- c("r1", "r0", "rd", "v1", "v0", "vd")
  missing_cols <- setdiff(required_cols, names(results))

  if (length(missing_cols) > 0) {
    warning(paste("Missing required columns in results:",
                  paste(missing_cols, collapse = ", "),
                  "Cannot calculate summary statistics."))
    return(list(
      individual_results = results,
      summary = NULL
    ))
  }

  # Calculate summary statistics
  # Group by method if it exists
  if ("method" %in% names(results)) {
    summary_stats <- results %>%
      dplyr::group_by(method) %>%
      dplyr::summarize(
        r1_mean = mean(r1, na.rm = TRUE),
        r1_sd = sd(r1, na.rm = TRUE),
        r0_mean = mean(r0, na.rm = TRUE),
        r0_sd = sd(r0, na.rm = TRUE),
        rd_mean = mean(rd, na.rm = TRUE),
        rd_sd = sd(rd, na.rm = TRUE),
        v1_mean = mean(v1, na.rm = TRUE),
        v0_mean = mean(v0, na.rm = TRUE),
        vd_mean = mean(vd, na.rm = TRUE),
        n_converged = sum(!is.na(rd))
      )
  } else {
    summary_stats <- results %>%
      dplyr::summarize(
        r1_mean = mean(r1, na.rm = TRUE),
        r1_sd = sd(r1, na.rm = TRUE),
        r0_mean = mean(r0, na.rm = TRUE),
        r0_sd = sd(r0, na.rm = TRUE),
        rd_mean = mean(rd, na.rm = TRUE),
        rd_sd = sd(rd, na.rm = TRUE),
        v1_mean = mean(v1, na.rm = TRUE),
        v0_mean = mean(v0, na.rm = TRUE),
        vd_mean = mean(vd, na.rm = TRUE),
        n_converged = sum(!is.na(rd))
      )
  }

  # Add coverage statistics if true values provided
  if (!is.null(true_values) && all(c("r1", "r0", "rd") %in% names(true_values))) {
    # Function to calculate 95% CI coverage
    calc_coverage <- function(estimate, var, true_value) {
      z <- 1.96
      lower <- estimate - z * sqrt(var)
      upper <- estimate + z * sqrt(var)
      covered <- lower <= true_value & upper >= true_value
      mean(covered, na.rm = TRUE)
    }

    # Check if we have variance-adjusted estimates
    has_mv <- all(c("mv1", "mv0", "mvd") %in% names(results))

    if (has_mv) {
      # Calculate coverage with variance-adjusted estimates if available
      if ("method" %in% names(results)) {
        coverage <- results %>%
          dplyr::group_by(method) %>%
          dplyr::summarize(
            r1_coverage = calc_coverage(r1, mv1, true_values$r1),
            r0_coverage = calc_coverage(r0, mv0, true_values$r0),
            rd_coverage = calc_coverage(rd, mvd, true_values$rd)
          )
      } else {
        coverage <- results %>%
          dplyr::summarize(
            r1_coverage = calc_coverage(r1, mv1, true_values$r1),
            r0_coverage = calc_coverage(r0, mv0, true_values$r0),
            rd_coverage = calc_coverage(rd, mvd, true_values$rd)
          )
      }
    } else {
      # Calculate coverage with regular variance estimates
      if ("method" %in% names(results)) {
        coverage <- results %>%
          dplyr::group_by(method) %>%
          dplyr::summarize(
            r1_coverage = calc_coverage(r1, v1, true_values$r1),
            r0_coverage = calc_coverage(r0, v0, true_values$r0),
            rd_coverage = calc_coverage(rd, vd, true_values$rd)
          )
      } else {
        coverage <- results %>%
          dplyr::summarize(
            r1_coverage = calc_coverage(r1, v1, true_values$r1),
            r0_coverage = calc_coverage(r0, v0, true_values$r0),
            rd_coverage = calc_coverage(rd, vd, true_values$rd)
          )
      }
    }

    if ("method" %in% names(summary_stats)) {
      summary_stats <- dplyr::left_join(summary_stats, coverage, by = "method")
    } else {
      summary_stats <- dplyr::bind_cols(summary_stats, coverage)
    }
  }

  return(list(
    individual_results = results,
    summary = summary_stats
  ))
}
