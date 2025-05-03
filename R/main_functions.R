#' @importFrom dplyr bind_rows %>% as_tibble mutate group_by summarize
#' @importFrom tidyr nest unnest
#' @importFrom purrr map
#' @importFrom utils write.csv
#' @importFrom stats setNames
NULL

#' Run TMLE Analysis on a Single Dataset or Multiple Simulations
#'
#' This function provides a unified interface for applying any TMLE variant
#' to a dataset or multiple simulation datasets, with options for parallelization.
#'
#' @param data A data frame containing analysis variables
#' @param outcome_var Name of the outcome variable
#' @param treatment_var Name of the treatment variable
#' @param covariates Vector of covariate names
#' @param tmle_variant TMLE variant to use (one of 'vanilla', 'cvq', 'cvq_multiple',
#'                    'fullcv', 'fullcv_multiple', 'singlecrossfit', 'doublecrossfit')
#' @param Q.SL.library SuperLearner library for outcome model
#' @param g.SL.library SuperLearner library for treatment model
#' @param V Number of cross-validation folds (for cross-validated methods)
#' @param family Distribution family ("binomial" or "gaussian")
#' @param num_repetitions Number of repetitions for multiple/crossfit methods
#' @param sim_id_col Name of the column containing simulation IDs (NULL for non-simulation data)
#' @param true_value True parameter value(s) for coverage calculations (scalar ATE or vector of r1, r0, rd)
#' @param parallel_strategy Parallelization strategy ('none', 'local', or 'slurm')
#' @param n_cores Number of cores to use for local parallelization
#' @param slurm_options List of SLURM options for SLURM parallelization
#' @param missing_handling Strategy for handling missing values
#' @param binary_outcome Logical indicating if outcome is binary
#' @param save_results Logical indicating whether to save results to file
#' @param results_dir Directory to save results (if save_results = TRUE)
#' @param ... Additional arguments passed to TMLE functions
#'
#' @return List containing TMLE analysis results with components:
#' \itemize{
#'   \item raw_final: Tibble with raw final results (r1, r0, rd, v1, v0, vd, mv1, mv0, mvd, method)
#'   \item raw_intermediate: Tibble with raw intermediate results for repetition-based methods
#'   \item interpreted_results: Tibble with formatted results (Estimand, Estimate, SE, CI, p-value)
#' }
#' @export
run_tmle_analysis <- function(
    data,
    outcome_var,
    treatment_var,
    covariates,
    tmle_variant = c("vanilla", "cvq", "cvq_multiple", "fullcv",
                     "fullcv_multiple", "singlecrossfit", "doublecrossfit"),
    Q.SL.library = NULL,
    g.SL.library = NULL,
    V = 5,
    family = "binomial",
    num_repetitions = 100,
    sim_id_col = NULL,
    true_value = NULL,
    parallel_strategy = c("none", "local", "slurm"),
    n_cores = 2,
    slurm_options = list(),
    missing_handling = "complete_case",
    binary_outcome = TRUE,
    save_results = FALSE,
    results_dir = "./results",
    ...
) {
  # Match arguments
  tmle_variant <- match.arg(tmle_variant)
  parallel_strategy <- match.arg(parallel_strategy)
  
  # Set up default SuperLearner libraries if not provided
  if (is.null(Q.SL.library)) {
    Q.SL.library <- create_sl_library()
  }
  if (is.null(g.SL.library)) {
    g.SL.library <- create_sl_library()
  }
  
  # Process true_value
  true_values_list <- NULL
  if (!is.null(true_value)) {
    if (length(true_value) == 1) {
      # Single ATE value provided
      true_values_list <- list(r1 = NA, r0 = NA, rd = true_value)
    } else if (length(true_value) == 3) {
      # Full r1, r0, rd values provided
      true_values_list <- list(r1 = true_value[1], r0 = true_value[2], rd = true_value[3])
    } else {
      warning("Invalid true_value format. Expected either a single ATE value or a vector of 3 values (r1, r0, rd).")
    }
  }
  
  # Determine if we're dealing with simulation data
  is_simulation <- !is.null(sim_id_col) && sim_id_col %in% names(data)
  
  # Initialize results containers
  raw_final_results <- NULL
  raw_intermediate_results <- NULL
  
  # Handle simulation data
  if (is_simulation) {
    # Split data by simulation ID
    split_data_list <- split_simulation_data(data, sim_id_col)
    sim_ids <- as.numeric(gsub("sim_", "", names(split_data_list)))
    
    # Process each simulation dataset
    sim_results_list <- list()
    
    message(paste("Processing", length(split_data_list), "simulation datasets..."))
    
    for (i in seq_along(split_data_list)) {
      sim_id <- sim_ids[i]
      sim_data <- split_data_list[[i]]
      
      message(paste("Processing simulation ID:", sim_id, "(", i, "of", length(split_data_list), ")"))
      
      # Run TMLE on this simulation dataset
      result <- run_single_tmle_analysis(
        data = sim_data,
        outcome_var = outcome_var,
        treatment_var = treatment_var,
        covariates = covariates,
        tmle_variant = tmle_variant,
        Q.SL.library = Q.SL.library,
        g.SL.library = g.SL.library,
        V = V,
        family = family,
        num_repetitions = num_repetitions,
        missing_handling = missing_handling,
        binary_outcome = binary_outcome,
        ...
      )
      
      # Add simulation ID to results
      if (!is.null(result$raw_final)) {
        result$raw_final$sim_id <- sim_id
      }
      
      if (!is.null(result$raw_intermediate)) {
        result$raw_intermediate$sim_id <- sim_id
      }
      
      sim_results_list[[i]] <- result
    }
    
    # Combine results from all simulations
    raw_final_results <- dplyr::bind_rows(lapply(sim_results_list, function(x) x$raw_final))
    
    # Only combine intermediate results if they exist
    has_intermediate <- any(sapply(sim_results_list, function(x) !is.null(x$raw_intermediate)))
    if (has_intermediate) {
      raw_intermediate_results <- dplyr::bind_rows(
        lapply(sim_results_list, function(x) x$raw_intermediate)
      )
    }
    
    # If we have true values, calculate coverage
    if (!is.null(true_values_list)) {
      coverage_results <- combine_simulation_results(
        lapply(sim_results_list, function(x) x$raw_final),
        true_values = true_values_list
      )
      
      # Add coverage to results
      combined_summary <- coverage_results$summary
    } else {
      # Just combine results without coverage
      combined_summary <- NULL
    }
    
  } else {
    # Single dataset analysis
    result <- run_single_tmle_analysis(
      data = data,
      outcome_var = outcome_var,
      treatment_var = treatment_var,
      covariates = covariates,
      tmle_variant = tmle_variant,
      Q.SL.library = Q.SL.library,
      g.SL.library = g.SL.library,
      V = V,
      family = family,
      num_repetitions = num_repetitions,
      missing_handling = missing_handling,
      binary_outcome = binary_outcome,
      ...
    )
    
    raw_final_results <- result$raw_final
    raw_intermediate_results <- result$raw_intermediate
    combined_summary <- NULL
  }
  
  # Create formatted/interpreted results from raw_final_results
  if (!is.null(raw_final_results)) {
    # Calculate confidence intervals and format results
    interpreted_results <- format_results(list(results = raw_final_results), confint = TRUE)
  } else {
    interpreted_results <- NULL
  }
  
  # Save results if requested
  if (save_results) {
    # Create results directory if it doesn't exist
    if (!dir.exists(results_dir)) {
      dir.create(results_dir, recursive = TRUE)
    }
    
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    
    # Save raw final results
    if (!is.null(raw_final_results)) {
      raw_final_file <- file.path(results_dir, paste0("tmle_raw_final_", timestamp, ".csv"))
      utils::write.csv(raw_final_results, raw_final_file, row.names = FALSE)
      message(paste("Raw final results saved to:", raw_final_file))
    }
    
    # Save raw intermediate results if they exist
    if (!is.null(raw_intermediate_results)) {
      raw_intermediate_file <- file.path(results_dir, paste0("tmle_raw_intermediate_", timestamp, ".csv"))
      utils::write.csv(raw_intermediate_results, raw_intermediate_file, row.names = FALSE)
      message(paste("Raw intermediate results saved to:", raw_intermediate_file))
    }
    
    # Save interpreted results
    if (!is.null(interpreted_results)) {
      interpreted_file <- file.path(results_dir, paste0("tmle_interpreted_", timestamp, ".csv"))
      utils::write.csv(interpreted_results, interpreted_file, row.names = FALSE)
      message(paste("Interpreted results saved to:", interpreted_file))
    }
    
    # Save combined summary if it exists
    if (!is.null(combined_summary)) {
      summary_file <- file.path(results_dir, paste0("tmle_summary_", timestamp, ".csv"))
      utils::write.csv(combined_summary, summary_file, row.names = FALSE)
      message(paste("Summary results saved to:", summary_file))
    }
  }
  
  # Return all results
  return(list(
    raw_final = raw_final_results,
    raw_intermediate = raw_intermediate_results,
    interpreted_results = interpreted_results,
    summary = combined_summary
  ))
}

#' Run TMLE Analysis on a Single Dataset
#'
#' Helper function to run TMLE on a single dataset. Used internally by run_tmle_analysis.
#'
#' @param data A data frame containing analysis variables
#' @param outcome_var Name of the outcome variable
#' @param treatment_var Name of the treatment variable
#' @param covariates Vector of covariate names
#' @param tmle_variant TMLE variant to use
#' @param Q.SL.library SuperLearner library for outcome model
#' @param g.SL.library SuperLearner library for treatment model
#' @param V Number of cross-validation folds
#' @param family Distribution family
#' @param num_repetitions Number of repetitions for multiple methods
#' @param missing_handling Strategy for handling missing values
#' @param binary_outcome Logical indicating if outcome is binary
#' @param ... Additional arguments passed to TMLE functions
#'
#' @return List containing TMLE results for a single dataset
run_single_tmle_analysis <- function(
    data,
    outcome_var,
    treatment_var,
    covariates,
    tmle_variant,
    Q.SL.library,
    g.SL.library,
    V,
    family,
    num_repetitions,
    missing_handling,
    binary_outcome,
    ...
) {
  # Prepare data for TMLE
  prepared_data <- prepare_tmle_data(
    data = data,
    outcome_var = outcome_var,
    treatment_var = treatment_var,
    covariates = covariates,
    binary_outcome = binary_outcome,
    missing_handling = missing_handling
  )
  
  Y <- prepared_data$Y
  A <- prepared_data$A
  W <- prepared_data$W
  
  # Initialize results
  raw_final <- NULL
  raw_intermediate <- NULL
  
  # Run the appropriate TMLE variant
  if (tmle_variant == "vanilla") {
    result <- vanilla_tmle(
      Y = Y,
      A = A,
      W = W,
      Q.SL.library = Q.SL.library,
      g.SL.library = g.SL.library,
      family = family,
      ...
    )
    
    if (!is.null(result)) {
      raw_final <- result$results
      raw_final$method <- "vanilla"
      # Add NA for mv* columns
      raw_final$mv1 <- NA
      raw_final$mv0 <- NA
      raw_final$mvd <- NA
    }
    
  } else if (tmle_variant == "cvq") {
    result <- cvq_tmle(
      Y = Y,
      A = A,
      W = W,
      Q.SL.library = Q.SL.library,
      g.SL.library = g.SL.library,
      V.Q = V,
      V.g = V,
      family = family,
      ...
    )
    
    if (!is.null(result)) {
      raw_final <- result$results
      raw_final$method <- "cvq"
      # Add NA for mv* columns
      raw_final$mv1 <- NA
      raw_final$mv0 <- NA
      raw_final$mvd <- NA
    }
    
  } else if (tmle_variant == "cvq_multiple") {
    result <- cvq_tmle_multiple(
      Y = Y,
      A = A,
      W = W,
      Q.SL.library = Q.SL.library,
      g.SL.library = g.SL.library,
      V.Q = V,
      V.g = V,
      family = family,
      num_repetitions = num_repetitions,
      ...
    )
    
    if (!is.null(result)) {
      # Final results
      raw_final <- dplyr::as_tibble(t(result$aggregated_results)) %>%
        dplyr::mutate(method = "cvq_multiple")
      
      # Intermediate results (repetitions)
      raw_intermediate <- result$repetitions %>%
        dplyr::mutate(
          iteration = 1:n(),
          method = "cvq_multiple"
        )
    }
    
  } else if (tmle_variant == "fullcv") {
    result <- fullcv_tmle(
      Y = Y,
      A = A,
      W = W,
      Q.SL.library = Q.SL.library,
      g.SL.library = g.SL.library,
      V = V,
      family = family,
      ...
    )
    
    if (!is.null(result)) {
      raw_final <- result$results
      raw_final$method <- "fullcv"
      # Add NA for mv* columns
      raw_final$mv1 <- NA
      raw_final$mv0 <- NA
      raw_final$mvd <- NA
    }
    
  } else if (tmle_variant == "fullcv_multiple") {
    result <- fullcv_tmle_multiple(
      Y = Y,
      A = A,
      W = W,
      Q.SL.library = Q.SL.library,
      g.SL.library = g.SL.library,
      V = V,
      family = family,
      num_repetitions = num_repetitions,
      ...
    )
    
    if (!is.null(result)) {
      # Final results
      raw_final <- dplyr::as_tibble(t(result$aggregated_results)) %>%
        dplyr::mutate(method = "fullcv_multiple")
      
      # Intermediate results (repetitions)
      raw_intermediate <- result$repetitions %>%
        dplyr::mutate(
          iteration = 1:n(),
          method = "fullcv_multiple"
        )
    }
    
  } else if (tmle_variant == "singlecrossfit") {
    result <- singlecrossfit_tmle(
      Y = Y,
      A = A,
      W = W,
      Q.SL.library = Q.SL.library,
      g.SL.library = g.SL.library,
      family = family,
      num_repetitions = num_repetitions,
      ...
    )
    
    if (!is.null(result)) {
      # Final results
      raw_final <- dplyr::as_tibble(t(result$aggregated_results)) %>%
        dplyr::mutate(method = "singlecrossfit")
      
      # Intermediate results (repetitions)
      raw_intermediate <- result$repetitions %>%
        dplyr::mutate(
          iteration = 1:n(),
          method = "singlecrossfit"
        )
    }
    
  } else if (tmle_variant == "doublecrossfit") {
    result <- doublecrossfit_tmle(
      Y = Y,
      A = A,
      W = W,
      Q.SL.library = Q.SL.library,
      g.SL.library = g.SL.library,
      family = family,
      num_repetitions = num_repetitions,
      ...
    )
    
    if (!is.null(result)) {
      # Final results
      raw_final <- dplyr::as_tibble(t(result$aggregated_results)) %>%
        dplyr::mutate(method = "doublecrossfit")
      
      # Intermediate results (repetitions)
      raw_intermediate <- result$repetitions %>%
        dplyr::mutate(
          iteration = 1:n(),
          method = "doublecrossfit"
        )
    }
  }
  
  # Return results
  return(list(
    raw_final = raw_final,
    raw_intermediate = raw_intermediate
  ))
}

#' Run Multiple TMLE Variants on a Dataset
#'
#' This function applies multiple TMLE variants to the same dataset and
#' combines the results for easy comparison.
#'
#' @param data A data frame containing analysis variables
#' @param outcome_var Name of the outcome variable
#' @param treatment_var Name of the treatment variable
#' @param covariates Vector of covariate names
#' @param tmle_variants Vector of TMLE variants to use (from 'vanilla', 'cvq', etc.)
#' @param Q.SL.library SuperLearner library for outcome model
#' @param g.SL.library SuperLearner library for treatment model
#' @param V Number of cross-validation folds
#' @param family Distribution family
#' @param num_repetitions Number of repetitions for multiple methods
#' @param sim_id_col Name of the column containing simulation IDs
#' @param true_value True parameter value(s) for coverage calculations
#' @param parallel_strategy Parallelization strategy
#' @param n_cores Number of cores to use for local parallelization
#' @param slurm_options List of SLURM options for SLURM parallelization
#' @param missing_handling Strategy for handling missing values
#' @param binary_outcome Logical indicating if outcome is binary
#' @param save_results Logical indicating whether to save results to file
#' @param results_dir Directory to save results
#' @param ... Additional arguments passed to TMLE functions
#'
#' @return List containing combined TMLE analysis results from multiple variants
#' @export
run_multiple_tmle_variants <- function(
    data,
    outcome_var,
    treatment_var,
    covariates,
    tmle_variants = c("vanilla", "cvq_multiple", "fullcv_multiple", "doublecrossfit"),
    Q.SL.library = NULL,
    g.SL.library = NULL,
    V = 5,
    family = "binomial",
    num_repetitions = 100,
    sim_id_col = NULL,
    true_value = NULL,
    parallel_strategy = "none",
    n_cores = 2,
    slurm_options = list(),
    missing_handling = "complete_case",
    binary_outcome = TRUE,
    save_results = FALSE,
    results_dir = "./results",
    ...
) {
  # Validate TMLE variant choices
  valid_variants <- c("vanilla", "cvq", "cvq_multiple", "fullcv",
                      "fullcv_multiple", "singlecrossfit", "doublecrossfit")
  invalid_variants <- setdiff(tmle_variants, valid_variants)
  
  if (length(invalid_variants) > 0) {
    stop(paste("Invalid TMLE variants:", paste(invalid_variants, collapse = ", "),
               ". Valid options are:", paste(valid_variants, collapse = ", ")))
  }
  
  # Initialize result containers
  all_results <- list()
  all_raw_final <- NULL
  all_raw_intermediate <- NULL
  
  # Run each TMLE variant
  for (variant in tmle_variants) {
    message(paste("Running TMLE variant:", variant))
    
    result <- run_tmle_analysis(
      data = data,
      outcome_var = outcome_var,
      treatment_var = treatment_var,
      covariates = covariates,
      tmle_variant = variant,
      Q.SL.library = Q.SL.library,
      g.SL.library = g.SL.library,
      V = V,
      family = family,
      num_repetitions = num_repetitions,
      sim_id_col = sim_id_col,
      true_value = true_value,
      parallel_strategy = parallel_strategy,
      n_cores = n_cores,
      slurm_options = slurm_options,
      missing_handling = missing_handling,
      binary_outcome = binary_outcome,
      save_results = FALSE,  # We'll save the combined results later
      ...
    )
    
    all_results[[variant]] <- result
    all_raw_final <- dplyr::bind_rows(all_raw_final, result$raw_final)
    
    if (!is.null(result$raw_intermediate)) {
      all_raw_intermediate <- dplyr::bind_rows(all_raw_intermediate, result$raw_intermediate)
    }
  }
  
  # Create combined interpreted results
  combined_interpreted <- format_results(list(results = all_raw_final), confint = TRUE)
  
  # Save combined results if requested
  if (save_results) {
    # Create results directory if it doesn't exist
    if (!dir.exists(results_dir)) {
      dir.create(results_dir, recursive = TRUE)
    }
    
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    variant_suffix <- paste(tmle_variants, collapse = "_")
    
    # Save raw final results
    if (!is.null(all_raw_final)) {
      raw_final_file <- file.path(results_dir,
                                  paste0("tmle_combined_raw_final_", variant_suffix, "_", timestamp, ".csv"))
      utils::write.csv(all_raw_final, raw_final_file, row.names = FALSE)
      message(paste("Combined raw final results saved to:", raw_final_file))
    }
    
    # Save raw intermediate results if they exist
    if (!is.null(all_raw_intermediate)) {
      raw_intermediate_file <- file.path(results_dir,
                                         paste0("tmle_combined_raw_intermediate_", variant_suffix, "_", timestamp, ".csv"))
      utils::write.csv(all_raw_intermediate, raw_intermediate_file, row.names = FALSE)
      message(paste("Combined raw intermediate results saved to:", raw_intermediate_file))
    }
    
    # Save interpreted results
    if (!is.null(combined_interpreted)) {
      interpreted_file <- file.path(results_dir,
                                    paste0("tmle_combined_interpreted_", variant_suffix, "_", timestamp, ".csv"))
      utils::write.csv(combined_interpreted, interpreted_file, row.names = FALSE)
      message(paste("Combined interpreted results saved to:", interpreted_file))
    }
  }
  
  # Return combined results
  return(list(
    raw_final = all_raw_final,
    raw_intermediate = all_raw_intermediate,
    interpreted_results = combined_interpreted,
    individual_results = all_results
  ))
}