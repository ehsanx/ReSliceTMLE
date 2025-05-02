#' Run vanilla TMLE (Targeted Maximum Likelihood Estimation)
#'
#' @param Y Outcome variable vector (binary 0/1)
#' @param A Treatment variable vector (binary 0/1)
#' @param W Covariate matrix or data frame
#' @param Q.SL.library SuperLearner library for outcome model
#' @param g.SL.library SuperLearner library for treatment model
#' @param family Distribution family ("binomial" or "gaussian")
#' @param ... Additional arguments passed to tmle::tmle
#'
#' @return A list containing TMLE results
#' @export
vanilla_tmle <- function(Y, A, W,
                         Q.SL.library = c("SL.glm", "SL.mean"),
                         g.SL.library = c("SL.glm", "SL.mean"),
                         family = "binomial",
                         ...) {

  # Input validation
  validate_inputs(Y, A, W, family)

  # Run TMLE
  tryCatch({
    tmle_fit <- tmle::tmle(Y = Y,
                           A = A,
                           W = W,
                           Q.SL.library = Q.SL.library,
                           g.SL.library = g.SL.library,
                           family = family,
                           cvQinit = FALSE,
                           ...)

    # Format results
    results <- format_tmle_results(tmle_fit)

    return(list(
      fit = tmle_fit,
      results = results
    ))
  }, error = function(e) {
    warning(paste("Error in vanilla TMLE:", e$message))
    return(NULL)
  })
}

#' Run Cross-validated on Output Model Targeted Maximum Likelihood Estimation (CVq-TMLE)
#'
#' @param Y Outcome variable vector (binary 0/1)
#' @param A Treatment variable vector (binary 0/1)
#' @param W Covariate matrix or data frame
#' @param Q.SL.library SuperLearner library for outcome model
#' @param g.SL.library SuperLearner library for treatment model
#' @param V.Q Number of CV folds for outcome model
#' @param V.g Number of CV folds for treatment model
#' @param family Distribution family ("binomial" or "gaussian")
#' @param ... Additional arguments passed to tmle::tmle
#'
#' @return A list containing CVq-TMLE results
#' @export
cvq_tmle <- function(Y, A, W,
                     Q.SL.library = c("SL.glm", "SL.mean"),
                     g.SL.library = c("SL.glm", "SL.mean"),
                     V.Q = 5,
                     V.g = 5,
                     family = "binomial",
                     ...) {

  # Input validation
  validate_inputs(Y, A, W, family)

  # Run CV-TMLE
  tryCatch({
    tmle_fit <- tmle::tmle(Y = Y,
                           A = A,
                           W = W,
                           Q.SL.library = Q.SL.library,
                           g.SL.library = g.SL.library,
                           V.Q = V.Q,
                           V.g = V.g,
                           family = family,
                           cvQinit = TRUE,
                           ...)

    # Format results
    results <- format_tmle_results(tmle_fit)

    return(list(
      fit = tmle_fit,
      results = results
    ))
  }, error = function(e) {
    warning(paste("Error in CVq-TMLE:", e$message))
    return(NULL)
  })
}

#' Run Multiple CVq-TMLE with Repetition
#'
#' Runs CVq-TMLE multiple times and aggregates the results
#'
#' @param Y Outcome variable vector (binary 0/1)
#' @param A Treatment variable vector (binary 0/1)
#' @param W Covariate matrix or data frame
#' @param Q.SL.library SuperLearner library for outcome model
#' @param g.SL.library SuperLearner library for treatment model
#' @param V.Q Number of CV folds for outcome model
#' @param V.g Number of CV folds for treatment model
#' @param family Distribution family ("binomial" or "gaussian")
#' @param num_repetitions Number of times to repeat CV-TMLE
#' @param max_attempts Maximum number of attempts for valid results
#' @param ... Additional arguments passed to tmle::tmle
#'
#' @return A list containing aggregated CV-TMLE results
#' @export
cvq_tmle_multiple <- function(Y, A, W,
                              Q.SL.library = c("SL.glm", "SL.mean"),
                              g.SL.library = c("SL.glm", "SL.mean"),
                              V.Q = 5,
                              V.g = 5,
                              family = "binomial",
                              num_repetitions = 100,
                              max_attempts = 500,
                              ...) {

  # Input validation
  validate_inputs(Y, A, W, family)

  # Function to run a single CV-TMLE iteration
  run_iteration <- function() {
    result <- cvq_tmle(Y, A, W,
                       Q.SL.library, g.SL.library,
                       V.Q, V.g, family, ...)
    if (!is.null(result)) {
      return(result$results)
    }
    return(NULL)
  }

  # Initialize storage for valid results
  valid_results <- list()
  total_attempts <- 0

  # Keep trying until we get enough results or hit max attempts
  while (length(valid_results) < num_repetitions &&
         total_attempts < max_attempts) {

    # Run a single iteration
    result <- run_iteration()
    total_attempts <- total_attempts + 1

    # Store valid result
    if (!is.null(result)) {
      valid_results[[length(valid_results) + 1]] <- result
    }

    # Print progress every 10 attempts
    if (total_attempts %% 10 == 0) {
      message(paste(length(valid_results),
                    "valid results after",
                    total_attempts, "attempts"))
    }
  }

  # Handle case where we didn't get enough valid results
  if (length(valid_results) < num_repetitions) {
    warning(paste("Only obtained", length(valid_results),
                  "valid results after", total_attempts, "attempts"))
    if (length(valid_results) == 0) {
      return(NULL)
    }
  }

  # Combine results
  runs <- dplyr::bind_rows(valid_results[1:min(length(valid_results), num_repetitions)])

  # Calculate medians
  medians <- apply(runs, 2, stats::median)

  # Calculate variance-adjusted estimates
  runs <- runs %>%
    dplyr::mutate(
      mv1 = v1 + (r1 - medians["r1"])^2,
      mv0 = v0 + (r0 - medians["r0"])^2,
      mvd = vd + (rd - medians["rd"])^2
    )

  # Final aggregated results
  results <- apply(runs, 2, stats::median)

  return(list(
    repetitions = runs,
    aggregated_results = results
  ))
}

#' Run Full Cross-validated TMLE
#'
#' Implements full cross-validated TMLE where initial Q and g are fitted on training
#' folds and then targeted on test folds. This avoids overfitting bias.
#'
#' @param Y Outcome variable vector (binary 0/1)
#' @param A Treatment variable vector (binary 0/1)
#' @param W Covariate matrix or data frame
#' @param Q.SL.library SuperLearner library for outcome model
#' @param g.SL.library SuperLearner library for treatment model
#' @param V Number of CV folds
#' @param family Distribution family ("binomial" or "gaussian")
#' @param ... Additional arguments
#'
#' @return A list containing full CV-TMLE results
#' @export
fullcv_tmle <- function(Y, A, W,
                        Q.SL.library = c("SL.glm", "SL.mean"),
                        g.SL.library = c("SL.glm", "SL.mean"),
                        V = 5,
                        family = "binomial",
                        ...) {

  # Input validation
  validate_inputs(Y, A, W, family)

  # Ensure W is a data frame
  if (is.matrix(W)) {
    W <- as.data.frame(W)
  }

  # Create full data frame
  data <- data.frame(Y = Y, A = A, W)

  tryCatch({
    # Create V random folds
    folds <- sample(rep(1:V, length.out = length(Y)))

    # Prepare storage for fold-specific estimates
    fold_estimates <- numeric(V)
    fold_ey1 <- numeric(V)
    fold_ey0 <- numeric(V)
    all_ic <- numeric(0)  # to store all influence curves

    # Loop through each fold
    for (k in 1:V) {
      # Split into training and test sets
      train_idx <- folds != k
      train_data <- data[train_idx, ]
      test_data <- data[!train_idx, ]

      # Extract components
      train_Y <- train_data$Y
      train_A <- train_data$A
      train_W <- train_data[, !(names(train_data) %in% c("Y", "A")), drop = FALSE]

      test_Y <- test_data$Y
      test_A <- test_data$A
      test_W <- test_data[, !(names(test_data) %in% c("Y", "A")), drop = FALSE]

      # Fit Q-model on training data
      # Include treatment in X for Q-model
      Q_covariates <- cbind(A = train_A, train_W)

      Q_fit <- SuperLearner::SuperLearner(
        Y = train_Y,
        X = Q_covariates,
        family = if(family == "binomial") binomial() else gaussian(),
        SL.library = Q.SL.library,
        ...
      )

      # Fit g-model on training data
      g_fit <- SuperLearner::SuperLearner(
        Y = train_A,
        X = train_W,
        family = binomial(),
        SL.library = g.SL.library,
        ...
      )

      # Create counterfactual datasets for test data
      test_data_A0 <- test_data
      test_data_A1 <- test_data
      test_data_A0$A <- 0
      test_data_A1$A <- 1

      # Extract counterfactual covariates
      test_A0_cov <- cbind(A = 0, test_W)
      test_A1_cov <- cbind(A = 1, test_W)

      # Make predictions for both potential outcomes
      Q0_pred <- predict(Q_fit, newdata = test_A0_cov)$pred
      Q1_pred <- predict(Q_fit, newdata = test_A1_cov)$pred

      # Bound predictions to avoid extreme values
      if (family == "binomial") {
        Q0_pred <- pmax(pmin(Q0_pred, 0.99), 0.01)
        Q1_pred <- pmax(pmin(Q1_pred, 0.99), 0.01)
      }

      # Get g predictions
      g_pred <- predict(g_fit, newdata = test_W)$pred
      g_pred <- pmax(pmin(g_pred, 0.99), 0.01)

      # Create Q matrix with counterfactual predictions
      Q_matrix <- cbind(Q0_pred, Q1_pred)

      # Run TMLE on test data
      tmle_fit <- tmle::tmle(
        Y = test_Y,
        A = test_A,
        W = test_W,
        Q = Q_matrix,
        g1W = g_pred,
        family = family
      )

      # Store estimates
      fold_estimates[k] <- tmle_fit$estimates$ATE$psi
      fold_ey1[k] <- tmle_fit$estimates$EY1$psi
      fold_ey0[k] <- tmle_fit$estimates$EY0$psi

      # Compute influence curves for this fold
      Qstar <- tmle_fit$Qstar
      Q0star <- Qstar[, 1]
      Q1star <- Qstar[, 2]
      Qstar_obs <- ifelse(test_A == 1, Q1star, Q0star)

      h_a <- test_A/g_pred - (1-test_A)/(1-g_pred)
      resid <- test_Y - Qstar_obs
      fold_ic <- h_a * resid + (Q1star - Q0star) - tmle_fit$estimates$ATE$psi
      all_ic <- c(all_ic, fold_ic)
    }

    # Compute final estimates
    final_est <- mean(fold_estimates)
    ey1_est <- mean(fold_ey1)
    ey0_est <- mean(fold_ey0)

    # Variance estimation based on influence curves
    var_est <- var(all_ic) / length(all_ic)

    # Format results
    results <- dplyr::tibble(
      r1 = ey1_est,
      r0 = ey0_est,
      rd = final_est,
      v1 = var_est,
      v0 = var_est,
      vd = var_est
    )

    return(list(
      fit = list(fold_estimates = fold_estimates,
                 fold_ey1 = fold_ey1,
                 fold_ey0 = fold_ey0,
                 influence_curves = all_ic),
      results = results
    ))
  }, error = function(e) {
    warning(paste("Error in fullCV TMLE:", e$message))
    return(NULL)
  })
}

#' Run Full Cross-validated TMLE with Repetition
#'
#' Runs fullCV TMLE multiple times and aggregates the results to reduce
#' the variability in the estimates from the cross-validation procedure.
#'
#' @param Y Outcome variable vector (binary 0/1)
#' @param A Treatment variable vector (binary 0/1)
#' @param W Covariate matrix or data frame
#' @param Q.SL.library SuperLearner library for outcome model
#' @param g.SL.library SuperLearner library for treatment model
#' @param V Number of CV folds
#' @param family Distribution family ("binomial" or "gaussian")
#' @param num_repetitions Number of times to repeat fullCV-TMLE
#' @param max_attempts Maximum number of attempts for valid results
#' @param ... Additional arguments
#'
#' @return A list containing aggregated fullCV-TMLE results
#' @export
fullcv_tmle_multiple <- function(Y, A, W,
                                 Q.SL.library = c("SL.glm", "SL.mean"),
                                 g.SL.library = c("SL.glm", "SL.mean"),
                                 V = 5,
                                 family = "binomial",
                                 num_repetitions = 100,
                                 max_attempts = 500,
                                 ...) {

  # Input validation
  validate_inputs(Y, A, W, family)

  # Function to run a single fullCV-TMLE iteration
  run_iteration <- function(i) {
    # Set a different seed for each iteration
    set.seed(123 + i)

    result <- fullcv_tmle(
      Y = Y, A = A, W = W,
      Q.SL.library = Q.SL.library,
      g.SL.library = g.SL.library,
      V = V,
      family = family,
      ...
    )

    if (!is.null(result)) {
      return(result$results)
    }
    return(NULL)
  }

  # Initialize storage for valid results
  valid_results <- list()
  total_attempts <- 0

  # Keep trying until we get enough results or hit max attempts
  while (length(valid_results) < num_repetitions &&
         total_attempts < max_attempts) {

    # Run a single iteration
    result <- run_iteration(total_attempts + 1)
    total_attempts <- total_attempts + 1

    # Store valid result
    if (!is.null(result) && !any(is.na(result))) {
      valid_results[[length(valid_results) + 1]] <- result
    }

    # Print progress every 10 attempts
    if (total_attempts %% 10 == 0 || length(valid_results) == num_repetitions) {
      message(paste(length(valid_results),
                    "valid results after",
                    total_attempts, "attempts"))
    }
  }

  # Handle case where we didn't get enough valid results
  if (length(valid_results) < num_repetitions) {
    warning(paste("Only obtained", length(valid_results),
                  "valid results after", total_attempts, "attempts"))
    if (length(valid_results) == 0) {
      return(NULL)
    }
  }

  # Combine results
  runs <- dplyr::bind_rows(valid_results[1:min(length(valid_results), num_repetitions)])

  # Calculate medians
  medians <- apply(runs, 2, stats::median)

  # Calculate variance-adjusted estimates
  runs <- runs %>%
    dplyr::mutate(
      mv1 = v1 + (r1 - medians["r1"])^2,
      mv0 = v0 + (r0 - medians["r0"])^2,
      mvd = vd + (rd - medians["rd"])^2
    )

  # Final aggregated results
  aggregated_results <- apply(runs, 2, stats::median)

  return(list(
    repetitions = runs,
    aggregated_results = aggregated_results
  ))
}


#' Run Single Crossfit TMLE
#'
#' Implements Single Crossfit TMLE where sample is split into two parts
#' and models trained on one part are used to make predictions on the other.
#' This function is intended for internal use within singlecrossfit_tmle_multiple.
#'
#' @param Y Outcome variable vector (binary 0/1)
#' @param A Treatment variable vector (binary 0/1)
#' @param W Covariate matrix or data frame
#' @param Q.SL.library SuperLearner library for outcome model
#' @param g.SL.library SuperLearner library for treatment model
#' @param family Distribution family ("binomial" or "gaussian")
#' @param ... Additional arguments
#'
#' @return A tibble containing Single Crossfit TMLE results
singlecrossfit_tmle_internal <- function(Y, A, W,
                                         Q.SL.library = c("SL.glm", "SL.mean"),
                                         g.SL.library = c("SL.glm", "SL.mean"),
                                         family = "binomial",
                                         ...) {

  # Input validation
  validate_inputs(Y, A, W, family)

  # Ensure W is a data frame
  if (is.matrix(W)) {
    W <- as.data.frame(W)
  }

  # Create full data frame
  data <- data.frame(Y = Y, A = A, W)

  tryCatch({
    # Split sample into two parts
    n <- nrow(data)
    splits <- sample(rep(1:2, diff(floor(n * c(0, 1/2, 2/2)))))
    data$s <- splits

    # Create nested dataset
    dat_nested <- data %>%
      dplyr::group_by(s) %>%
      tidyr::nest()

    # Fit propensity score model
    pi_fitter <- function(df) {
      SuperLearner::SuperLearner(
        Y = as.matrix(df[, "A"]),
        X = df[, !(names(df) %in% c("Y", "A", "s")), drop = FALSE],
        family = binomial(),
        SL.library = g.SL.library,
        ...
      )
    }

    dat_nested$pi_fit <- purrr::map(dat_nested$data, pi_fitter)

    # Fit outcome model
    mu_fitter <- function(df) {
      SuperLearner::SuperLearner(
        Y = as.matrix(df[, "Y"]),
        X = df[, !(names(df) %in% c("Y", "s")), drop = FALSE],
        family = if(family == "binomial") binomial() else gaussian(),
        SL.library = Q.SL.library,
        ...
      )
    }

    dat_nested$mu_fit <- purrr::map(dat_nested$data, mu_fitter)

    # Calculate propensity scores and outcomes
    data <- data %>%
      dplyr::mutate(
        pi = case_when(
          s == 1 ~ predict(dat_nested$pi_fit[[2]], newdata = data[, !(names(data) %in% c("Y", "A", "s")), drop = FALSE])$pred,
          s == 2 ~ predict(dat_nested$pi_fit[[1]], newdata = data[, !(names(data) %in% c("Y", "A", "s")), drop = FALSE])$pred
        ),
        mu = case_when(
          s == 1 ~ predict(dat_nested$mu_fit[[2]], newdata = data[, !(names(data) %in% c("Y", "s")), drop = FALSE])$pred,
          s == 2 ~ predict(dat_nested$mu_fit[[1]], newdata = data[, !(names(data) %in% c("Y", "s")), drop = FALSE])$pred
        )
      )

    # Create datasets with A set to 0 and 1
    dat1 <- data
    dat0 <- data
    dat1$A <- 1
    dat0$A <- 0

    # Calculate mu1 and mu0
    data <- data %>%
      dplyr::mutate(
        mu1 = case_when(
          s == 1 ~ predict(dat_nested$mu_fit[[2]], newdata = dat1[, !(names(dat1) %in% c("Y", "s")), drop = FALSE])$pred,
          s == 2 ~ predict(dat_nested$mu_fit[[1]], newdata = dat1[, !(names(dat1) %in% c("Y", "s")), drop = FALSE])$pred
        ),
        mu0 = case_when(
          s == 1 ~ predict(dat_nested$mu_fit[[2]], newdata = dat0[, !(names(dat0) %in% c("Y", "s")), drop = FALSE])$pred,
          s == 2 ~ predict(dat_nested$mu_fit[[1]], newdata = dat0[, !(names(dat0) %in% c("Y", "s")), drop = FALSE])$pred
        )
      )

    # Truncate propensity scores and calculate clever covariates
    data <- data %>%
      dplyr::mutate(
        pi = ifelse(pi < 0.025, 0.025, ifelse(pi > 0.975, 0.975, pi)),
        H1 = A / pi,
        H0 = (1 - A) / (1 - pi)
      )

    # Adjust for zero probabilities
    data <- data %>%
      dplyr::mutate(
        mu = ifelse(mu == 0, 1e-17, ifelse(mu == 1, 1 - 1e-17, mu)),
        mu1 = ifelse(mu1 == 0, 1e-17, ifelse(mu1 == 1, 1 - 1e-17, mu1)),
        mu0 = ifelse(mu0 == 0, 1e-17, ifelse(mu0 == 1, 1 - 1e-17, mu0))
      )

    # Targeting step
    epsilon <- stats::coef(glm(Y ~ -1 + H0 + H1 + offset(qlogis(mu)), data = data, family = binomial))

    # Update initial estimates
    data <- data %>%
      dplyr::mutate(
        mu0_1 = plogis(qlogis(mu0) + epsilon[1] / (1 - pi)),
        mu1_1 = plogis(qlogis(mu1) + epsilon[2] / pi)
      )

    # Calculate final estimates
    r1 <- mean(data$mu1_1)
    r0 <- mean(data$mu0_1)
    rd <- r1 - r0

    # Calculate influence functions
    data <- data %>%
      dplyr::mutate(
        if1 = A / pi * (Y - mu1_1) + mu1_1 - r1,
        if0 = (1 - A) / (1 - pi) * (Y - mu0_1) + mu0_1 - r0,
        ifd = if1 - if0
      )

    # Calculate variances
    v1 <- stats::var(data$if1) / n
    v0 <- stats::var(data$if0) / n
    vd <- stats::var(data$ifd) / n

    # Compile results
    results <- dplyr::tibble(r1, r0, rd, v1, v0, vd)
    return(results)

  }, error = function(e) {
    warning(paste("Error in Single Crossfit TMLE:", e$message))
    return(NULL)
  })
}

#' Run Single Crossfit TMLE with Multiple Repetitions
#'
#' Implements Single Crossfit TMLE with multiple repetitions to reduce
#' the variability from the sample splitting procedure.
#'
#' @param Y Outcome variable vector (binary 0/1)
#' @param A Treatment variable vector (binary 0/1)
#' @param W Covariate matrix or data frame
#' @param Q.SL.library SuperLearner library for outcome model
#' @param g.SL.library SuperLearner library for treatment model
#' @param family Distribution family ("binomial" or "gaussian")
#' @param num_repetitions Number of times to repeat the procedure
#' @param max_attempts Maximum number of attempts for valid results
#' @param ... Additional arguments
#'
#' @return A list containing aggregated Single Crossfit TMLE results
#' @export
singlecrossfit_tmle <- function(Y, A, W,
                                Q.SL.library = c("SL.glm", "SL.mean"),
                                g.SL.library = c("SL.glm", "SL.mean"),
                                family = "binomial",
                                num_repetitions = 100,
                                max_attempts = 1000,
                                ...) {

  # Input validation
  validate_inputs(Y, A, W, family)

  # Initialize results
  runs <- dplyr::tibble(r1=double(), r0=double(), rd=double(), v1=double(), v0=double(), vd=double())
  valid_results <- 0
  total_attempts <- 0

  # Run until we get num_repetitions valid results or hit max attempts
  while(valid_results < num_repetitions && total_attempts < max_attempts) {
    total_attempts <- total_attempts + 1

    # Set a different seed for each attempt
    set.seed(123 + total_attempts)

    # Try to get a valid result
    tryCatch({
      result <- singlecrossfit_tmle_internal(
        Y = Y, A = A, W = W,
        Q.SL.library = Q.SL.library,
        g.SL.library = g.SL.library,
        family = family,
        ...
      )

      # Check if result contains any NA values
      if(!is.null(result) && !any(is.na(result))) {
        valid_results <- valid_results + 1

        # Append result to runs
        runs <- dplyr::bind_rows(runs, result)

        # Print progress periodically
        if(valid_results %% 10 == 0 || valid_results == num_repetitions) {
          message(paste("Valid result", valid_results, "of", num_repetitions, "obtained"))
        }
      }
    }, error = function(e) {
      if(total_attempts %% 50 == 0) {
        warning(paste("Error in singlecrossfit_tmle attempt", total_attempts, ":", e$message))
      }
    })
  }

  if(valid_results < num_repetitions) {
    warning(paste("Could only obtain", valid_results, "valid results out of", num_repetitions, "requested"))
    if(valid_results == 0) {
      return(NULL)
    }
  }

  # Medians of splits
  medians <- apply(runs, 2, stats::median)

  # Corrected variance terms
  runs <- runs %>%
    dplyr::mutate(
      mv1 = v1 + (r1-medians["r1"])^2,
      mv0 = v0 + (r0-medians["r0"])^2,
      mvd = vd + (rd-medians["rd"])^2
    )

  # Calculate final results
  aggregated_results <- apply(runs, 2, stats::median)

  return(list(
    repetitions = runs,
    aggregated_results = aggregated_results
  ))
}


#' Run Double Crossfit TMLE (Internal function)
#'
#' Implements Double Crossfit TMLE where sample is split into three parts,
#' with models trained on different parts and cross-fitted.
#' This function is intended for internal use within doublecrossfit_tmle.
#'
#' @param Y Outcome variable vector (binary 0/1)
#' @param A Treatment variable vector (binary 0/1)
#' @param W Covariate matrix or data frame
#' @param Q.SL.library SuperLearner library for outcome model
#' @param g.SL.library SuperLearner library for treatment model
#' @param family Distribution family ("binomial" or "gaussian")
#' @param ... Additional arguments
#'
#' @return A tibble containing Double Crossfit TMLE results
doublecrossfit_tmle_internal <- function(Y, A, W,
                                         Q.SL.library = c("SL.glm", "SL.mean"),
                                         g.SL.library = c("SL.glm", "SL.mean"),
                                         family = "binomial",
                                         ...) {

  # Input validation
  validate_inputs(Y, A, W, family)

  # Ensure W is a data frame
  if (is.matrix(W)) {
    W <- as.data.frame(W)
  }

  # Create full data frame
  data <- data.frame(Y = Y, A = A, W)

  tryCatch({
    # Split sample into three parts
    n <- nrow(data)
    splits <- sample(rep(1:3, diff(floor(n * c(0, 1/3, 2/3, 3/3)))))
    data$s <- splits

    # Create nested dataset
    dat_nested <- data %>%
      dplyr::group_by(s) %>%
      tidyr::nest()

    # Fit propensity score model
    pi_fitter <- function(df) {
      SuperLearner::SuperLearner(
        Y = as.matrix(df[, "A"]),
        X = df[, !(names(df) %in% c("Y", "A", "s")), drop = FALSE],
        family = binomial(),
        SL.library = g.SL.library,
        ...
      )
    }

    dat_nested$pi_fit <- purrr::map(dat_nested$data, pi_fitter)

    # Fit outcome model
    mu_fitter <- function(df) {
      SuperLearner::SuperLearner(
        Y = as.matrix(df[, "Y"]),
        X = df[, !(names(df) %in% c("Y", "s")), drop = FALSE],
        family = if(family == "binomial") binomial() else gaussian(),
        SL.library = Q.SL.library,
        ...
      )
    }

    dat_nested$mu_fit <- purrr::map(dat_nested$data, mu_fitter)

    # Calculate p-scores using each split
    data <- data %>%
      dplyr::mutate(
        pi1 = predict(dat_nested$pi_fit[[1]], newdata = data[, !(names(data) %in% c("Y", "A", "s")), drop = FALSE])$pred,
        pi2 = predict(dat_nested$pi_fit[[2]], newdata = data[, !(names(data) %in% c("Y", "A", "s")), drop = FALSE])$pred,
        pi3 = predict(dat_nested$pi_fit[[3]], newdata = data[, !(names(data) %in% c("Y", "A", "s")), drop = FALSE])$pred
      )

    # Calculate outcome model predictions
    data <- data %>%
      dplyr::mutate(
        mu_1 = predict(dat_nested$mu_fit[[1]], newdata = data[, !(names(data) %in% c("Y", "s")), drop = FALSE])$pred,
        mu_2 = predict(dat_nested$mu_fit[[2]], newdata = data[, !(names(data) %in% c("Y", "s")), drop = FALSE])$pred,
        mu_3 = predict(dat_nested$mu_fit[[3]], newdata = data[, !(names(data) %in% c("Y", "s")), drop = FALSE])$pred
      )

    # Create counterfactual datasets with A set to 0 and 1
    dat1 <- data
    dat0 <- data
    dat1$A <- 1
    dat0$A <- 0

    # Calculate counterfactual predictions
    data <- data %>%
      dplyr::mutate(
        mu1_1 = predict(dat_nested$mu_fit[[1]], newdata = dat1[, !(names(dat1) %in% c("Y", "s")), drop = FALSE])$pred,
        mu1_2 = predict(dat_nested$mu_fit[[2]], newdata = dat1[, !(names(dat1) %in% c("Y", "s")), drop = FALSE])$pred,
        mu1_3 = predict(dat_nested$mu_fit[[3]], newdata = dat1[, !(names(dat1) %in% c("Y", "s")), drop = FALSE])$pred,
        mu0_1 = predict(dat_nested$mu_fit[[1]], newdata = dat0[, !(names(dat0) %in% c("Y", "s")), drop = FALSE])$pred,
        mu0_2 = predict(dat_nested$mu_fit[[2]], newdata = dat0[, !(names(dat0) %in% c("Y", "s")), drop = FALSE])$pred,
        mu0_3 = predict(dat_nested$mu_fit[[3]], newdata = dat0[, !(names(dat0) %in% c("Y", "s")), drop = FALSE])$pred
      )

    # Truncate propensity scores
    data <- data %>%
      dplyr::mutate(
        pi1 = ifelse(pi1 < 0.025, 0.025, ifelse(pi1 > 0.975, 0.975, pi1)),
        pi2 = ifelse(pi2 < 0.025, 0.025, ifelse(pi2 > 0.975, 0.975, pi2)),
        pi3 = ifelse(pi3 < 0.025, 0.025, ifelse(pi3 > 0.975, 0.975, pi3))
      )

    # Calculate clever covariates
    data <- data %>%
      dplyr::mutate(
        H1_1 = A / pi1,
        H1_2 = A / pi2,
        H1_3 = A / pi3,
        H0_1 = (1 - A) / (1 - pi1),
        H0_2 = (1 - A) / (1 - pi2),
        H0_3 = (1 - A) / (1 - pi3)
      )

    # Adjust for zero probabilities
    data <- data %>%
      dplyr::mutate(
        mu_1 = ifelse(mu_1 == 0, 1e-17, ifelse(mu_1 == 1, 1 - 1e-17, mu_1)),
        mu_2 = ifelse(mu_2 == 0, 1e-17, ifelse(mu_2 == 1, 1 - 1e-17, mu_2)),
        mu_3 = ifelse(mu_3 == 0, 1e-17, ifelse(mu_3 == 1, 1 - 1e-17, mu_3)),
        mu1_1 = ifelse(mu1_1 == 0, 1e-17, ifelse(mu1_1 == 1, 1 - 1e-17, mu1_1)),
        mu1_2 = ifelse(mu1_2 == 0, 1e-17, ifelse(mu1_2 == 1, 1 - 1e-17, mu1_2)),
        mu1_3 = ifelse(mu1_3 == 0, 1e-17, ifelse(mu1_3 == 1, 1 - 1e-17, mu1_3)),
        mu0_1 = ifelse(mu0_1 == 0, 1e-17, ifelse(mu0_1 == 1, 1 - 1e-17, mu0_1)),
        mu0_2 = ifelse(mu0_2 == 0, 1e-17, ifelse(mu0_2 == 1, 1 - 1e-17, mu0_2)),
        mu0_3 = ifelse(mu0_3 == 0, 1e-17, ifelse(mu0_3 == 1, 1 - 1e-17, mu0_3))
      )

    # Targeting steps (one for each fold)
    epsilon_1 <- stats::coef(glm(Y ~ -1 + H0_2 + H1_2 + offset(qlogis(mu_3)),
                                 data = data %>% dplyr::filter(s == 1),
                                 family = binomial))

    epsilon_2 <- stats::coef(glm(Y ~ -1 + H0_3 + H1_3 + offset(qlogis(mu_1)),
                                 data = data %>% dplyr::filter(s == 2),
                                 family = binomial))

    epsilon_3 <- stats::coef(glm(Y ~ -1 + H0_1 + H1_1 + offset(qlogis(mu_2)),
                                 data = data %>% dplyr::filter(s == 3),
                                 family = binomial))

    # Update initial estimates
    data <- data %>%
      dplyr::mutate(
        mu0_1_1 = plogis(qlogis(mu0_3) + epsilon_1[1] / (1 - pi2)),
        mu0_1_2 = plogis(qlogis(mu0_1) + epsilon_2[1] / (1 - pi3)),
        mu0_1_3 = plogis(qlogis(mu0_2) + epsilon_3[1] / (1 - pi1)),
        mu1_1_1 = plogis(qlogis(mu1_3) + epsilon_1[2] / pi2),
        mu1_1_2 = plogis(qlogis(mu1_1) + epsilon_2[2] / pi3),
        mu1_1_3 = plogis(qlogis(mu1_2) + epsilon_3[2] / pi1)
      )

    # Calculate fold-specific estimates
    r1_1 = mean((data %>% dplyr::filter(s == 1))$mu1_1_1)
    r1_2 = mean((data %>% dplyr::filter(s == 2))$mu1_1_2)
    r1_3 = mean((data %>% dplyr::filter(s == 3))$mu1_1_3)
    r0_1 = mean((data %>% dplyr::filter(s == 1))$mu0_1_1)
    r0_2 = mean((data %>% dplyr::filter(s == 2))$mu0_1_2)
    r0_3 = mean((data %>% dplyr::filter(s == 3))$mu0_1_3)

    # Calculate influence functions
    data <- data %>%
      dplyr::mutate(
        if1_1 = A / pi2 * (Y - mu1_1_1) + mu1_1_1 - r1_1,
        if1_2 = A / pi3 * (Y - mu1_1_2) + mu1_1_2 - r1_2,
        if1_3 = A / pi1 * (Y - mu1_1_3) + mu1_1_3 - r1_3,
        if0_1 = (1 - A) / (1 - pi2) * (Y - mu0_1_1) + mu0_1_1 - r0_1,
        if0_2 = (1 - A) / (1 - pi3) * (Y - mu0_1_2) + mu0_1_2 - r0_2,
        if0_3 = (1 - A) / (1 - pi1) * (Y - mu0_1_3) + mu0_1_3 - r0_3,
        ifd_1 = if1_1 - if0_1,
        ifd_2 = if1_2 - if0_2,
        ifd_3 = if1_3 - if0_3
      )

    # Calculate final estimates by averaging across folds
    r1 <- (r1_1 + r1_2 + r1_3) / 3
    r0 <- (r0_1 + r0_2 + r0_3) / 3
    rd <- (r1 - r0)

    # Calculate variances
    v1 <- (stats::var((data %>% dplyr::filter(s == 1))$if1_1) +
             stats::var((data %>% dplyr::filter(s == 2))$if1_2) +
             stats::var((data %>% dplyr::filter(s == 3))$if1_3)) / (3 * n)

    v0 <- (stats::var((data %>% dplyr::filter(s == 1))$if0_1) +
             stats::var((data %>% dplyr::filter(s == 2))$if0_2) +
             stats::var((data %>% dplyr::filter(s == 3))$if0_3)) / (3 * n)

    vd <- (stats::var((data %>% dplyr::filter(s == 1))$ifd_1) +
             stats::var((data %>% dplyr::filter(s == 2))$ifd_2) +
             stats::var((data %>% dplyr::filter(s == 3))$ifd_3)) / (3 * n)

    # Compile results
    results <- dplyr::tibble(r1, r0, rd, v1, v0, vd)
    return(results)

  }, error = function(e) {
    warning(paste("Error in Double Crossfit TMLE:", e$message))
    return(NULL)
  })
}

#' Run Double Crossfit TMLE with Multiple Repetitions
#'
#' Implements Double Crossfit TMLE with multiple repetitions to reduce
#' the variability from the sample splitting procedure.
#'
#' @param Y Outcome variable vector (binary 0/1)
#' @param A Treatment variable vector (binary 0/1)
#' @param W Covariate matrix or data frame
#' @param Q.SL.library SuperLearner library for outcome model
#' @param g.SL.library SuperLearner library for treatment model
#' @param family Distribution family ("binomial" or "gaussian")
#' @param num_repetitions Number of times to repeat the procedure
#' @param max_attempts Maximum number of attempts for valid results
#' @param ... Additional arguments
#'
#' @return A list containing aggregated Double Crossfit TMLE results
#' @export
doublecrossfit_tmle <- function(Y, A, W,
                                Q.SL.library = c("SL.glm", "SL.mean"),
                                g.SL.library = c("SL.glm", "SL.mean"),
                                family = "binomial",
                                num_repetitions = 100,
                                max_attempts = 1000,
                                ...) {

  # Input validation
  validate_inputs(Y, A, W, family)

  # Initialize results
  runs <- dplyr::tibble(r1=double(), r0=double(), rd=double(), v1=double(), v0=double(), vd=double())
  valid_results <- 0
  total_attempts <- 0

  # Run until we get num_repetitions valid results or hit max attempts
  while(valid_results < num_repetitions && total_attempts < max_attempts) {
    total_attempts <- total_attempts + 1

    # Set a different seed for each attempt
    set.seed(123 + total_attempts)

    # Try to get a valid result
    tryCatch({
      result <- doublecrossfit_tmle_internal(
        Y = Y, A = A, W = W,
        Q.SL.library = Q.SL.library,
        g.SL.library = g.SL.library,
        family = family,
        ...
      )

      # Check if result contains any NA values
      if(!is.null(result) && !any(is.na(result))) {
        valid_results <- valid_results + 1

        # Append result to runs
        runs <- dplyr::bind_rows(runs, result)

        # Print progress periodically
        if(valid_results %% 10 == 0 || valid_results == num_repetitions) {
          message(paste("Valid result", valid_results, "of", num_repetitions, "obtained"))
        }
      }
    }, error = function(e) {
      if(total_attempts %% 50 == 0) {
        warning(paste("Error in doublecrossfit_tmle attempt", total_attempts, ":", e$message))
      }
    })
  }

  if(valid_results < num_repetitions) {
    warning(paste("Could only obtain", valid_results, "valid results out of", num_repetitions, "requested"))
    if(valid_results == 0) {
      return(NULL)
    }
  }

  # Medians of splits
  medians <- apply(runs, 2, stats::median)

  # Corrected variance terms
  runs <- runs %>%
    dplyr::mutate(
      mv1 = v1 + (r1-medians["r1"])^2,
      mv0 = v0 + (r0-medians["r0"])^2,
      mvd = vd + (rd-medians["rd"])^2
    )

  # Calculate final results
  aggregated_results <- apply(runs, 2, stats::median)

  return(list(
    repetitions = runs,
    aggregated_results = aggregated_results
  ))
}

# Helper Functions

#' Validate TMLE inputs
#'
#' @param Y Outcome variable vector
#' @param A Treatment variable vector
#' @param W Covariate matrix or data frame
#' @param family Distribution family
#'
#' @return Invisible NULL if validation passes, otherwise throws an error
validate_inputs <- function(Y, A, W, family) {
  # Check Y
  if (!is.numeric(Y)) {
    stop("Y must be a numeric vector")
  }

  # Check A
  if (!is.numeric(A)) {
    stop("A must be a numeric vector")
  }

  # For binary outcomes/treatments
  if (family == "binomial") {
    if (!all(Y %in% c(0, 1))) {
      stop("For binomial family, Y must contain only 0s and 1s")
    }
    if (!all(A %in% c(0, 1))) {
      stop("A must contain only 0s and 1s")
    }
  }

  # Check W
  if (!is.data.frame(W) && !is.matrix(W)) {
    stop("W must be a data frame or matrix")
  }

  # Check lengths
  if (length(Y) != length(A) || length(Y) != nrow(W)) {
    stop("Y, A, and W must have the same number of observations")
  }

  # Check family
  if (!family %in% c("binomial", "gaussian")) {
    stop("family must be either 'binomial' or 'gaussian'")
  }

  invisible(NULL)
}

#' Format TMLE results
#'
#' @param tmle_fit A fitted tmle object
#'
#' @return A tibble with formatted results
format_tmle_results <- function(tmle_fit) {
  # Extract results
  dplyr::tibble(
    r1 = tmle_fit$estimates$EY1$psi,
    r0 = tmle_fit$estimates$EY0$psi,
    rd = tmle_fit$estimates$ATE$psi,
    v1 = tmle_fit$estimates$EY1$var.psi,
    v0 = tmle_fit$estimates$EY0$var.psi,
    vd = tmle_fit$estimates$ATE$var.psi
  )
}
