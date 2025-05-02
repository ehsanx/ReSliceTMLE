test_that("vanilla_tmle runs without errors", {
  data(tmle_example_data)
  
  # Extract a single simulation
  test_data <- tmle_example_data[tmle_example_data$sim_id == 1, ]
  
  # Define variables
  Y <- test_data$Y
  A <- test_data$A
  W <- test_data[, c("X1", "X2", "X3", "X4")]
  
  # Run vanilla TMLE
  result <- vanilla_tmle(Y, A, W)
  
  # Check if result has expected components
  expect_true(!is.null(result))
  expect_true("fit" %in% names(result))
  expect_true("results" %in% names(result))
  
  # Check if estimates are reasonable
  expect_true(result$results$r1 >= 0 && result$results$r1 <= 1)
  expect_true(result$results$r0 >= 0 && result$results$r0 <= 1)
  expect_true(abs(result$results$rd) < 1)  # ATE should be between -1 and 1
})