# ReSliceTMLE: Resampling-Based Targeted Maximum Likelihood Estimation

This R package implements robust variants of Targeted Maximum Likelihood Estimation (TMLE) using repeated cross-validation and sample splitting approaches to improve variance estimation in causal inference settings.

---

## 📋 Overview

ReSliceTMLE provides a unified framework for implementing various TMLE variants that leverage cross-validation and data-splitting techniques to obtain more robust causal estimates. These methods are particularly useful when:

- Standard TMLE produces unstable variance estimates
- Coverage of confidence intervals is suboptimal
- Dealing with complex confounding scenarios
- Working with high-dimensional data or complex model specifications

---

## 🚀 TMLE Variants

The package implements multiple TMLE variants:

| **Method**              | **Description**                                                           |
|-------------------------|---------------------------------------------------------------------------|
| `vanilla_tmle`          | Standard TMLE implementation                                              |
| `cvq_tmle`              | Cross-validated TMLE (Q model cross-validated)                            |
| `cvq_tmle_multiple`     | CV-TMLE with multiple repetitions for more stable inference               |
| `fullcv_tmle`           | Full cross-validated TMLE with both Q and g models cross-validated        |
| `fullcv_tmle_multiple`  | Full CV-TMLE with multiple repetitions                                    |
| `singlecrossfit_tmle`   | Single crossfit TMLE (sample-splitting with two-fold cross-fitting)       |
| `doublecrossfit_tmle`   | Double crossfit TMLE (enhanced sample-splitting with three-fold cross-fitting) |

---

## 📦 Installation

You can install ReSliceTMLE from GitHub:

```r
# install.packages("devtools")
devtools::install_github("ehsanx/ReSliceTMLE")
```

---

## 🔍 Usage Examples

### Basic TMLE Analysis

```r
library(ReSliceTMLE)

# Load example data
data("tmle_example_data")

# Extract a single simulation
test_data <- tmle_example_data[tmle_example_data$sim_id == 1, ]

# Run vanilla TMLE
result <- run_tmle_analysis(
  data = test_data,
  outcome_var = "Y",
  treatment_var = "A",
  covariates = c("X1", "X2", "X3", "X4"),
  tmle_variant = "vanilla"
)

# View formatted results
result$interpreted_results
```

### Comparing Multiple TMLE Variants

```r
# Compare different TMLE variants on the same dataset
multi_result <- run_multiple_tmle_variants(
  data = test_data,
  outcome_var = "Y",
  treatment_var = "A",
  covariates = c("X1", "X2", "X3", "X4"),
  tmle_variants = c("vanilla", "cvq", "doublecrossfit"),
  num_repetitions = 50
)

# View comparative results
multi_result$interpreted_results
```

### Simulation Study

```r
# Run TMLE on all simulations
sim_result <- run_tmle_analysis(
  data = tmle_example_data,
  outcome_var = "Y",
  treatment_var = "A",
  covariates = c("X1", "X2", "X3", "X4"),
  tmle_variant = "doublecrossfit",
  sim_id_col = "sim_id",
  true_value = 0.25  # True ATE for evaluation
)

# View simulation summary
sim_result$summary
```

---

## 🧰 Key Functions

| **Function**              | **Purpose**                                                              |
|---------------------------|--------------------------------------------------------------------------|
| `run_tmle_analysis`       | Main user interface to apply any TMLE variant to a dataset               |
| `run_multiple_tmle_variants` | Compare multiple TMLE variants on the same dataset                    |
| `vanilla_tmle`            | Implementation of standard TMLE                                          |
| `cvq_tmle`                | Implementation of CV-TMLE                                                |
| `fullcv_tmle`             | Implementation of full CV-TMLE                                           |
| `singlecrossfit_tmle`     | Implementation of single crossfit TMLE                                   |
| `doublecrossfit_tmle`     | Implementation of double crossfit TMLE                                   |
| `prepare_tmle_data`       | Process data for TMLE analysis (missing data handling, etc.)             |
| `format_results`          | Format TMLE results with confidence intervals and p-values               |
| `create_sl_library`       | Create SuperLearner library for fitting Q and g models                   |

---

## 🔧 Dependencies

The package relies on several R packages:

- **Core functionality**: dplyr, tidyr, purrr
- **Statistical methods**: tmle, SuperLearner
- **Utilities**: stats, utils, rlang

---

## 📊 Model Specification

ReSliceTMLE supports flexible model specification through SuperLearner:

```r
# Create a custom SuperLearner library
custom_library <- create_sl_library(
  include_glm = TRUE,
  include_gam = TRUE,
  include_rf = TRUE,
  include_gbm = FALSE
)

# Run TMLE with custom models
result <- run_tmle_analysis(
  data = test_data,
  outcome_var = "Y",
  treatment_var = "A",
  covariates = c("X1", "X2", "X3", "X4"),
  tmle_variant = "doublecrossfit",
  Q.SL.library = custom_library,
  g.SL.library = custom_library
)
```

---

## 📄 License

TBD

---

## 📚 Citation

TBD

```
[Citation information placeholder]
```