% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{prepare_tmle_data}
\alias{prepare_tmle_data}
\title{Prepare data for TMLE analysis}
\usage{
prepare_tmle_data(
  data,
  outcome_var,
  treatment_var,
  covariates,
  binary_outcome = TRUE,
  missing_handling = "complete_case"
)
}
\arguments{
\item{data}{A data frame containing the analysis variables}

\item{outcome_var}{Name of the outcome variable}

\item{treatment_var}{Name of the treatment variable}

\item{covariates}{Names of covariates to include in the analysis}

\item{binary_outcome}{Logical indicating if outcome is binary}

\item{missing_handling}{Strategy for handling missing values ("complete_case", "impute_mean", or "impute_model")}
}
\value{
A list containing prepared Y, A, and W components
}
\description{
This function prepares a dataset for TMLE analysis by handling missing values,
converting variable types, and extracting outcome, treatment, and covariates.
}
