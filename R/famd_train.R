#' Factor Analysis of Mixed Data (FAMD) 'training' to attribute model weights based on the similiarty of a patient's covariates to the model development cohorts
#'
#' @param target_variable A secondary PK parameter (CMAX - maximal concentration, AUC - exposure or CMIN - trough concentration).
#' @param train A training dataset with ID to indicate a subject (set of covariates), covariates, observed target variable (for example CMAX), predicted target variable (for example CMAX_PRED), MODEL name used for the simulation.
#' @param continuous_cov List of continuous covariates.
#' @param categorical_cov List of categorical covariates.
#' @returns The FAMD model with the principal components, centroids and covariance matrices.
#' @export
#' @import dplyr
#' @import tidyr
#' @import ggplot2
#' @import FNN
#' @import factoextra
#' @import here
#' @import probably
#' @import RSpectra
#' @import stringr
#' @import purrr
#' @import FactoMineR
#' @examples
#' result<- famd_train(train = train, continuous_cov = c("WT", "CREAT", "AGE"), categorical_cov = c("BURN", "SEX"), target_variable = "CMIN")

famd_train <- function(train, continuous_cov, categorical_cov, target_variable) {

  target_column <- paste0(target_variable, "_IND") # Individual concentrations end with _IND

  # Define covariates
  covariates <- c(continuous_cov, categorical_cov)

  # Get the different models
  model_groups <- unique(train$MODEL_COHORT)

  # Select relevant columns
  train <- train |>
    dplyr::select(all_of(c("MODEL_COHORT", covariates)))

  train <- train |>
    mutate(across(all_of(categorical_cov), as.factor))

  # Split into model-specific datasets
  latent_data <- lapply(model_groups, function(model) {
    train |>
      dplyr::filter(MODEL_COHORT == model) |>
      dplyr::select(all_of(covariates)) |>
      distinct()
  })
  names(latent_data) <- model_groups

  # Train FAMD model on combined dataset
  combined_dataset <- bind_rows(latent_data)
  n_cov <- length(covariates)
  pca_result <- FAMD(combined_dataset, ncp = n_cov, graph = FALSE)

  # Get the first n number of principal components which explain 90 % of the variance
  explained_var <- pca_result$eig[, 2] / 100
  optimal_components <- which(cumsum(explained_var) >= 0.9)[1]

  # Apply FAMD transformation to training data
  latent_spaces <- lapply(latent_data, function(data) {
    predict(pca_result, newdata = data)$coord[, 1:optimal_components]
  })

  # Compute centroids and regularized covariance matrices
  centroids <- lapply(latent_spaces, colMeans)
  cov_matrices <- lapply(latent_spaces, function(latent) {
    cov_matrix <- cov(latent)
    cov_matrix + diag(rep(1e-6, ncol(cov_matrix)))  # Regularization
  })

  # Contribution of different variables to dimension1
  contrib_dim1 <- fviz_contrib(pca_result, "var", axes = 1) +
    theme(
      axis.text = element_text(size = 16),
      axis.title = element_text(size = 20),
      plot.title = element_text(size=24)
    ) +
    ylim(0,100)

  contrib_dim2 <- fviz_contrib(pca_result, "var", axes = 2) +
    theme(
      axis.text = element_text(size = 16),
      axis.title = element_text(size = 20),
      plot.title = element_text(size=24)
    ) +
    ylim(0,100)

  # Contribution of quantitative variables
  contrib_quant_var <- fviz_famd_var(pca_result, "quanti.var", col.var = "contrib",
                                     gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                     repel = TRUE,
                                     arrowsize = 2,
                                     labelsize = 8)  +
    theme(
      axis.text = element_text(size = 16),
      axis.title = element_text(size = 20),
      plot.title = element_text(size=24)
    )

  principal_components = pca_result$eig[1:optimal_components, ]

  return(list(pca_result = pca_result, optimal_components = optimal_components, target_variable = target_variable, categorical_cov = categorical_cov, continuous_cov = continuous_cov,
              centroids = centroids, cov_matrices = cov_matrices, contrib_dim1 = contrib_dim1, contrib_dim2 = contrib_dim2, contrib_quant_var = contrib_quant_var, principal_components = principal_components))
}
