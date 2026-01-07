#' Train classification tree informed model ensembling to calculate weights.
#'
#' @param target_variable A secondary PK parameter (CMAX - maximal concentration, AUC - exposure or CMIN - trough concentration).
#' @param data A dataset with ID to indicate a subject (set of covariates), covariates, observed target variable (for example CMAX), predicted target variable (for example CMAX_PRED), MODEL name used for the simulation.
#' @param continuous_cov List of continuous covariates.
#' @param categorical_cov List of categorical covariates.
#' @param seed Random seed for reproducibility.
#' @param conservative_pruning If FALSE (the default) take the highest CP with the lowest cross-validation error else take the highest CP whose cross-validation error is within 1 SD of the lowest cross-validation error.
#' @returns A table with the decision trees' leaves for different subgroups and models with the probability of correct predictions (prob_yes = model weight before standardization) and a list of rpart objects with the classification trees, plots of confusion matrices for each model.
#' @export
#' @import tidyr
#' @import rpart
#' @import rpart.plot
#' @import dplyr
#' @import purrr
#' @import stringr
#' @examples
#' leaf_results <- classification_tree_model_ensembling_train(data = data, target_variable = "CMAX", continuous_cov = c("WT", "CRCL"), categorical_cov = c("SEX", "BURN"))

classification_tree_model_ensembling_train <- function(data, target_variable, continuous_cov, categorical_cov, seed = 1991,conservative_pruning = FALSE) {

  set.seed(seed)

  target_column <- paste0(target_variable, "_IND") # Individual concentrations end with _IND

  # Define the predicted/observed ratio column
  data <- data |>
    group_by(ID) |>
    dplyr::mutate(ratio = get(paste0(target_variable, "_PRED")) /
             get(target_column)) |>
    ungroup() |>
    dplyr::mutate(CORRECT = if_else(ratio >= 0.8 & ratio <= 1.25, "YES", "NO"))
  data$CORRECT <- as.factor(data$CORRECT)

  # Convert categorical variables to factors
  for (cat_var in categorical_cov) {
    data[[cat_var]] <- as.factor(data[[cat_var]])
  }

  # Extract unique models
  models <- unique(data$MODEL)

  covariates <- c(continuous_cov, categorical_cov)

  formula <- as.formula(paste("CORRECT ~", paste(covariates, collapse = " + ")))

  leaf_results <- list()
  cfit_list <- NULL
  confusion_matrices <- list()
  confusion_matrix_plots <- list()

  # Loop over models
  for (model in models) {
    sub <- data |>
      dplyr::filter(MODEL == model) |>
      dplyr::select(all_of(continuous_cov), all_of(categorical_cov), CORRECT)

    # Convert target variable to factor
    sub$CORRECT <- as.factor(sub$CORRECT)

    # Cross-validation to find best cp
    cfit_big <- rpart(
      formula,
      data = sub,
      method = 'class',
      parms = list(split = "information"),
      control = rpart.control(
        cp = 0.001,
        xval = 10,
        minsplit = 4,
        minbucket = 4
      )
    )


    cp.select <- function(big.tree, conservative_pruning = FALSE) {
      cptable <- big.tree$cptable |> as_tibble()


      # Either take the highest CP with the lowest cross-validation error
      # Or be more conservative and take the highest CP whose cross-validation
      # error is within 1 SD of the lowest cross-validation error
      if(conservative_pruning == FALSE) {
        best_cp_value <- cptable |>
         dplyr::filter(xerror == min(xerror)) |>
          arrange(-CP) |> slice(1) |> pull(CP)
      } else {
        best_cp_value <- cptable |>
          mutate(lowest_xerror_cp = xerror == min(xerror)) |>
         dplyr::filter(xerror <= xerror[lowest_xerror_cp] + xstd[lowest_xerror_cp]) |>
          arrange(-CP) |> slice(1) |> pull(CP)
      }

      return(best_cp_value)
    }
    best_cp <- cp.select(cfit_big)
    cfit <- rpart::prune(cfit_big, cp = best_cp)

    # Plot tree
    plot <- rpart.plot(cfit, main = paste("Decision Tree for", model, "(cp =", best_cp, ")"))

    # Predictions for the confusion matrix
    preds <- predict(cfit, newdata = sub, type = "class")
    confusion_matrix <- table(Observed = sub$CORRECT, Predicted = preds)
    print(paste("Confusion matrix for model:", model))
    confusion_matrices[[model]] <- confusion_matrix
    confusion_matrix_df <- as.data.frame(confusion_matrix)
    # Plot confusion matrix
    confusion_matrix_plot <- ggplot(confusion_matrix_df, aes(x = Predicted, y = Observed, fill = Freq)) +
      geom_tile(color = "white") +
      geom_text(aes(label = Freq), vjust = 1) +
      scale_fill_gradient(low = "white", high = "steelblue") +
      labs(title = paste("Confusion Matrix -", model),
           x = "Predicted", y = "Observed") +
      theme_minimal()

    confusion_matrix_plots[[model]] <- confusion_matrix_plot # Save confusion matrix plot

    # Extract leaf data
    leaf_data <- cfit$frame |>
      tibble::rownames_to_column("node") |>
      dplyr::filter(var == "<leaf>") |>
      dplyr::mutate(
        prob_yes = yval2[, 5],
        n_obs = n
      )

    # Extract paths
    paths <- path.rpart(cfit, nodes = as.numeric(leaf_data$node), pretty = 0)

    leaf_table <- leaf_data |>
      dplyr::mutate(
        path = sapply(paths, function(path) paste(path, collapse = " -> ")),
        MODEL = model
      )

    # For the following part categorical covariates with more than two categories should be treated as continuous covariates
    multi_category_covs <- categorical_cov[sapply(sub[categorical_cov], function(x) length(unique(x)) > 2)]
    categorical_cov <- setdiff(categorical_cov, multi_category_covs)

    # Process continuous variables
    for (cont_var in continuous_cov) {
      leaf_table[[paste0(cont_var, "_lower")]] <- ifelse(
        str_detect(leaf_table$path, paste0(cont_var, ">=")),
        str_extract_all(leaf_table$path, paste0(cont_var, ">=\\s?\\d+\\.?\\d*")) |>
          lapply(function(x) max(as.numeric(str_extract(x, "\\d+\\.?\\d*")))) |> unlist(),
        NA
      )

      leaf_table[[paste0(cont_var, "_upper")]] <- ifelse(
        str_detect(leaf_table$path, paste0(cont_var, "<")),
        str_extract_all(leaf_table$path, paste0(cont_var, "<\\s?\\d+\\.?\\d*")) |>
          lapply(function(x) min(as.numeric(str_extract(x, "\\d+\\.?\\d*")))) |> unlist(),
        NA
      )
    }

    # Process categorical variables
    for (cat_var in categorical_cov) {
      leaf_table[[cat_var]] <- ifelse(
        str_detect(leaf_table$path, paste0(cat_var, "=1")), "1",
        ifelse(str_detect(leaf_table$path, paste0(cat_var, "=0")), "0", NA)
      )
    }

    for (multi_cat_var in multi_category_covs) {
      leaf_table[[paste0(multi_cat_var, "_lower")]] <- ifelse(
        str_detect(leaf_table$path, paste0(multi_cat_var, "=")),
        str_extract_all(leaf_table$path, paste0(multi_cat_var, "=\\d+(,\\d+)*")) |>
          lapply(function(x) min(as.numeric(unlist(strsplit(x, "[^0-9]+"))), na.rm = TRUE)) |>
          unlist(),
        NA
      )

      leaf_table[[paste0(multi_cat_var, "_upper")]] <- ifelse(
        str_detect(leaf_table$path, paste0(multi_cat_var, "=")),
        str_extract_all(leaf_table$path, paste0(multi_cat_var, "=\\d+(,\\d+)*")) |>
          lapply(function(x) max(as.numeric(unlist(strsplit(x, "[^0-9]+"))), na.rm = TRUE)) |>
          unlist(),
        NA
      )
    }

    # Bind results
    cfit_list[[length(cfit_list)+1]] <- cfit

    leaf_results <- bind_rows(leaf_results, leaf_table)
  }

  names(cfit_list) <- models

  # Convert numeric columns
  for (cont_var in continuous_cov) {
    leaf_results[[paste0(cont_var, "_lower")]] <- as.numeric(leaf_results[[paste0(cont_var, "_lower")]])
    leaf_results[[paste0(cont_var, "_upper")]] <- as.numeric(leaf_results[[paste0(cont_var, "_upper")]])
  }

  leaf_results <- leaf_results |>
    dplyr::select(-n, -var, -wt, -dev, -yval, -complexity, -starts_with("yval"))

  return(list(leaf_results = leaf_results, plot = plot, cfit_list = cfit_list, confusion_matrix = confusion_matrix_plots, target_variable = target_variable, continuous_cov = continuous_cov, categorical_cov = categorical_cov))
}

