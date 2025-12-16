library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(gridExtra)

plot_correlations <- function(df, title) {
  
  df_corr <- df %>%
    dplyr::select(-any_of(c("Site", "Sex", "Group")))  # drop if present
  
  corr_matrix <- cor(df_corr, use = "pairwise.complete.obs")
  
  corr_df <- as.data.frame(as.table(corr_matrix))
  colnames(corr_df) <- c("Var1", "Var2", "Correlation")
  
  print(
    ggplot(corr_df, aes(Var1, Var2, fill = Correlation)) +
    geom_tile() +
    scale_fill_gradient2(
      low = "blue", mid = "white", high = "red",
      midpoint = 0, limits = c(-1, 1)
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, size = 8),
      axis.text.y = element_text(size = 8),
      plot.title = element_text(hjust = 0.5)
    ) +
    labs(title = title, x = NULL, y = NULL)
  )
}

remove_collinearity <- function(df, threshold = 0.9) {
  
  df_corr <- df %>%
    dplyr::select(-any_of(c("Site", "Sex", "Group")))  # drop if present
  
  corr_matrix <- abs(cor(df_corr, use = "pairwise.complete.obs"))
  
  upper_triangle <- corr_matrix
  upper_triangle[lower.tri(upper_triangle, diag = TRUE)] <- NA
  
  to_drop <- colnames(upper_triangle)[
    apply(upper_triangle, 2, function(x) any(x > threshold, na.rm = TRUE))
  ]
  
  df_reduced <- df %>% select(-all_of(to_drop))
  
  message("Dropped columns due to collinearity: ", paste(to_drop, collapse=", "))
  
  return(df_reduced)
}

write_modeling_data_and_outcome_to_file <- function(
    quick_run, metric, params, set_parameters_manually, target,
    X, r2_train, r2_test, best_params, bootstrap, elapsed_time
) {
  
  filename <- paste0(target, "_", metric, "_xgboost_run_results_summary.txt")
  
  con <- file(filename, open = "a")
  on.exit(close(con))
  
  now <- Sys.time()
  
  writeLines(paste("Run Date and Time:", format(now, "%Y-%m-%d %H:%M:%S")), con)
  writeLines("####### Model performance summary ######", con)
  
  if (quick_run == 1)
    writeLines("This was a quick run to check code, not fit model", con)
  
  if (set_parameters_manually == 0)
    writeLines("Used hyperparameter optimization", con)
  else
    writeLines("Set parameter manually", con)
  
  if (bootstrap == 1)
    writeLines("Bootstrapped", con)
  
  writeLines("Used harmonization by site", con)
  writeLines(paste("Metric:", metric), con)
  writeLines(paste("Target:", target), con)
  
  feature_names <- colnames(X)
  writeLines(paste("Features:", paste(feature_names, collapse = ", ")), con)
  
  writeLines("Parameter specified", con)
  for (nm in names(params)) {
    writeLines(paste(nm, ":", params[[nm]]), con)
  }
  
  writeLines(paste("Best Parameters:", paste(best_params, collapse = ", ")), con)
  writeLines("Performance metrics:", con)
  writeLines(sprintf("R2 train = %.4f", r2_train), con)
  writeLines(sprintf("R2 test = %.4f", r2_test), con)
  writeLines(sprintf("Run completion time: %.2f", elapsed_time), con)
}

plot_xgb_actual_vs_pred <- function(
    metric, target, r2_train, r2_test, df, best_params, show_plot
) {
  
  plot_df <- bind_rows(
    data.frame(
      actual = df[[target]],
      predicted = df$test_predictions,
      type = "Test",
      r2 = r2_test
    ),
    data.frame(
      actual = df[[target]],
      predicted = df$train_predictions,
      type = "Train",
      r2 = r2_train
    )
  )
  
  p <- ggplot(plot_df, aes(actual, predicted)) +
    geom_point(alpha = 0.6) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    facet_wrap(~ type, ncol = 1) +
    theme_minimal() +
    labs(
      x = paste("Actual", metric, target),
      y = paste("Predicted", metric, target),
      title = paste(
        metric, "Predictions\n",
        sprintf("colsample_bytree=%.2f n_estimators=%d min_child_weight=%d gamma=%.2f\neta=%.2e subsample=%.2f max_depth=%d harmonized by Site",
                best_params$colsample_bytree,
                best_params$n_estimators,
                best_params$min_child_weight,
                best_params$gamma,
                best_params$eta,
                best_params$subsample,
                best_params$max_depth)
      )
    )
  
  ggsave(
    filename = paste0(target, "_", metric, "_xgboost_actual_vs_predicted.png"),
    plot = p, width = 10, height = 12
  )
  
  if (show_plot == 1) print(p)
}

generate_bootstrap_indices <- function(n_rows, n_iterations = 100) {
  
  set.seed(42)
  
  lapply(seq_len(n_iterations), function(i) {
    sample(seq_len(n_rows), size = n_rows, replace = TRUE)
  })
}

calculate_percentile <- function(r2_test_array, alpha) {
  
  lower_bound <- quantile(r2_test_array, probs = alpha)
  
  if (lower_bound > 0) {
    result_text <- sprintf(
      "R² is significantly greater than 0 at the %d%% level (one-tailed)",
      alpha * 100
    )
  } else {
    result_text <- "R² is NOT significant at the 5% level (one-tailed)"
  }
  
  cat(sprintf("%d%% lower confidence bound: %.4f\n",
              (1 - alpha) * 100, lower_bound))
  
  list(
    result_text = result_text,
    percentile_value = lower_bound
  )
}

plot_r2_distribution <- function(
    r2_test_array, result_text, percentile_value,
    target, metric, alpha, n_bootstraps, alg
) {
  
  df <- data.frame(r2 = r2_test_array)
  
  ggplot(df, aes(r2)) +
    geom_histogram(bins = 30, fill = "skyblue", color = "black") +
    geom_vline(
      xintercept = percentile_value,
      color = "red", linetype = "dashed", linewidth = 1
    ) +
    theme_minimal() +
    labs(
      title = paste0(
        target, " predicted from ", metric, " with ", alg,
        "\nBootstrap r² test Distribution\n",
        alpha * 100, "% Percentile Marked\n",
        result_text, "\n",
        "nbootstraps=", n_bootstraps
      ),
      x = "r² test",
      y = "Frequency"
    )
}

aggregate_feature_importances <- function(
    feature_importance_list,
    feature_names,
    n_boot,
    outputfilename,
    top_n = 10,
    plot = TRUE
) {
  
  importance_matrix <- do.call(rbind, feature_importance_list)
  
  importance_df <- data.frame(
    feature = feature_names,
    mean_importance = colMeans(importance_matrix),
    std_importance = apply(importance_matrix, 2, sd)
  ) %>%
    arrange(desc(mean_importance))
  
  write_csv(importance_df, outputfilename)
  
  if (plot) {
    top_df <- importance_df %>% head(top_n)
    
    ggplot(top_df, aes(reorder(feature, mean_importance), mean_importance)) +
      geom_col(fill = "skyblue") +
      geom_errorbar(
        aes(ymin = mean_importance - std_importance,
            ymax = mean_importance + std_importance),
        width = 0.2
      ) +
      coord_flip() +
      theme_minimal() +
      labs(
        title = paste0(
          "Top ", top_n,
          " Feature Importances (Bootstrapped XGBoost ",
          n_boot, " bootstraps)"
        ),
        x = NULL,
        y = "Mean Importance"
      )
  }
  
  return(importance_df)
}

divide_columns_by_tottiss <- function(df_all_brain_behav, df, mystr) {
  
  suffix <- paste0("_", mystr)
  
  if (mystr == "VSA") {
    tot_col <- paste0("total_Tissue_vol_", toupper(mystr))
    icv_col <- paste0("ICV_vol_", toupper(mystr))
  } else {
    tot_col <- paste0("totTiss_", toupper(mystr))
    icv_col <- paste0("ICV_", toupper(mystr))
  }
  
  # Safety check
  if (!tot_col %in% names(df_all_brain_behav)) {
    stop(
      sprintf(
        "Cannot divide by totTiss. Column '%s' not found in df_all_brain_behav.",
        tot_col
      ),
      call. = FALSE
    )
  }
  
  # Find shared columns ending in _mystr (excluding CandID)
  shared_cols <- intersect(names(df), names(df_all_brain_behav))
  
  shared_cols <- shared_cols[
    shared_cols != "CandID" &
      endsWith(shared_cols, suffix)
  ]
  
  # Remove totTiss and ICV columns
  shared_cols <- setdiff(shared_cols, c(tot_col, icv_col))
  
  # Perform division
  df_all_brain_behav <- df_all_brain_behav %>%
    mutate(across(
      all_of(shared_cols),
      ~ .x / .data[[tot_col]]
    ))
  
  return(df_all_brain_behav)
}
