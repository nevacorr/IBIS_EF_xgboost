library(xgboost)
library(dplyr)
library(ggplot2)
library(neuroCombat)

# --- User-defined functions ---

source("load_all_data.R") 
source("create_predictor_target_vars.R")
source("predict_SA_xgboost.R") 
source("UtilityFunctions.R")
source("load_brain_data.R")

debug(load_all_data)
debug(create_predictor_target_vars)
debug(predict_SA_xgboost)
debug(reshape_dataframe)
debug(load_infant_subcortical_data)
debug(load_vsa_subcortical_data)
debug(load_vsa_ct_sa_data)
debug(load_and_clean_vsa_dti_data)
# debug(load_and_clean_infant_volume_data_and_all_behavior)
debug(load_and_clean_vsa_volume_data)
debug(plot_correlations)
debug(remove_collinearity)
debug(write_modeling_data_and_outcome_to_file)
debug(plot_xgb_actual_vs_pred)
debug(generate_bootstrap_indices)
debug(calculate_percentile)
debug(plot_r2_distribution)
debug(aggregate_feature_importances)


# ---- Parameters ----
target <- "BRIEF2_GEC_T_score"
metric <- "subcort_infant"
#options 'volume_infant', 'volume_VSA', 'subcort_VSA', 'subcort_infant', 'ad_VSA', 'rd_VSA', 'md_VSA', 'fa_VSA'
#        'surface_area_VSA', 'cortical_thickness_VSA', 'subcort_infant+volume_infant'
include_group <- 0
bootstrap <- TRUE
n_bootstraps <- 100
show_heat_map <- FALSE
remove_colinear <- FALSE
run_dummy_quick_fit_xgb <- FALSE
alpha <- 0.05

run_ridge_regression_fit <- FALSE
run_xgboost_fit <- TRUE
set_xgb_params_man <- TRUE

if (set_xgb_params_man) {
  params <- list(
    n_estimators = 50,
    min_child_weight = 1,
    gamma = 0.01,
    eta = 0.01542,
    subsample = 0.576,
    colsample_bytree = 0.2,
    max_depth = 6
  )
} else {
  params <- list(
    n_estimators = c(50, 2001),
    min_child_weight = c(1, 11),
    gamma = c(0.01, 5.0),  # Log-uniform tuning handled separately
    eta = c(0.005, 0.5),
    subsample = c(0.2, 1.0),
    colsample_bytree = c(0.2, 1.0),
    max_depth = c(2, 6)
  )
}

# ---- Load Data ----
df <- load_all_data()

data_list <- create_predictor_target_vars(
  df, target, metric, include_group,
  run_dummy_quick_fit_xgb, show_heat_map, remove_colinear
)

X <- data_list$X
y <- data_list$y
group_vals <- data_list$group_vals
sex_vals <- data_list$sex_vals

cat(sprintf("Running with target = %s, metric = %s, include_group = %d, quick fit = %d\n",
            target, metric, include_group, as.integer(run_dummy_quick_fit_xgb)))

# ---- Run XGBoost ----
if (run_xgboost_fit) {
  result <- predict_SA_xgboost_covbat(
    X, y, group_vals, sex_vals, target, metric, params,
    run_dummy_quick_fit_xgb, set_xgb_params_man, 
    show_results_plot = FALSE,
    bootstrap, n_bootstraps
  )
  
  r2_test_array_xgb <- result$r2_test_array
  feature_importance_df <- result$feature_importance_list
  
  # ---- Calculate Percentile ----
  percentile_result <- calculate_percentile(r2_test_array_xgb, alpha)
  result_text_xgb <- percentile_result$result_text
  percentile_value_xgb <- percentile_result$percentile_value
  
  # ---- Plot R2 Distribution ----
  plot_r2_distribution(r2_test_array_xgb, result_text_xgb, percentile_value_xgb,
                       target, metric, alpha, n_bootstraps, alg = "XGBoost")
}
