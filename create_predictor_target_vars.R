create_predictor_target_vars <- function(
    dforig,
    target,
    metric,
    include_group,
    run_dummy_quick_fit_xgb,
    show_heat_map,
    remove_colinear
) {
  
  library(dplyr)
  library(stringr)
  
  df <- dforig %>% as.data.frame()
  
  # Extract first three characters of Identifiers -> Site
  df <- df %>%
    mutate(Site = str_sub(Identifiers, 1, 3))
  
  # Move Site to be the third column
  df <- df %>%
    select(1:2, Site, everything(), -Identifiers)
  
  # Select brain predictor columns based on metric
  if (metric == "volume_infant") {
    
    pred_brain_cols <- names(df)[
      str_detect(names(df), "WM_V12|WM_V24|GM_V12|GM_V24")
    ]
    
  } else if (metric == "volume_VSA") {
    
    pred_brain_cols <- names(df)[
      str_detect(names(df), "WM_VSA|GM_VSA")
    ]
    
  } else if (metric == "subcort_infant") {
    
    pred_brain_cols <- names(df)[
      str_detect(
        names(df),
        paste(
          "Amygdala_v12|Putamen_v12|Caudate_v12|Thalamus_v12|",
          "GlobusPall_v12|Hippocampus_v12|",
          "Amygdala_v24|Putamen_v24|Caudate_v24|Thalamus_v24|",
          "GlobusPall_v24|Hippocampus_v24",
          sep = ""
        )
      )
    ]
    
  } else if (metric == "subcort_infant+volume_infant") {
    
    pred_brain_cols <- names(df)[
      str_detect(
        names(df),
        paste(
          "Amygdala_v12|Putamen_v12|Caudate_v12|Thalamus_v12|",
          "GlobusPall_v12|Hippocampus_v12|",
          "Amygdala_v24|Putamen_v24|Caudate_v24|Thalamus_v24|",
          "GlobusPall_v24|Hippocampus_v24|",
          "WM_V12|WM_V24|GM_V12|GM_V24",
          sep = ""
        )
      )
    ]
    
  } else if (metric == "subcort_VSA") {
    
    pred_brain_cols <- names(df)[
      str_detect(
        names(df),
        "(Amygdala|Caudate|Putamen|Thalamus|Hippocampus|Globus_Pall).*VSA"
      )
    ]
    
  } else if (metric == "cortical_thickness_VSA") {
    
    pred_brain_cols <- names(df)[
      str_detect(names(df), "CT_VSA")
    ]
    
  } else if (metric == "surface_area_VSA") {
    
    pred_brain_cols <- names(df)[
      str_detect(names(df), "SA_VSA")
    ]
    
  } else if (metric %in% c("fa_VSA", "rd_VSA", "md_VSA", "ad_VSA")) {
    
    prefix <- toupper(str_split(metric, "_")[[1]][1])
    pred_brain_cols <- names(df)[
      str_starts(names(df), paste0(prefix, "_"))
    ]
  }
  
  # Non-brain predictors
  if (include_group) {
    pred_non_brain_cols <- c(
      "Site", "Sex", "Group", "Group_HR+", "Group_HR-", "Group_LR-"
    )
  } else {
    pred_non_brain_cols <- c("Site", "Sex", "Group")
  }
  
  predictor_list <- c(pred_non_brain_cols, pred_brain_cols)
  
  # Keep only rows with non-missing target
  df <- df %>% filter(!is.na(.data[[target]]))
  
  if (run_dummy_quick_fit_xgb == 1) {
    set.seed(42)
    df <- df %>% sample_frac(0.1)
    n_iter <- 5
  }
  
  # Group and sex vectors
  group_vals <- df$Group %>% as.vector()
  sex_vals   <- df$Sex %>% as.vector()
  
  # Predictor matrix
  X <- df %>%
    select(all_of(predictor_list)) %>%
    select(-Group)
  
  # Target vector
  y <- df[[target]]
  
  # if (show_heat_map) {
  #   plot_title <- paste0("Correlation between regional ", metric)
  #   corr_matrix <- plot_correlations(X, plot_title)
  # }
  # 
  # if (remove_colinear) {
  #   X <- remove_collinearity(X, threshold = 0.9)
  #   plot_title <- "After removing colinear features"
  #   corr_matrix <- plot_correlations(X, plot_title)
  # }
  
  return(list(
    X = X,
    y = y,
    group_vals = group_vals,
    sex_vals = sex_vals
  ))
}