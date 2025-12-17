library(xgboost)
library(dplyr)
library(tidyr)
library(caret)
library(parallel)
library(iml)

predict_SA_xgboost <- function(X, y, params,
                               run_dummy_quick_fit = FALSE,
                               set_params_man = TRUE,
                               bootstrap = TRUE,
                               n_bootstraps = 100) {
  
  r2_test_all_bootstraps <- c()
  n_iter <- ifelse(run_dummy_quick_fit, 5, 100)
  
  # Feature importance storage
  feature_importance_list <- list()
  
  for (b in 1:n_bootstraps) {
    
    cat(sprintf("\rBootstrap %d of %d", b, n_bootstraps))
    
    # Initialize predictions and counts
    train_predictions <- rep(0, length(y))
    test_predictions <- rep(0, length(y))
    train_counts <- rep(0, length(y))
    
    # 10-fold cross-validation
    folds <- createFolds(y, k = 10, list = TRUE, returnTrain = TRUE)
    
    for (i in seq_along(folds)) {
      train_index <- folds[[i]]
      test_index <- setdiff(seq_along(y), train_index)
      
      X_train <- X[train_index, ]
      y_train <- y[train_index]
      X_test <- X[test_index, ]
      y_test <- y[test_index]
      
      if (bootstrap) {
        boot_indices <- sample(seq_along(y_train), length(y_train), replace = TRUE)
        X_train_boot <- X_train[boot_indices, ]
        y_train_boot <- y_train[boot_indices]
      } else {
        X_train_boot <- X_train
        y_train_boot <- y_train
      }
      
      # Record locations of nans
      na_mask_train <- is.na(X_train_boot%>% select(-Site))
      na_mask_test  <- is.na(X_test%>% select(-Site))
      
      # Perform median imputation by site for every column of train set
      X_train_imputed <- impute_median_by_site(X_train_boot, site_col="Site")
      
      feature_cols <- setdiff(
        names(X_train_boot),
        c("Site", "Sex")
      )
      
      # Apply train median per site to test data
      train_site_medians <- X_train_boot %>%
        group_by(Site) %>%
        summarise(
          across(all_of(feature_cols), ~ median(.x, na.rm = TRUE)),
          .groups = "drop"
        )
      X_test_imputed <- X_test
      
      for (col in feature_cols) {
        X_test_imputed[[col]] <- ifelse(
          is.na(X_test_imputed[[col]]),
          train_site_medians[[col]][match(X_test_imputed$Site, train_site_medians$Site)],
          X_test_imputed[[col]]
        )
      }
      
      # Perform covbat harmonization
      
      batch_train <- as.factor(X_train_imputed$Site)
      batch_test <- factor(X_test_imputed$Site, levels = levels(batch_train))
      
      train_features <- X_train_imputed %>% select(-Site, -Sex)
      test_features <- X_test_imputed %>% select(-Site, -Sex)
      
      covbat_fit <- covfam_edited(
        data = train_features,
        bat = batch_train,
        robust.LS = FALSE,    # set explicitly
      )
      test_corrected <- predict(covbat_fit, test_features, batch_test)
      
      # Add sex back into final harmonized dataframe
      X_train_harmonized <- cbind(Sex = X_train_imputed$Sex, covbat_fit$dat.covbat)
      X_test_harmonized  <- cbind(Sex = X_test_imputed$Sex, test_corrected$dat.covbat)
      
      stopifnot(identical(colnames(X_train_harmonized),
                          colnames(na_mask_train)))
      
      stopifnot(identical(colnames(X_test_harmonized),
                          colnames(na_mask_test)))
      
      # Restore Nans in data
      X_train_harmonized[na_mask_train] <- NA
      X_test_harmonized[na_mask_test]   <- NA
      
      # Convert to DMatrix
      dtrain <- xgb.DMatrix(data = as.matrix(X_train_harmonized), label = y_train_boot)
      dtest <- xgb.DMatrix(data = as.matrix(X_test_harmonized), label = y_test)
      
      if (!set_params_man) {
        # Hyperparameter tuning (Bayesian optimization) can go here
        # Using ParBayesianOptimization or caret::train with tuneLength
      } else {
        xgb_model <- xgb.train(params = list(
          objective = "reg:squarederror",
          eta = params$eta,
          max_depth = params$max_depth,
          min_child_weight = params$min_child_weight,
          gamma = params$gamma,
          subsample = params$subsample,
          colsample_bytree = params$colsample_bytree
        ),
        data = dtrain,
        nrounds = params$n_estimators,
        verbose = 0)
        
        # Predictions
        test_predictions[test_index] <- predict(xgb_model, dtest)
        train_predictions[train_index] <- train_predictions[train_index] + predict(xgb_model, dtrain)
        
        # Store feature importance
        feature_importance_list[[length(feature_importance_list) + 1]] <- xgb.importance(model = xgb_model)
      }
      
      train_counts[train_index] <- train_counts[train_index] + 1
    }
    
    # Average train predictions across folds
    train_predictions <- train_predictions / train_counts
    
    r2_test <- caret::R2(test_predictions, y)
    r2_train <- caret::R2(train_predictions, y)
    
    print(paste0("Bootstrap ", b, " R2 test: ", round(r2_test, 3)))
    r2_test_all_bootstraps <- c(r2_test_all_bootstraps, r2_test)
  }
  
  # Aggregate feature importances
  # (Combine feature_importance_list into a summary dataframe)
  
  return(list(r2_test_array = r2_test_all_bootstraps,
              feature_importance_list = feature_importance_list))
}