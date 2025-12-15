library(xgboost)
library(dplyr)
library(tidyr)
library(neuroCombat)
library(caret)
library(parallel)
library(iml)

predict_SA_xgboost <- function(X, y, group_vals, sex_vals, target, metric, params,
                               run_dummy_quick_fit = FALSE,
                               set_params_man = TRUE,
                               bootstrap = TRUE,
                               n_bootstraps = 100) {
  
  r2_test_all_bootstraps <- c()
  n_iter <- ifelse(run_dummy_quick_fit, 5, 100)
  
  # Feature importance storage
  feature_importance_list <- list()
  
  for (b in 1:n_bootstraps) {
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
      
      # NeuroCombat harmonization
      batch_train <- as.factor(X_train_boot$Site)
      batch_test <- factor(X_test$Site, levels = levels(batch_train))
      
      X_train_harmonized <- neuroCombat(dat = t(as.matrix(X_train_boot %>% select(-Site, -Sex))),
                                        batch = batch_train)$dat.combat
      X_test_harmonized <- neuroCombat(dat = t(as.matrix(X_test %>% select(-Site, -Sex))),
                                       batch = batch_test,
                                       mod = NULL)$dat.combat
      
      # Add Sex back
      X_train_harmonized <- cbind(Sex = X_train_boot$Sex, t(X_train_harmonized))
      X_test_harmonized <- cbind(Sex = X_test$Sex, t(X_test_harmonized))
      
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