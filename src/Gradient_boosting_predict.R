# Define function for gradient boosting
gbm_func <- function(train_pheno, train_geno, test_geno, n_trees = 1000, learning_rate = 0.01, max_depth = 3) {
  # convert train_geno and test_geno to data frames
  train_geno <- as.data.frame(train_geno)
  test_geno <- as.data.frame(test_geno)
  # Define the GBM model
  gbm_model <- gbm(
    formula = train_pheno ~ .,
    distribution = "gaussian",
    data = train_geno,
    n.trees = n_trees,
    interaction.depth = max_depth,
    shrinkage = learning_rate,
    verbose = FALSE
  )
  
  # Predict the phenotypes for training and test data
  train_pred <- predict(gbm_model, newdata = train_geno, n.trees = n_trees)
  test_pred <- predict(gbm_model, newdata = test_geno, n.trees = n_trees)
  
  # Return predicted values and model object
  return_value <- list("train_predicted" = train_pred, "test_predicted" = test_pred, "model" = gbm_model)
  return(return_value)
}


# sol_VL
#gbm_sol_VL <- gbm_func(sol_VL_train_phenotype, sol_VL_train_marker, sol_VL_test_marker)
gbm_sol_VL <- gbm_func(sol_VL_train_phenotype, sol_VL_train_marker, as.data.frame(sol_VL_test_marker))


test_predicted <- gbm_sol_VL$test_predicted
train_predicted <- gbm_sol_VL$train_predicted

summary(test_predicted)
summary(sol_VL_test_phenotype)


cor(sol_VL_train_phenotype,train_predicted)
cor(sol_VL_test_phenotype, test_predicted)

sol_VL_train_test_gbm <- cbind(sol_VL_test_phenotype, test_predicted)
