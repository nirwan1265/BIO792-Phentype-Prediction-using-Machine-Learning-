lasso_func <- function(train_pheno, train_geno, test_geno) {
  
  # Center and scale the genotype matrix
  train_geno <- scale(train_geno, center = TRUE, scale = TRUE)
  test_geno <- scale(test_geno, center = TRUE, scale = TRUE)
  
  # Perform k-fold cross-validation to find optimal lambda value
  cv_model <- cv.glmnet(train_geno, train_pheno, alpha = 1)
  
  # Find optimal lambda value that minimizes test MSE
  best_lambda <- cv_model$lambda.min
  
  # Fit LASSO model with best lambda
  best_model <- glmnet(train_geno, train_pheno, alpha = 1, lambda = best_lambda)
  
  # Predict phenotypes for test and train data
  train_predicted <- predict(best_model, newx = train_geno)
  val_predicted <- predict(best_model, newx = test_geno)
  
    # Return predicted values and model object
  return_value <- list("train_predicted" = train_predicted, "val_predicted"=val_predicted, "model" = best_model)
  return(return_value)
}

# sol_VL
lasso_sol_VL <- lasso_func(sol_VL_train_phenotype, sol_VL_train_marker, sol_VL_test_marker)

test_predicted <- lasso_sol_VL$test_predicted
train_predicted <- lasso_sol_VL$train_predicted

summary(test_predicted)
summary(sol_VL_test_phenotype)


cor(sol_VL_train_phenotype,train_predicted)
cor(sol_VL_test_phenotype, test_predicted)

sol_VL_train_test_lasso <- cbind(sol_VL_test_phenotype, test_predicted)
