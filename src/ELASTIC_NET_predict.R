# Function for performing elastic net regression
# train_pheno: vector of training set phenotypes
# train_geno: matrix of training set genotypes
# test_geno: matrix of test set genotypes
# alpha: mixing parameter for elastic net regularization
# lambda: regularization parameter for elastic net regression
elastic_net_func <- function(train_pheno, train_geno, test_geno, alpha = 0.5, lambda = 0.1) {
  
  # Center and scale the genotype matrix
  train_geno <- scale(train_geno, center = TRUE, scale = TRUE)
  test_geno <- scale(test_geno, center = TRUE, scale = TRUE)
  
  # Fit the elastic net model
  library(glmnet)
  en <- glmnet(train_geno, train_pheno, alpha = alpha, lambda = lambda, standardize = FALSE)
  
  # Choose optimal lambda value based on cross-validation
  cv_model <- cv.glmnet(train_geno, train_pheno, alpha = alpha)
  best_lambda <- cv_model$lambda.min
  
  # Predict phenotypes for training and test data
  train_pred <- predict(en, newx = train_geno, s = best_lambda)
  test_pred <- predict(en, newx = test_geno, s = best_lambda)
  
  # Return predicted values and model object
  return_value <- list("train_predicted" = train_pred, "test_predicted" = test_pred, "model" = en)
  return(return_value)
}


# sol_VL
elastic_net_sol_VL <- elastic_net_func(sol_VL_train_phenotype, sol_VL_train_marker, sol_VL_test_marker)

test_predicted <- elastic_net_sol_VL$test_predicted
train_predicted <- elastic_net_sol_VL$train_predicted

summary(test_predicted)
summary(sol_VL_test_phenotype)


cor(sol_VL_train_phenotype,train_predicted)
cor(sol_VL_test_phenotype, test_predicted)
